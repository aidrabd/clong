import os
import sys
import numpy as np
import csv
from keras.models import load_model
from keras.utils import to_categorical

MAX_SEQ_LENGTH = 500
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
MODEL_PATH = "clong.h5"
OUTPUT_DIR = "Output"
MECHANISM_GENE_INFO = {
    "hns": ("Transcription & Translation", "global DNA-binding transcriptional dual regulator", "40%"),
    "ihfB": ("Transcription & Translation", "integration host factor; DNA-binding protein", "35%"),
    "hyfR": ("Transcription & Translation", "DNA-binding transcriptional activator", "16%"),
    "rplY": ("Transcription & Translation", "50S ribosomal subunit protein", "11%"),
    "aroG": ("Metabolism & Respiration", "3-deoxy-D-arabino-heptulosonate-7-phosphate synthase", "29%"),
    "aroD": ("Metabolism & Respiration", "3-dehydroquinate dehydratase", "24%"),
    "lipB": ("Metabolism & Respiration", "lipoyl-protein ligase", "23%"),
    "purE": ("Metabolism & Respiration", "N5-carboxyaminoimidazole ribonucleotide mutase", "21%"),
    "pdxA": ("Metabolism & Respiration", "4-hydroxy-L-threonine phosphate dehydrogenase", "21%"),
    "gmhA*": ("Metabolism & Respiration", "D-sedoheptulose 7-phosphate isomerase", "19%"),
    "pabB": ("Metabolism & Respiration", "aminodeoxychorismate synthase", "18%"),
    "ynjE": ("Metabolism & Respiration", "thiosulfate sulfur transferase", "18%"),
    "nrfG": ("Metabolism & Respiration", "heme lyase", "17%"),
    "psuK": ("Metabolism & Respiration", "pseudouridine kinase", "10%"),
    "lpp": ("Membrane & Transport", "murein lipoprotein", "27%"),
    "yfiB": ("Membrane & Transport", "outer membrane lipoprotein", "20%"),
    "sapD": ("Membrane & Transport", "antimicrobial peptide transporter", "17%"),
    "uidC": ("Membrane & Transport", "outer membrane porin protein", "17%"),
    "ygiV": ("Membrane & Transport", "inner membrane protein", "12%"),
    "secB": ("Protease & Chaperone", "protein export chaperone", "29%"),
    "lon": ("Protease & Chaperone", "DNA-binding ATP-dependent protease", "25%"),
    "skp": ("Protease & Chaperone", "periplasmic chaperone", "16%"),
    "pbl": ("Others", "lytic transglycosylase", "21%"),
    "ycbJ": ("Others", "unknown", "19%"),
    "trxA": ("Others", "thioredoxin", "15%"),
    "ycgL": ("Others", "unknown", "14%"),
    "ycgN": ("Others", "unknown", "12%"),
    "recC": ("Others", "exonuclease V", "11%"),
    "yfdG": ("Others", "prophage; bactoprenol-linked glucose translocase", "10%"),
}

def one_hot_encode(seq, max_length=MAX_SEQ_LENGTH):
    aa_to_index = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}
    one_hot = np.zeros((max_length, len(AMINO_ACIDS)), dtype=int)
    for i, aa in enumerate(seq[:max_length]):
        if aa in aa_to_index:
            one_hot[i, aa_to_index[aa]] = 1
    return one_hot

def load_fasta_sequences(filename):
    sequences = []
    with open(filename, 'r') as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                seq = ''
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences

def predict_sequences(model, sequences):
    X = np.array([one_hot_encode(seq) for seq in sequences])
    preds = model.predict(X)
    pred_labels = np.argmax(preds, axis=1)
    return preds, pred_labels

def main(fasta_file):
    if not os.path.exists(MODEL_PATH):
        print(f"Error: Model file '{MODEL_PATH}' not found in current directory.")
        sys.exit(1)

    if not os.path.exists(fasta_file):
        print(f"Error: Input FASTA file '{fasta_file}' not found.")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Loading model from {MODEL_PATH}...")
    model = load_model(MODEL_PATH)

    print(f"Loading sequences from {fasta_file}...")
    sequences = load_fasta_sequences(fasta_file)
    if len(sequences) == 0:
        print("No sequences found in input file.")
        sys.exit(1)
    print(f"Loaded {len(sequences)} sequences.")

    preds_prob, preds_label_indices = predict_sequences(model, sequences)

    labels = sorted(MECHANISM_GENE_INFO.keys())
    index_to_label = {i: label for i, label in enumerate(labels)}

    output_file = os.path.join(OUTPUT_DIR, "clong_predictions.csv")
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sequence', 'Predicted Gene', 'Predicted Probability', 'Mechanism', 'Protein', 'Lifespan Extension'])
        for seq, pred_idx, prob_vector in zip(sequences, preds_label_indices, preds_prob):
            pred_gene = index_to_label.get(pred_idx, 'Unknown')
            prob = float(prob_vector[pred_idx])
            info = MECHANISM_GENE_INFO.get(pred_gene, ('', '', ''))
            writer.writerow([seq, pred_gene, prob, *info])

    print(f"Prediction complete. Results saved to '{output_file}'")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clong_predict.py <input_fasta_file>")
        sys.exit(1)
    fasta_path = sys.argv[1]

    main(fasta_path)

