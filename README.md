# CLONG (C. elegans Longevity-related Gene Predictor)
CLONG is a model for predicting bacterial genes that influence _ C. elegans_ lifespan when used as food sources, particularly focusing on probiotic strains.

# Overview
CLONG analyzes bacterial proteome to identify genes that may affect _C. elegans_ longevity. The deletion of those genes could extend the lifespan of _C. elegans_. It was validated using _Lacticaseibacillus rhamnosus_ as a model probiotic strain, known for extending _ C. elegans_ lifespan and providing resistance to oxidative stress and pathogens.

# Features
- Predicts genes that potentially influence _ C. elegans_ lifespan
- Assigns confidence scores (predicted probability) to gene predictions
- Classifies genes into functional categories
- Identifies candidates for gene deletion to enhance probiotic effects

## Installation

### Prerequisites

- Ubuntu/Linux terminal

### Installation

```bash
# Clone and setup
git clone https://github.com/aidrabd/aidimin.git
cd aidimin

# Make prediction script executable
chmod +x predict.py
```

First, make sure you have conda installed:

```bash
1. Install  Miniconda (if not installed)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

2. Activate conda base environment

conda init

Then, restart your terminal or run:

source ~/.bashrc

After that, activate the base environment with:

conda activate
```

Second, make sure you have Python specific version installed:

```bash
conda create -n py312 python=3.12.9
conda activate py312
python --version
```

Third, make sure you have specific Tensorflow, Keras, numpy, scikit-learn versions installed:

```bash

tensorflow>=2.8.0
keras>=2.8.0
numpy>=1.21.0
scikit-learn>=1.0.0

conda install -c conda-forge "tensorflow>=2.8.0"
conda install -c conda-forge "keras>=2.8.0"
conda install -c conda-forge "numpy>=1.21.0"
conda install -c conda-forge "scikit-learn>=1.0.0"
```

## Usage

### Basic Usage

```bash
# Activate Python 3.12
conda activate py312

# Predict a single FASTA file
python predict.py -f your_sequences.fasta

# Predict with custom output directory
python predict.py -f sequences.fa -o results/

# Auto-detect and predict all FASTA files in current directory
python predict.py --auto
```

### Advanced Usage

```bash
# Custom model path and batch size
python predict.py -f input.fasta -m path/to/model.h5 -b 64

# Quiet mode (less verbose output)
python predict.py -f sequences.fasta --quiet

# Help
python predict.py --help
```


# Model Output Example
High Confidence Predictions (≥0.90):
- skp (0.9923) - Protease/chaperone factor
- aroD (0.9884) - Metabolism/respiration factor
- purE (0.9407, 0.9433) - N5-carboxyaminoimidazole ribonucleotide mutase

Lower Confidence Predictions:
- aroG (0.8347)
- hns (0.6074)

# Applications
- Probiotic strain optimization
- Identification of lifespan-modulating bacterial genes
- Assessment of potential probiotic safety
- Prediction of bacterial resistance mechanisms

# Interpretation
- Fewer predicted genes suggest fewer factors that could hinder C. elegans lifespan extension
- Gene deletions based on predictions may enhance probiotic effects
- Absence of certain genes may indicate pre-existing resistance mechanisms

# Reference Implementation
The model uses the Escherichia coli BW25113 strain as a reference, which contains 29 gene classes. Predictions in other strains are compared against this reference set.

# Usage in Probiotic Development
- Screening potential probiotic strains
- Identifying genetic modifications for enhanced probiotic effects
- Evaluating strain safety and resistance mechanisms
