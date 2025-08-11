# SOLVIA
**S**tructure-based M**OL**ecular dynamics-driven hemolysis predictor **V**ia **I**ntegrative **A**I

MD-driven hemolysis predictor for antimicrobial peptides

## Overview
SOLVIA integrates coarse-grained molecular dynamics (CG-MD) simulations with interpretable machine learning to predict hemolytic toxicity in antimicrobial peptides (AMPs), enabling high-throughput screening for safer antibiotic development.

## Project Status
- [x] Project setup and data preparation
- [ ] Module 1: Structure prediction (ColabFold) - In Progress
- [ ] Module 2: Coarse-graining (martinize2)
- [ ] Module 3: MD simulations (GROMACS)
- [ ] Module 4: Feature extraction
- [ ] Module 5: ML model training

## Installation

### Prerequisites
- Linux (Ubuntu 22.04+ recommended) 
- Python 3.9+
- NVIDIA GPU with CUDA 12.1+ (recommended)
- At least 32GB RAM

### Setup Steps

1. Clone the repository:
```bash
git clone https://github.com/[username]/solvia.git
cd solvia
```

2. Create and activate virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate
```

3. Install base dependencies:
```bash
pip install -r requirements.txt
```

4. Install specialized tools:
- ColabFold (for structure prediction)
- Martinize2 (for coarse-graining)
- GROMACS (for MD simulations)
- INSANE (for membrane building)

Detailed installation instructions coming soon...

## Data Structure
```
solvia/
├── data/
│   ├── raw/                 # Original data files
│   │   ├── fasta_split/     # Individual peptide FASTA files
│   │   ├── peptides.csv     # Peptide metadata
│   │   └── peptides.fasta   # Combined FASTA
│   └── processed/           # Processed data
│       ├── pdb/             # ColabFold structure predictions
│       └── cg_pdb/          # Coarse-grained structures
├── src/                     # Source code
├── scripts/                 # Processing scripts
├── tests/                   # Unit tests
└── logs/                    # Process logs
```

## Usage
Coming soon...

## Contributing
This project is under active development. Please check back for contribution guidelines.

## License
[License information to be added]
