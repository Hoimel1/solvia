# SOLVIA
**S**tructure-based M**OL**ecular dynamics-driven hemolysis predictor **V**ia **I**ntegrative **A**I

MD-driven hemolysis predictor for antimicrobial peptides

## Overview
SOLVIA integrates coarse-grained molecular dynamics (CG-MD) simulations with interpretable machine learning to predict hemolytic toxicity in antimicrobial peptides (AMPs), enabling high-throughput screening for safer antibiotic development.

## Project Status
- [x] Project setup and data preparation
- [x] Module 1: Structure prediction (ColabFold) - Complete
- [x] Module 2: Coarse-graining (martinize2) - Complete
- [x] Module 3: MD simulations (GROMACS) - Setup Complete
- [ ] Module 4: Feature extraction - Next
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
- ColabFold (for structure prediction) - ✓ Installed
- Martinize2 (for coarse-graining) - ✓ Installed
- GROMACS (for MD simulations) - ✓ Installed
- INSANE (for membrane building) - Pending (alternative setup available)

## Module 1: Structure Prediction (ColabFold)

ColabFold has been installed locally in `localcolabfold/`. To run structure predictions:

### Single sequence:
```bash
export PATH="/home/michelhuller/solvia/localcolabfold/localcolabfold/colabfold-conda/bin:$PATH"
colabfold_batch input.fasta output_dir
```

### Batch processing for all peptides:
```bash
# Test run (first 5 sequences only)
python scripts/batch_colabfold.py --test-run

# Full run (all 2012 peptides)
python scripts/batch_colabfold.py

# Custom parameters
python scripts/batch_colabfold.py \
    --num-models 5 \
    --num-recycle 3 \
    --batch-size 10
```

The script features:
- Automatic progress tracking and resume capability
- MSA queries with rate limiting
- GPU support (when available)
- Saves best-ranked PDB to `data/processed/pdb/`
- Detailed logging to `logs/colabfold_batch.log`

**Note**: Processing 2012 peptides will take significant time (~10-20 hours depending on GPU availability).

## Module 2: Coarse-Graining (Martinize2)

Martinize2 is installed via vermouth package. It converts atomistic PDB structures to coarse-grained representations for Martini simulations.

### Single peptide:
```bash
martinize2 -f input.pdb -o topol.top -x output_cg.pdb -ff martini3001
```

### Batch processing:
```bash
# Test run (processes available PDB files)
python scripts/batch_martinize2.py --test-run

# Full run (all PDB files in data/processed/pdb/)
python scripts/batch_martinize2.py

# Custom parameters
python scripts/batch_martinize2.py \
    --force-field martini3001 \
    --secondary-structure auto
```

The script features:
- Automatic secondary structure determination for short peptides
- Progress tracking and resume capability
- Organized topology output
- Saves CG-PDB to `data/processed/cg_pdb/`
- Saves topologies to `data/processed/topologies/`

**Results**: Average ~49 beads per peptide (from ~350 atoms)

## Module 3: MD Simulations (GROMACS)

GROMACS 2021.4 is installed for running molecular dynamics simulations.

### Force Field Setup
**Important**: Upload your Martini force field files to:
- `force_fields/martini3/` - For Martini 3 force field files
- `force_fields/martini2/` - For Martini 2 force field files (if needed)

Required files:
- `martini_v3.0.0_solvents_v1.itp` - Water and ion topologies
- `martini_v3.0.0_phospholipids_v1.itp` - Lipid topologies
- Additional ITP files as needed

### Single peptide simulation setup:
```bash
python scripts/setup_membrane_simulation.py \
    --peptide-id SOLVIA_1 \
    --cg-pdb data/processed/cg_pdb/SOLVIA_1_cg.pdb \
    --topology-dir data/processed/topologies \
    --output-dir simulations/systems
```

### Batch simulation setup:
```bash
# Test run (first 5 peptides)
python scripts/batch_simulation_setup.py --test-run

# Full run (all CG peptides)
python scripts/batch_simulation_setup.py

# Parallel processing
python scripts/batch_simulation_setup.py --num-workers 4
```

The setup creates:
- Box generation around peptide
- MDP files for energy minimization, equilibration (NVT/NPT), and production
- Job submission scripts
- Organized simulation directories in `simulations/systems/`

### Running simulations:
```bash
cd simulations/systems/SOLVIA_1
bash submit.sh  # Or submit to cluster queue
```

**Note**: Full membrane insertion requires INSANE script and force field files.

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
