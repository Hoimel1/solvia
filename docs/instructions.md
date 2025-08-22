# SOLVIA Repository Instructions

Last updated: August 18, 2025

## 🚀 Quick Start

This repository contains the SOLVIA framework for predicting hemolytic toxicity of antimicrobial peptides using coarse-grained molecular dynamics simulations.

### GPU-Enabled GROMACS

The system is configured with GROMACS 2023.3 with CUDA support:
- **Binary location**: `/home/michelhuller/solvia/software/gromacs-gpu-89/bin/gmx`
- **Symlink**: `/usr/local/bin/gmx` (GPU-enabled version)
- **GPU**: NVIDIA L4 (23GB VRAM)

## 📁 Repository Structure

```
solvia/
├── data/                    # Input/output data
│   ├── input/              # Input sequences and labels
│   ├── processed/          # Processed structures (PDB, CG-PDB)
│   └── raw/                # Raw peptide sequences
├── docs/                   # Documentation
├── force_fields/           # Martini 3 force field files
├── scripts/                # Universal Python scripts (cleaned & organized)
├── simulations/            # Simulation runs
│   ├── solvia2_run1_4pep_1nm/  # Active production run (100ns)
│   ├── archive/            # Old test runs
│   └── templates/          # MDP templates
├── software/               # Installed software
│   └── gromacs-gpu-89/    # GPU-enabled GROMACS
└── src/                    # Source code
    ├── api/               # FastAPI REST interface
    ├── ml/                # Machine learning models
    └── simulation/        # Simulation modules
```

## 🎮 Running Simulations

### 1. Structure Prediction (if needed)

```bash
python scripts/colabfold.py \
    --fasta data/raw/fasta/peptide.fasta \
    --output data/processed/pdb/peptide
```

### 2. Coarse-Graining

```bash
python scripts/coarse_grain.py \
    --pdb data/processed/pdb/peptide/peptide_rank_001.pdb \
    --output data/processed/cg_pdb \
    --posres  # Include position restraints
```

### 3. Create Peptide-Membrane System

```bash
python scripts/create_system.py \
    --config config/variant1_baseline.yaml \
    --peptide-gro data/processed/cg_pdb/peptide.gro \
    --peptide-itp data/processed/cg_pdb/peptide.itp \
    --output simulations/new_run
```

### 4. Run Equilibration

```bash
python scripts/equilibrate.py \
    --system-dir simulations/new_run \
    --gpu  # Use GPU acceleration
```

This runs:
- Energy minimization
- NVT equilibration (500 ps) with position restraints
- NPT equilibration (500 ps)

### 5. Run Production

```bash
python scripts/production.py \
    --system-dir simulations/new_run \
    --time-ns 100 \
    --gpu  # Use GPU acceleration
```

### 6. Analyze Trajectory

```bash
python scripts/analyze.py \
    --tpr simulations/new_run/prod.tpr \
    --xtc simulations/new_run/prod.xtc \
    --output simulations/new_run/analysis.json
```

### 7. Monitor Progress

```bash
# Check GPU usage
nvidia-smi

# Watch simulation log
tail -f simulations/new_run/prod.log

# Check trajectory size
ls -lh simulations/new_run/prod.xtc
```

## 🎯 Interactive Console

Run simulations interactively with the SOLVIA console:

```bash
cd /home/michelhuller/solvia
python solvia_console.py
```

The console provides:
- Interactive FASTA file selection
- Pre-configured simulation variants (Baseline, Escalation, Cholesterol)
- Automated pipeline execution via Snakemake
- Progress tracking and error handling

## 🔧 Key Scripts (Universal & Parametrized)

### Core Pipeline Scripts
- `scripts/colabfold.py` - Structure prediction with ColabFold
- `scripts/coarse_grain.py` - Martinize2 coarse-graining
- `scripts/create_system.py` - Create peptide-membrane systems
- `scripts/equilibrate.py` - Run equilibration (NVT/NPT)
- `scripts/production.py` - Run production simulations
- `scripts/analyze.py` - Analyze trajectories
- `scripts/aggregate_features.py` - Aggregate ML features
- `scripts/train_hemolysis_model.py` - Train ML models

### Interactive Tools
- `solvia_console.py` - Interactive simulation console
- `Snakefile` - Snakemake workflow definition

## 🐍 Python Environment

Activate the environment with required packages:
```bash
pip install -r requirements.txt
```

Key packages:
- MDAnalysis for trajectory analysis
- XGBoost for ML predictions
- Snakemake for workflow management
- FastAPI for REST API

## 🧬 Input Format

Place peptide sequences in FASTA format:
```
data/input/sequences/peptide_name.fasta
```

Example:
```
>SOLVIA_2
ALWKTLLKKVLKAAA
```

## 📊 Current Status

### Active Simulations
- **solvia2_run1_4pep_1nm**: 1.1 μs production run with 4 peptides at 1nm distance
  - Status: Completed
  - Started: Aug 14, 2025
  - Finished: Aug 18, 2025

### Naming Convention
All production runs follow the pattern:
```
solvia2_run{N}_{peptides}_{distance}
```

Example: `solvia2_run1_4pep_1nm` = Run 1 with 4 peptides at 1nm distance

## ⚠️ Important Notes

1. **Always use GPU GROMACS**: The system `gmx` (without path) defaults to GPU version
2. **Force fields**: All Martini 3 files are in `/home/michelhuller/solvia/force_fields/`
3. **Memory**: Large membrane systems require significant RAM
4. **Storage**: 100ns simulations generate ~10-20GB trajectory files

## 🐛 Troubleshooting

### GPU not detected
```bash
# Check GPU status
nvidia-smi

# Verify GROMACS GPU support
gmx --version | grep GPU
```

### Missing force field errors
Ensure topology includes all required Martini 3 files:
- martini_v3.0.0.itp
- martini_v3.0.0_phospholipids_v1.itp
- martini_v3.0_sterols_v1.0.itp
- etc.

### Simulation crashes
Check logs in simulation directory:
- `em.log` - Energy minimization
- `nvt_pr.log` - NVT equilibration
- `npt_pr.log` - NPT equilibration
- `prod.log` - Production run

## 📞 Support

For detailed scientific background, see:
- `docs/solvia_concept.md` - Framework overview
- `docs/troubleshooting.md` - Parameter optimization guide

For implementation details, check the source code in `src/`.
