# 🚀 SOLVIA Toxprediction Engine - PMF Implementation

## Overview

The SOLVIA repository has been transformed into a focused **PMF-based Toxprediction Engine** for predicting hemolytic toxicity of antimicrobial peptides. This implementation addresses all critical issues identified in `docs/aenderungen.md` and provides a robust, publication-ready pipeline.

## 🎯 Key Improvements Implemented

### 1. **Local Patch Reference** ✅
- **Problem**: Global leaflet COM caused systematic bias
- **Solution**: Dynamic local patch (2nm radius) around peptide COM
- **File**: `scripts/universal/08_run_pmf_enhanced.py`

### 2. **Standardized RBC Membrane** ✅
- **Problem**: Too many membrane options, poor reproducibility
- **Solution**: Single standard membrane (45% cholesterol)
- **File**: `scripts/universal/04_build_membrane_standard.py`

### 3. **Deterministic Replicates** ✅
- **Problem**: No systematic orientation/seed management
- **Solution**: R1=parallel, R2=tilted(30°), seeds from peptide ID
- **File**: `scripts/universal/05_insert_peptides_pmf.py`

### 4. **SMD Initialization** ✅
- **Problem**: Direct snap caused poor sampling
- **Solution**: SMD ladder for window preparation (5ns)
- **File**: `scripts/universal/08_run_pmf_enhanced.py`

### 5. **MBAR/WHAM Analysis** ✅
- **Problem**: No analysis pipeline
- **Solution**: Complete MBAR with bootstrap CIs
- **File**: `scripts/analysis/pmf_mbar_analysis.py`

### 6. **Quality Control System** ✅
- **Problem**: No QC gates
- **Solution**: Comprehensive QC with auto-correction suggestions
- **File**: `scripts/analysis/pmf_qc_system.py`

## 📁 New Files Created

### Core Pipeline Scripts
```
scripts/universal/
├── 04_build_membrane_standard.py  # Standardized membrane builder
├── 05_insert_peptides_pmf.py      # PMF-optimized peptide insertion
├── 06_equilibrate_pmf.py          # Controlled equilibration protocol
└── 08_run_pmf_enhanced.py         # Enhanced PMF with local patch

scripts/analysis/
├── pmf_mbar_analysis.py           # MBAR/WHAM analysis pipeline
└── pmf_qc_system.py               # Quality control system
```

### Configuration
```
config/
└── pmf_standard_config.yaml       # PMF-specific configuration

docs/
└── usage_pmf.md                   # Complete PMF usage guide
```

### Master Pipeline
```
run_pmf_pipeline.sh                # Automated complete workflow
```

## 🚀 Quick Start

Run complete PMF analysis for a single peptide:

```bash
# Automated pipeline for peptide ID 1
./run_pmf_pipeline.sh 1

# With options
./run_pmf_pipeline.sh 1 --replicates 2 --skip-colabfold
```

## 📊 Expected Results

For each peptide, the pipeline generates:

### Thermodynamic Features
- **ΔG_ads**: Free energy of adsorption
- **ΔG_insert**: Free energy of insertion
- **ΔG‡**: Activation barrier

### Toxicity Classification
Peptides are classified as toxic if:
- ΔG_ads < -8 kJ/mol AND
- (ΔG_insert ≤ -3 kJ/mol OR ΔG‡ < 12 kJ/mol)

### Quality Metrics
- Window overlap matrix
- Convergence analysis
- Effective sample size
- Replicate consistency

## 🔬 Scientific Validation

The pipeline targets:
- **Spearman ρ ≥ 0.5** for HC50 correlation
- **ROC-AUC ≥ 0.80** for binary classification
- **95% CI** through bootstrap (1000 samples)

## 📈 Performance

### Per Peptide
- Simulation: ~0.64-0.96 μs
- Time: 24-48 hours on GPU
- Storage: ~50 GB

### Full Study (42 peptides)
- Total: ~27-40 μs
- Timeline: 2-3 weeks with 4 GPUs

## 🛠️ Technical Requirements

### Software
- GROMACS 2023.3 (GPU-enabled)
- Python 3.8+
- pymbar (for MBAR analysis)

### Hardware
- NVIDIA GPU (L4 or better)
- 32+ GB RAM
- 2+ TB storage

## 📋 Workflow Steps

1. **Build standard membrane** (once)
2. **Process peptide structure**
3. **Insert with orientations** (2 replicates)
4. **Equilibrate** (6.5 ns protocol)
5. **Run PMF** (16-22 windows, 20 ns each)
6. **Quality control** (automatic)
7. **MBAR analysis** (with bootstrap)

## 🔍 Quality Control

The pipeline includes automatic QC checks:
- Minimum neighbor overlap: 10%
- Target overlap: 20% (auto-densify)
- ESS: ≥200 frames
- Convergence: ≤2 kJ/mol half-time
- Replicate consistency: ≤2 kJ/mol

## 📝 Publications

This implementation is designed for publication-ready results with:
- Reproducible seeds
- Comprehensive QC
- Bootstrap confidence intervals
- Automated reporting

## 🚨 Important Notes

1. **Always use the enhanced PMF script** (`08_run_pmf_enhanced.py`)
2. **Enable lipid z-restraints** for membrane stability
3. **Run 2 replicates minimum** for validation
4. **Check QC reports** before analysis

## 📧 Usage

For detailed usage, see:
- Complete guide: `docs/usage_pmf.md`
- Configuration: `config/pmf_standard_config.yaml`
- Master pipeline: `./run_pmf_pipeline.sh --help`

## ✅ Summary

The SOLVIA repository is now a focused, robust PMF-based toxprediction engine with:
- **Local patch reference** (eliminates bias)
- **Standardized protocols** (ensures reproducibility)
- **Comprehensive QC** (guarantees quality)
- **Automated pipeline** (reduces errors)
- **Publication-ready output** (with CIs)

This implementation addresses all issues from `docs/aenderungen.md` and provides a reliable workflow for predicting hemolytic toxicity of antimicrobial peptides.
