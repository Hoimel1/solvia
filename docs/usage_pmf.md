# 🎯 SOLVIA PMF-Optimized Pipeline for Hemolytic Toxicity Prediction

## Executive Summary

SOLVIA has been transformed into a focused **Toxprediction Engine** using robust PMF (Potential of Mean Force) calculations to predict hemolytic toxicity of peptides. This pipeline implements all critical improvements from the scientific analysis to ensure reproducible, publication-ready results.

## 🚀 Quick Start - Complete PMF Workflow

```bash
# 1. Setup peptide structure
python3 scripts/universal/01_setup_run.py data/raw/fasta/SOLVIA_1.fasta

# 2. Coarse-grain the peptide
python3 scripts/universal/03_coarse_grain.py simulations/solvia_1_run_1

# 3. Build standard RBC membrane (only once)
python3 scripts/universal/04_build_membrane_standard.py --output-dir membrane_standard

# 4. Insert peptide with proper orientation (2 replicates)
python3 scripts/universal/05_insert_peptides_pmf.py simulations/solvia_1_run_1 --replicate 1
python3 scripts/universal/05_insert_peptides_pmf.py simulations/solvia_1_run_1 --replicate 2

# 5. Equilibrate with lipid restraints
python3 scripts/universal/06_equilibrate_pmf.py simulations/solvia_1_run_1 --replicate 1
python3 scripts/universal/06_equilibrate_pmf.py simulations/solvia_1_run_1 --replicate 2

# 6. Run PMF with local patch reference
python3 scripts/universal/08_run_pmf_enhanced.py simulations/solvia_1_run_1 --replicate 1
python3 scripts/universal/08_run_pmf_enhanced.py simulations/solvia_1_run_1 --replicate 2

# 7. Analyze with MBAR
python3 scripts/analysis/pmf_mbar_analysis.py simulations/solvia_1_run_1/pmf/replicate_1
python3 scripts/analysis/pmf_mbar_analysis.py simulations/solvia_1_run_1/pmf/replicate_2
```

## 📋 Key Features of the Enhanced Pipeline

### 1. **Local Patch Reference** (Critical Improvement)
- Uses phosphates within 2.0 nm radius of peptide COM
- Eliminates systematic bias from global leaflet COM
- Ensures z-coordinate reflects true local membrane interaction

### 2. **Standardized RBC Membrane**
- Fixed composition: 45% cholesterol (physiological RBC level)
- Reproducible setup with seed=42
- Pre-generated leaflet indices and z-restraints

### 3. **Deterministic Replicates**
- R1: Parallel orientation (0°)
- R2: Tilted orientation (30°)
- Seeds derived from peptide ID + replicate number

### 4. **Adaptive Window Spacing**
- Dense sampling at interface (-0.6 to +0.6 nm): 0.15 nm steps
- Automatic densification for overlap ≥ 0.20
- Extended range to -2.0 nm for complete insertion profile

### 5. **Quality Control Gates**
- Minimum neighbor overlap: 10%
- Effective sample size: ≥200 frames
- Half-time convergence: ≤2 kJ/mol
- Replicate consistency: ≤2 kJ/mol

## 🔬 Scientific Basis

### PMF Features for Toxicity Prediction

The pipeline extracts three key thermodynamic features:

1. **ΔG_ads**: Free energy of adsorption (surface binding)
2. **ΔG_insert**: Free energy of insertion (membrane penetration)
3. **ΔG‡**: Activation barrier for insertion

### Toxicity Criteria

Peptides are classified as hemolytically toxic if:
- ΔG_ads < -8 kJ/mol AND
- (ΔG_insert ≤ -3 kJ/mol OR ΔG‡ < 12 kJ/mol)

Expected performance:
- Spearman ρ ≥ 0.5 for correlation with HC50
- ROC-AUC ≥ 0.80 for binary classification

## 🏗️ Standard Operating Procedure (SOP)

### Step 1: Prepare Standard Membrane

```bash
# Run only once - creates reusable membrane template
python3 scripts/universal/04_build_membrane_standard.py

# Output files:
# membrane_standard/
#   ├── membrane.gro          # Structure
#   ├── membrane.top          # Topology
#   ├── index_leaflets.ndx    # OuterPO4/InnerPO4 groups
#   └── posre_lipid_z.itp     # Z-only restraints
```

### Step 2: Process Peptide

```bash
# Setup and coarse-grain
RUN_DIR="simulations/solvia_${ID}_run_1"
python3 scripts/universal/01_setup_run.py data/raw/fasta/SOLVIA_${ID}.fasta
python3 scripts/universal/03_coarse_grain.py ${RUN_DIR}
```

### Step 3: Create PMF Systems (2 Replicates)

```bash
# Replicate 1 - Parallel orientation
python3 scripts/universal/05_insert_peptides_pmf.py ${RUN_DIR} \
  --replicate 1 --peptide-id SOLVIA_${ID} --posres-lipid-z

# Replicate 2 - Tilted orientation (30°)
python3 scripts/universal/05_insert_peptides_pmf.py ${RUN_DIR} \
  --replicate 2 --peptide-id SOLVIA_${ID} --posres-lipid-z
```

### Step 4: Equilibration Protocol

```bash
# For each replicate
for REP in 1 2; do
  python3 scripts/universal/06_equilibrate_pmf.py ${RUN_DIR} --replicate ${REP}
done

# Equilibration sequence (6.5 ns total):
# 1. Energy minimization with restraints
# 2. NVT (0.5 ns) - gentle heating from 100K
# 3. NPT pre-equilibration (1 ns) - with restraints
# 4. NPT production (5 ns) - restraints on lipids only
```

### Step 5: PMF Calculation with QC

```bash
# Run umbrella sampling for each replicate
for REP in 1 2; do
  python3 scripts/universal/08_run_pmf_enhanced.py ${RUN_DIR} --replicate ${REP}
done

# Features:
# - Local patch reference (2 nm radius)
# - ~16-22 windows from +2.8 to -2.0 nm
# - 20 ns production per window
# - Automatic overlap checking
# - Deterministic seeds
```

### Step 6: MBAR Analysis

```bash
# Analyze each replicate
for REP in 1 2; do
  python3 scripts/analysis/pmf_mbar_analysis.py \
    ${RUN_DIR}/pmf/replicate_${REP} \
    --method mbar --bootstrap 1000
done

# Outputs:
# - pmf_analysis_results.yaml  # Features and metrics
# - analysis_plots/
#   ├── pmf_profile.png        # PMF with 95% CI
#   ├── overlap_matrix.png     # Window overlap heatmap
#   └── convergence.png        # Time convergence check
```

## 📊 Batch Processing for Multiple Peptides

### Automated Pipeline for 42-Peptide Study

```bash
#!/bin/bash
# process_all_peptides.sh

# Peptide lists from study
TOXIC_PEPTIDES=(1 8 14 215 164 126 68 32 482 490 515 524 527 617 624 850 858 941 974 1023 1045)
NONTOXIC_PEPTIDES=(12 398 1051 1125 1219 1315 1343 1363 1564 1587 1663 1680 1684 1743 1844 1941 1952 1962 2012 1115 794)

# Process all peptides
for ID in "${TOXIC_PEPTIDES[@]}" "${NONTOXIC_PEPTIDES[@]}"; do
    echo "Processing SOLVIA_${ID}"
    
    # Full pipeline
    ./run_pmf_pipeline.sh ${ID}
    
    # Collect results
    cp simulations/solvia_${ID}_run_1/pmf/*/pmf_analysis_results.yaml \
       results/SOLVIA_${ID}_results.yaml
done

# Aggregate features
python3 scripts/analysis/aggregate_pmf_features.py results/
```

## 🔍 Quality Control Checkpoints

### Window Overlap Matrix
```python
# Minimum acceptable overlap: 0.10
# Target overlap: 0.20
# If overlap < 0.10: Add intermediate window automatically
```

### Convergence Criteria
```python
# ESS (Effective Sample Size): ≥ 200 frames
# Half-time check: |ΔG_half1 - ΔG_half2| ≤ 2 kJ/mol
# If not converged: Extend by 10 ns
```

### Replicate Consistency
```python
# |ΔG_ads(R1) - ΔG_ads(R2)| ≤ 2 kJ/mol
# |ΔG_insert(R1) - ΔG_insert(R2)| ≤ 2 kJ/mol
# If inconsistent: Review and possibly re-run
```

## 📈 Expected Outputs

### Per-Peptide Results

```yaml
# pmf_analysis_results.yaml
features:
  delta_g_ads: -12.3      # kJ/mol
  delta_g_insert: -5.7    # kJ/mol
  delta_g_barrier: 8.2    # kJ/mol
  z_surf_min: 0.35        # nm
  z_head_min: -0.82       # nm
  z_barrier: -0.15        # nm

quality_metrics:
  mean_overlap: 0.23
  min_overlap: 0.12
  convergence: -0.8       # kJ/mol

classification:
  toxic_prediction: true  # Based on thresholds
  confidence: 0.89
```

### Study-Level Analysis

```yaml
# study_results.yaml
correlation:
  spearman_rho: 0.62
  p_value: 0.001
  confidence_interval: [0.45, 0.78]

classification:
  roc_auc: 0.85
  sensitivity: 0.90
  specificity: 0.75
  
validated_thresholds:
  delta_g_ads: -8.0
  delta_g_insert: -3.0
  delta_g_barrier: 12.0
```

## 💻 Computational Requirements

### Per Peptide
- **Umbrella sampling**: ~0.64 μs (16 windows × 20 ns × 2 replicates)
- **With QC extensions**: up to 0.96 μs
- **GPU time**: ~24-48 hours on NVIDIA L4

### Full Study (42 peptides)
- **Total simulation**: ~27-40 μs
- **Storage**: ~2 TB for trajectories
- **Timeline**: 2-3 weeks with 4 GPUs

## 🚨 Troubleshooting

### Common Issues and Solutions

#### 1. Low Window Overlap
```bash
# Check overlap matrix
grep "overlap" pmf/replicate_1/qc_report.yaml

# Solution: Script automatically adds windows
# Or manually add intermediate window:
python3 scripts/universal/run_single_window.py --center -0.25
```

#### 2. Poor Convergence
```bash
# Extend simulation time
python3 scripts/universal/extend_window.py --window z_+0.000 --add-ns 10
```

#### 3. Replicate Divergence
```bash
# Compare PMF profiles
python3 scripts/analysis/compare_replicates.py ${RUN_DIR}/pmf

# Check for:
# - Different binding modes
# - Sampling issues
# - Need for third replicate
```

## 📚 Key Configuration Files

### pmf_standard_config.yaml
```yaml
pmf:
  umbrella:
    ref_mode: "patch"       # Critical: local reference
    patch_radius: 2.0       # nm
    force_constant: 900     # kJ/mol/nm²
    production_ns: 20.0
    
  qc:
    min_neighbor_overlap: 0.10
    target_overlap: 0.20
    min_ess_frames: 200
    
  analysis:
    method: "mbar"
    bootstrap:
      enabled: true
      n_bootstrap: 1000
```

## ✅ Publication Readiness Checklist

- [ ] All peptides processed with 2 replicates
- [ ] QC gates passed for all windows
- [ ] Replicate consistency verified
- [ ] Features extracted with confidence intervals
- [ ] Correlation analysis completed
- [ ] ROC curves generated
- [ ] Git commit hash recorded
- [ ] GROMACS version documented
- [ ] Seeds recorded for reproducibility

## 🎓 Scientific References

This pipeline implements the methodology from:
- PMF concept document: `docs/pmf_concept.md`
- Critical improvements: `docs/aenderungen.md`
- Validation data: `data/pmf_validation/pmf_study_metadata.yaml`

## 📧 Support

For issues or questions about the PMF pipeline:
1. Check error logs in `simulation_dir/logs/`
2. Review QC reports in `pmf/replicate_*/qc_report.yaml`
3. Verify configuration in `config/pmf_standard_config.yaml`
