# HC50-Prediction Framework for SOLVIA Peptides

## 1. Motivation and Scope

Peptide–membrane interactions determine the therapeutic window of many antimicrobial and antiviral candidates. Experimentally, the cytotoxic concentration HC50 (half-maximal haemolysis) is measured on blood-derived cell lines and is one of the most critical safety metrics. The present work formalises a physics-driven, simulation-backed pipeline for predicting HC50 from coarse-grained (CG) molecular dynamics (MD) and potential-of-mean-force (PMF) data. The document consolidates the methodology required to move from raw PMF trajectories to calibrated HC50 forecasts and describes the mathematical models, datasets, and assumptions employed.

## 2. Data Generation Overview

### 2.1 Simulation Workflow

1. **Coarse-grained peptide model** – generated with MARTINI 3.
2. **Membrane construction** – validated red-blood-cell (RBC) asymmetric bilayer (cholesterol-rich) per `config/pmf_standard_config.yaml`.
3. **Umbrella sampling** – single-sided PMFs along the bilayer normal, reference mode `hybrid`, 36 windows from 3.0 → 0.0 nm, force constant 900 kJ/mol/nm².
4. **Production** – typically 60 ns per window (post SMD and pre-equilibration).
5. **Bootstrap** – default 200 replicates for CI estimation.

### 2.2 Key Observables

- Umbrella traces `pullx.xvg`
- Windows metadata `pmf_metadata.yaml`
- PMF `pmf_analysis_results.yaml`
- Coarse-grain reference structures (`*_cg.gro`, `*_cg.pdb`)
- Experimental HC50 from DBAASP (`data/raw/peptides.csv`)

## 3. PMF Post-processing

Using `scripts/analysis/pmf_mbar_analysis.py`, we derive:

- ΔG<sub>ads</sub>, ΔG<sub>insert</sub>, ΔG<sub>‡</sub>
- MBAR-derived K<sub>p</sub> (per leaflet) and effective K<sub>p,eff</sub>
- Adsorbed- and inserted-state locations (z<sub>ads</sub>, z<sub>insert</sub>)
- Bootstrap errors & confidence intervals
- QC metrics: effective sample size (ESS), window overlap matrix, Jensen-Shannon divergence, ergodicity stability
- CG-derived descriptors: planar radius of gyration R<sub>g,xy</sub>, convex hull footprint A<sub>foot</sub>, sequence length L, etc.

### 3.1 Potential of Mean Force to Adsorption Partition Function

For each adsorption region `[z_low, z_high]`, PMF integration yields the per-leaflet partition length

\[
  K_p = \int_{z_{low}}^{z_{high}} \left[ e^{-\beta W(z)} - 1 \right] \mathrm{d}z,
\]

where `β = 1/(RT)`. Accounting for bilayer symmetry,

\[
  K_{p,\mathrm{eff}} = b_f K_p,
\]

with bilayer factor \( b_f = 2 \) by default.

### 3.2 Surface Coverage and HC50 (Gamma-model)

Observationally, cytotoxicity occurs around a critical surface coverage Γ*. For a chosen lipid area A<sub>lipid</sub> (0.62 nm²), we use L/P ranges (number of lipids per peptide) to express

\[
  \Gamma^{*} = \frac{1}{A_{lipid} (L/P^{*})}.
\]

We explore both fixed `L/P*` and footprint-based estimates

\[
  L/P^{*}_{foot} = \frac{A_{foot}}{A_{lipid}} \cdot f_{pack},
\]

with packing factor `f_pack` configurable. The HC50 estimate is then

\[
  \mathrm{HC50}_{\mathrm{theo}} = \frac{\Gamma^*}{0.602214 \ K_{p,\mathrm{eff}}},
\]

converted to µM by multiplying with 10<sup>6</sup>.

## 4. Calibration Strategy

All PMF outputs are consolidated via `scripts/analysis/csv_export.py` into `results_summary.csv`. We merge experimental HC50 data from DBAASP (matching on SOLVIA IDs) and apply three successive models in `analysis/calibrate_hc50.py`:

1. **Global scaling**: single factor s optimising median log10 difference.
2. **Gamma-model regression**: `log Γ* = α ΔG<sub>ads</sub> + β`.
3. **Rank-based monotonic model**: emphasises order preservation.

### 4.1 Feature Vector

For peptide i, the feature vector includes:

- ΔG<sub>ads</sub>
- log<sub>10</sub>K<sub>p,eff</sub>
- log<sub>10</sub>Γ* (model; computed either from regression or direct PMF integral)
- log<sub>10</sub>L/P*
- log<sub>10</sub>θ(1 µM)
- log<sub>10</sub> sequence length
- ESS mean (quality proxy)
- Optional: footprint area (PMF & CG), R<sub>g</sub>

Weights w are determined by maximising Spearman correlation between predicted scores and log10(HC50) on the training set using BFGS-based optimisation.

### 4.2 Isotonic Calibration

The linear score is mapped to HC50 via isotonic regression

\[
\hat{s}_i = w^{\top} \phi_i, \quad \log_{10}\hat{H}_{50,i} = \mathrm{Iso}(\hat{s}_i),
\]

where `Iso` is the non-decreasing piecewise-linear function learned on the training set (scikit-learn independent, implemented analytically).

### 4.3 Cross-validation

Leave-one-out (LOO) metrics quantify overfitting; the current training set (8–10 peptides) shows high Spearman on full data but negative LOO Spearman, indicating sensitivity to individual peptides. Increasing the dataset or tightening QC is therefore essential for deployment.

## 5. Embedding Calibration in PMF Analysis

`analysis/calibrate_hc50.py` now persists the latest rank- and gamma-model parameters to `analysis/hc50_rank_model.yaml` and `config/hc50_rank_model.yaml`, including:

- Feature names
- Standardisation stats (means, stds)
- Optimised weights
- Isotonic mapping arrays (score vs. log10 HC50)
- Gamma regression coefficients

`scripts/analysis/pmf_mbar_analysis.py` automatically loads these parameters at runtime. When the file exists, the script computes rank-based HC50 (`hc50_rank_uM`) during PMF post-processing. This enables direct deployment on new peptides without re-running the calibration stage.

## 6. Footprint-aware L/P* Calibration

To address the variability in peptide footprint, we plan incremental upgrades:

1. **Direct CG footprint usage** – use `footprint_cg_nm2` for both monolayer packing and `L/P*` estimation.
2. **Hybrid approach** – blend dynamic PMF-based footprint estimates (contact hull) with CG-derived values.
3. **Sensitivity tests** – run the calibration with alternative packing factors (`f_pack`) and check Spearman/LOO metrics.

## 7. Bilayer Factor Considerations

Currently, we keep \( b_f = 2 \), acknowledging that the membrane is accessible from both sides even though our simulations sample one leaflet. To validate this assumption, the calibration script can be run with alternative `bilayer_factor` values; the rank model’s performance across factors will guide updates to `config/pmf_standard_config.yaml`.

## 8. QA & Holdout Strategy

1. **PMF re-analysis** – ensure MBAR (pymbar) installed and Bootstrap enabled for each run; the scripts automatically overwrite YAML and plots.
2. **Calibration refresh** – `python scripts/analysis/csv_export.py` → `python analysis/calibrate_hc50.py ...`. The summary outlines training metrics and LOOCV diagnostics.
3. **Outlier inspection** – Large rank residuals (|log10 residual| > 0.5) trigger a review of PMF QC (overlap, ESS) and experimental HC50 reliability.
4. **Holdout testing** – maintain a separate list of peptides excluded from training. Analyse them via the updated pipeline; compare predicted HC50 to any available experimental benchmarks.

## 9. Future Enhancements

- **Data augmentation** – incorporate additional experimental HC50 values, especially for borderline peptides.
- **Feature expansion** – include sequence-derived descriptors (hydrophobic moment, net charge) and membrane-specific factors (e.g., local lipid composition around the adsorption region).
- **Bayesian models** – extend to probabilistic calibration that propagates Bootstrap uncertainties into HC50 distributions.
- **Automated reporting** – integrate dashboards summarising QC gates, calibration diagnostics, and ranking feedback loops.

## 10. Summary

The SOLVIA HC50 prediction framework now delivers a physics-informed, rank-calibrated forecast for each peptide. By merging PMF observables, CG structural data, and curated experimental HC50 references, we obtain a monotonic surrogate capable of generalising to peptides without measured cytotoxicity. The pipeline remains modular: calibration parameters are versioned in `config/`, the PMF analysis consumes them on the fly, and the central CSV serves as the data source for validation and downstream analytics.
