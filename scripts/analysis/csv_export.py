#!/usr/bin/env python3
import math
import numbers
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from scripts.features import compute_pca_components

HC50_CALIBRATION_CUTOFF_UM = 150.0

PMF_FEATURE_COLUMNS = [
    "DeltaG_ads_kJ_per_mol",
    "DeltaG_ads_per_area_kJ_per_mol_nm2",
    "DeltaG_ads_per_res_kJ_per_mol",
    "DeltaG_barrier_kJ_per_mol",
    "Kp_ads_nm",
    "Kp_eff_nm",
    "theta_at_1uM",
]

PCA_VARIANCE_THRESHOLD = 0.95


def _clip_hc50(value):
    if value is None:
        return None
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(val):
        return None
    return min(val, HC50_CALIBRATION_CUTOFF_UM)

# Alle Runs finden
run_dirs = list(Path("simulations").rglob("pmf_analysis_results.yaml"))

records = []

for f in run_dirs:
    try:
        with open(f, "r") as fh:
            res = yaml.safe_load(fh)
    except Exception:
        continue

    # Metadaten
    run_dir = f.parent
    peptide_id = run_dir.parts[1]  # e.g. "solvia_8_run_1"
    replicate_id = run_dir.name    # e.g. "pmf_midplane" or replicate_x
    
    feat = res.get("features", {})
    ads = res.get("adsorption", {})
    qc  = (res.get("quality_metrics") or {}).get("qc_status", {})

    hc50_range = ads.get("hc50_uM_range") or [None, None]
    hc50_low_raw = hc50_range[0] if len(hc50_range) > 0 else None
    hc50_high_raw = hc50_range[1] if len(hc50_range) > 1 else None
    if hc50_low_raw is not None and hc50_high_raw is not None:
        hc50_mean_raw = (float(hc50_low_raw) + float(hc50_high_raw)) / 2.0
    else:
        hc50_mean_raw = None
    hc50_low = _clip_hc50(hc50_low_raw)
    hc50_high = _clip_hc50(hc50_high_raw)
    hc50_mean = _clip_hc50(hc50_mean_raw)
    hc50_exceeds = any(
        (
            val is not None
            and isinstance(val, numbers.Real)
            and not math.isnan(float(val))
            and float(val) > HC50_CALIBRATION_CUTOFF_UM
        )
        for val in (hc50_low_raw, hc50_high_raw, hc50_mean_raw)
    )

    records.append({
        "Peptide_ID": peptide_id,
        "Replicate_ID": replicate_id,
        "Rg_xy_nm": (res.get("peptide") or {}).get("rg_xy_nm"),
        "N_residues": feat.get("n_residues"),
        "DeltaG_ads_kJ_per_mol": feat.get("delta_g_ads"),
        "DeltaG_ads_per_area_kJ_per_mol_nm2": feat.get("delta_g_ads_per_area_kj_per_mol_nm2"),
        "DeltaG_ads_per_res_kJ_per_mol": feat.get("delta_g_ads_per_res_kj_per_mol"),
        "z_ads_nm": feat.get("z_ads"),
        "DeltaG_barrier_kJ_per_mol": feat.get("delta_g_barrier"),
        "z_barrier_nm": feat.get("z_barrier"),
        "Kp_ads_nm": ads.get("kp_ads_nm"),
        "Kp_eff_nm": ads.get("kp_eff_nm"),
        "HC50_pred_low_uM": hc50_low,
        "HC50_pred_high_uM": hc50_high,
        "HC50_pred_mean_uM": hc50_mean,
        "HC50_pred_low_uM_raw": hc50_low_raw,
        "HC50_pred_high_uM_raw": hc50_high_raw,
        "HC50_pred_mean_uM_raw": hc50_mean_raw,
        "HC50_pred_clipped_at_uM": HC50_CALIBRATION_CUTOFF_UM,
        "HC50_pred_exceeds_cutoff": hc50_exceeds,
        "HC50_rank_uM": ads.get("hc50_rank_uM") or feat.get("hc50_rank_uM"),
        "HC50_rank_score": feat.get("hc50_rank_score"),
        "HC50_rank_enh_uM": feat.get("hc50_rank_enh_uM"),
        "HC50_rank_enh_score": feat.get("hc50_rank_enh_score"),
        "HC50_rank_knn_uM": feat.get("hc50_rank_knn_uM"),
        "HC50_rank_knn_score": feat.get("hc50_rank_knn_score"),
        "HC50_rank_combo_uM": feat.get("hc50_rank_combo_uM"),
        "HC50_rank_combo_score": feat.get("hc50_rank_combo_score"),
        "HC50_rank_v3_uM": feat.get("hc50_rank_v3_uM"),
        "HC50_rank_v3_score": feat.get("hc50_rank_v3_score"),
        "theta_at_1uM": ads.get("theta_at_1uM"),
        "Footprint_area_nm2": ads.get("footprint_nm2"),
        "Footprint_cg_nm2": feat.get("footprint_cg_nm2"),
        "Rg_xy_cg_nm": feat.get("rg_xy_cg_nm"),
        "Sequence_length": feat.get("sequence_length"),
        "Sequence_length_cg": feat.get("sequence_length_cg"),
        "Sequence_charge": feat.get("sequence_charge"),
        "Sequence_charge_density": feat.get("sequence_charge_density"),
        "Sequence_isoelectric_point": feat.get("sequence_isoelectric_point"),
        "Sequence_hydrophobicity": feat.get("sequence_hydrophobicity"),
        "Sequence_hydrophobic_moment": feat.get("sequence_hydrophobic_moment"),
        "Sequence_logP": feat.get("sequence_logp"),
        "Sequence_logP_per_residue": feat.get("sequence_logp_per_residue"),
        "Sequence_helix_fraction": feat.get("sequence_helix_fraction"),
        "Sequence_helix_propensity": feat.get("sequence_helix_propensity"),
        "Sequence_positive_fraction": feat.get("sequence_positive_fraction"),
        "Sequence_negative_fraction": feat.get("sequence_negative_fraction"),
        "Sequence_aromatic_fraction": feat.get("sequence_aromatic_fraction"),
        "QC_pass": qc.get("passed"),
        "ESS_mean": (res.get("quality_metrics") or {}).get("ess_mean"),
        "ESS_min": (res.get("quality_metrics") or {}).get("ess_min"),
        "Overlap_min_adjacent": (res.get("quality_metrics") or {}).get("overlap_summary", {}).get("min_overlap_adjacent"),
        "JS_mean_adjacent": (res.get("quality_metrics") or {}).get("js_divergence_summary", {}).get("mean_js_adjacent"),
        "Ergodicity_pass_fraction": (res.get("quality_metrics") or {}).get("ergodicity_pass_fraction"),
        "Drift_max_nm_per_ns": None,  # Platzhalter
        "HC50_exp_uM": None,          # für experimentelle Werte
        "HC50_exp_CI_low_uM": None,
        "HC50_exp_CI_high_uM": None,
        "HC50_rank_resid_log10": feat.get("hc50_rank_resid_log10"),
        "HC50_rank_enh_resid_log10": feat.get("hc50_rank_enh_resid_log10"),
        "HC50_rank_knn_resid_log10": feat.get("hc50_rank_knn_resid_log10"),
        "HC50_rank_combo_resid_log10": feat.get("hc50_rank_combo_resid_log10"),
        "HC50_rank_v3_resid_log10": feat.get("hc50_rank_v3_resid_log10"),
    })

# DataFrame bauen und als CSV (und optional XLSX) speichern
df = pd.DataFrame.from_records(records)

# Ensure PMF feature columns exist even if missing for some runs
for col in PMF_FEATURE_COLUMNS:
    if col not in df.columns:
        df[col] = np.nan

pca_metadata = {}
numeric = df[PMF_FEATURE_COLUMNS].apply(pd.to_numeric, errors='coerce')
mask = numeric.notna().all(axis=1)
if mask.sum() >= 2:
    try:
        scores, loadings, variance_ratio, means = compute_pca_components(
            numeric.loc[mask].to_numpy(dtype=float),
            variance_threshold=PCA_VARIANCE_THRESHOLD,
        )
    except ValueError:
        pass
    else:
        component_cols = [f"PMF_PC{i+1}" for i in range(scores.shape[1])]
        for col in component_cols:
            df[col] = np.nan
        df.loc[mask, component_cols] = scores

        pca_metadata = {
            "source_columns": PMF_FEATURE_COLUMNS,
            "variance_threshold": PCA_VARIANCE_THRESHOLD,
            "means": means.tolist(),
            "components": [],
        }
        for idx, col in enumerate(component_cols):
            load_dict = {
                PMF_FEATURE_COLUMNS[row_idx]: float(loadings[row_idx, idx])
                for row_idx in range(loadings.shape[0])
            }
            pca_metadata["components"].append(
                {
                    "name": col,
                    "variance_ratio": float(variance_ratio[idx]),
                    "loadings": load_dict,
                }
            )

        out_dir = Path("analysis")
        out_dir.mkdir(exist_ok=True)
        try:
            decor_df = df.loc[:, ["Peptide_ID", "Replicate_ID", *component_cols]].copy()
            decor_df.to_csv(out_dir / "pmf_features_decorrelated.csv", index=False)
            with (out_dir / "pmf_feature_pca.yaml").open("w", encoding="utf-8") as handle:
                yaml.safe_dump(pca_metadata, handle, sort_keys=False)
            print(f"✓ Wrote decorrelated PMF features to {out_dir / 'pmf_features_decorrelated.csv'}")
        except Exception:
            pass

df.to_csv("results_summary.csv", index=False)

# Optional weiterhin Excel bereitstellen, falls ältere Pipelines darauf zugreifen
try:
    df.to_excel("results_summary.xlsx", index=False)
except Exception:
    pass

print(f"✅ Exported {len(records)} runs to results_summary.csv")
