#!/usr/bin/env python3
"""HC50 calibration pipeline with replicate aggregation and Gamma-model fit."""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy import stats, optimize
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from cycler import cycler
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import Pipeline

PUBLICATION_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "axes.linewidth": 1.2,
    "axes.grid": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 5.5,
    "ytick.major.size": 5.5,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "legend.frameon": False,
    "legend.fontsize": 11,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}

matplotlib.rcParams.update(PUBLICATION_STYLE)
plt.rcParams["axes.prop_cycle"] = cycler(color=[
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd",
    "#ff7f0e", "#17becf", "#8c564b", "#bcbd22",
])


def _format_publication_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='both', direction='out', length=5.0, width=1.0)


CAL_DATA_COLUMNS = [
    "Peptide_ID",
    "Replicate_ID",
    "Rg_xy_nm",
    "N_residues",
    "DeltaG_ads_kJ_per_mol",
    "DeltaG_ads_per_area_kJ_per_mol_nm2",
    "DeltaG_ads_per_res_kJ_per_mol",
    "z_ads_nm",
    "DeltaG_barrier_kJ_per_mol",
    "z_barrier_nm",
    "Kp_ads_nm",
    "Kp_eff_nm",
    "HC50_pred_low_uM",
    "HC50_pred_high_uM",
    "HC50_pred_mean_uM",
    "theta_at_1uM",
    "Footprint_area_nm2",
    "QC_pass",
    "ESS_mean",
    "ESS_min",
    "Overlap_min_adjacent",
    "JS_mean_adjacent",
    "Ergodicity_pass_fraction",
    "Drift_max_nm_per_ns",
    "HC50_exp_uM",
    "HC50_exp_CI_low_uM",
    "HC50_exp_CI_high_uM",
]

PEPTIDE_BASE_RE = re.compile(r"(.+)_run_\d+$")
CONV_MOLAR = 0.602214  # 1 M = 0.602214 nm^-3
CONFIG_DIR = Path(__file__).parent.parent / "config"

HC50_CENSOR_MAX = 150.0


def _resolve_hc50_cap(value: Optional[float]) -> float:
    if value is None or not np.isfinite(value):
        return HC50_CENSOR_MAX
    return float(value)


def _to_numeric_array(values: Any) -> np.ndarray:
    if isinstance(values, pd.Series):
        return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)


def _censor_exp(values: Any, cap: float) -> np.ndarray:
    arr = _to_numeric_array(values)
    return np.minimum(arr, cap)


def _censor_pair(exp_values: Any, pred_values: Any, cap: float) -> Tuple[np.ndarray, np.ndarray]:
    exp_raw = _to_numeric_array(exp_values)
    pred_raw = _to_numeric_array(pred_values)
    exp_clip = np.minimum(exp_raw, cap)
    pred_adj = np.minimum(pred_raw, cap)
    mask_high = np.isfinite(exp_raw) & (exp_raw > cap)
    pred_adj[mask_high] = cap
    return exp_clip, pred_adj


@dataclass
class CalibrationResults:
    scale_factor: float
    shift_median_log10: float
    metrics_before: Dict[str, Any]
    metrics_after: Dict[str, Any]
    weighted_scale_factor: Optional[float] = None
    weighted_shift_log10: Optional[float] = None
    calibration_df: Optional[pd.DataFrame] = None


@dataclass
class PowerLawResults:
    a: float
    b: float
    a_ci95: Tuple[float, float]
    b_ci95: Tuple[float, float]
    metrics: Dict[str, Any]
    n_points: int


@dataclass
class GammaModelResults:
    alpha: float
    beta: float
    alpha_ci95: Tuple[float, float]
    beta_ci95: Tuple[float, float]
    r_value: float
    metrics: Dict[str, Any]
    n_points: int


@dataclass
class RankModelResults:
    feature_names: Tuple[str, ...]
    weights: np.ndarray
    feature_means: np.ndarray
    feature_stds: np.ndarray
    score_sorted: np.ndarray
    log_sorted: np.ndarray
    spearman: float
    pearson: float
    rmse_log10: float
    mae_log10: float
    loo_spearman: Optional[float]
    loo_rmse_log10: Optional[float]
    n_points: int


@dataclass
class RankModelEnhancedResults:
    feature_names: Tuple[str, ...]
    coefficients: np.ndarray
    intercept: float
    feature_means: np.ndarray
    feature_stds: np.ndarray
    alpha: float
    spearman: float
    pearson: float
    rmse_log10: float
    mae_log10: float
    loo_spearman: Optional[float]
    loo_rmse_log10: Optional[float]
    n_points: int


@dataclass
class RankModelKNNResults:
    feature_names: Tuple[str, ...]
    k_neighbors: int
    feature_means: np.ndarray
    feature_stds: np.ndarray
    spearman: float
    pearson: float
    rmse_log10: float
    mae_log10: float
    loo_spearman: Optional[float]
    loo_rmse_log10: Optional[float]
    n_points: int


@dataclass
class RankModelEnsembleResults:
    component_names: Tuple[str, ...]
    coefficients: np.ndarray
    intercept: float
    spearman: float
    pearson: float
    rmse_log10: float
    mae_log10: float
    loo_spearman: Optional[float]
    loo_rmse_log10: Optional[float]
    n_points: int


@dataclass
class RankModelV3Results:
    feature_names: Tuple[str, ...]
    coefficients: np.ndarray
    intercept: float
    alpha: float
    l1_ratio: float
    scaler_mean: np.ndarray
    scaler_scale: np.ndarray
    iso_x: np.ndarray
    iso_y: np.ndarray
    spearman: float
    pearson: float
    rmse_log10: float
    mae_log10: float
    cv_spearman: Optional[float]
    cv_rmse_log10: Optional[float]
    cv_mae_log10: Optional[float]
    n_points: int
    sample_weight_mean: float


def _base_peptide(peptide_id: str) -> str:
    if not isinstance(peptide_id, str):
        return str(peptide_id)
    m = PEPTIDE_BASE_RE.match(peptide_id)
    return m.group(1) if m else peptide_id


def _to_bool(val: Any) -> Optional[bool]:
    if pd.isna(val):
        return None
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    if isinstance(val, (int, np.integer, float, np.floating)):
        return bool(val)
    text = str(val).strip().lower()
    if text in {"true", "t", "yes", "y", "pass", "passed", "ok"}:
        return True
    if text in {"false", "f", "no", "n", "fail", "failed"}:
        return False
    return None


def _geom_mean(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    vals = vals[vals > 0]
    if vals.empty:
        return float("nan")
    return float(np.exp(np.mean(np.log(vals))))


def _median(series: pd.Series) -> float:
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return float(vals.median()) if not vals.empty else float("nan")


def _safe_fraction(values: Iterable[Optional[bool]]) -> float:
    vals = [bool(v) for v in values if v is not None]
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def _safe_flag_any(series: Optional[pd.Series]) -> Optional[bool]:
    if series is None:
        return None
    s = pd.Series(series)
    if s.empty:
        return None
    numeric = pd.to_numeric(s, errors="coerce")
    if not numeric.notna().any():
        return None
    return bool(numeric.fillna(0).astype(bool).any())


def _build_rank_features_enhanced(df: pd.DataFrame) -> pd.DataFrame:
    """Construct expanded feature set for the enhanced rank model."""
    index = df.index

    def _series(name: str, default: float = np.nan) -> pd.Series:
        if name in df:
            return pd.to_numeric(df[name], errors='coerce')
        return pd.Series(default, index=index, dtype=float)

    feats = pd.DataFrame(index=index)
    delta_g = _series('DeltaG_ads_kJ_per_mol', 0.0)
    feats['delta_g'] = delta_g
    feats['delta_g_sq'] = delta_g ** 2

    def _log_clip(series: pd.Series, floor: float) -> pd.Series:
        arr = series.astype(float)
        return np.log10(np.clip(arr, floor, None))

    log_kp = _log_clip(_series('Kp_eff_nm', 1.0), 1e-9)
    feats['log_kp'] = log_kp
    feats['log_kp_sq'] = log_kp ** 2
    log_theta = _log_clip(_series('theta_at_1uM', 1.0), 1e-12)
    feats['log_theta'] = log_theta
    feats['log_theta_sq'] = log_theta ** 2
    log_gamma = _log_clip(_series('Gamma_star_model', 1.0), 1e-12)
    feats['log_gamma'] = log_gamma
    log_lp = _log_clip(_series('LP_star_model', 1.0), 1e-12)
    feats['log_lp'] = log_lp

    rank_score = _series('HC50_rank_score', 0.0).fillna(0.0)
    feats['rank_score'] = rank_score
    ess_mean = _series('ESS_mean_median', 0.0)
    feats['ess_mean'] = ess_mean
    feats['log_ess_mean'] = np.log10(np.clip(ess_mean.astype(float), 1e-6, None))
    seq = _series('sequence_length', 0.0)
    feats['log_seq_len'] = np.where(seq > 0, np.log10(seq.astype(float)), 0.0)
    feats['z_ads'] = _series('z_ads_nm', 0.0)
    feats['overlap_min'] = _series('Overlap_min_adjacent_median', 0.0)
    feats['js_mean'] = _series('JS_mean_adjacent_median', 0.0)

    # Interaction terms
    feats['delta_g_log_kp'] = delta_g * log_kp
    feats['delta_g_log_theta'] = delta_g * log_theta
    feats['delta_g_log_gamma'] = delta_g * log_gamma
    feats['log_kp_log_theta'] = log_kp * log_theta
    feats['log_kp_log_gamma'] = log_kp * log_gamma
    feats['rank_score_delta_g'] = rank_score * delta_g

    return feats.fillna(0.0)


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return float("nan")
    x_use = x[mask]
    y_use = y[mask]
    if np.allclose(x_use, x_use[0]) or np.allclose(y_use, y_use[0]):
        return float("nan")
    if method == "pearson":
        corr, _ = stats.pearsonr(x_use, y_use)
    else:
        corr, _ = stats.spearmanr(x_use, y_use)
    return float(corr)


def calc_metrics(log_exp: np.ndarray, log_pred: np.ndarray) -> Dict[str, Any]:
    mask = np.isfinite(log_exp) & np.isfinite(log_pred)
    if mask.sum() == 0:
        return {
            "n": 0,
            "rmse_log10": float("nan"),
            "mae_log10": float("nan"),
            "bias_median_log10": float("nan"),
            "pearson_r": float("nan"),
            "spearman_r": float("nan"),
        }
    resid = log_pred[mask] - log_exp[mask]
    return {
        "n": int(mask.sum()),
        "rmse_log10": float(np.sqrt(np.mean(resid**2))),
        "mae_log10": float(np.mean(np.abs(resid))),
        "bias_median_log10": float(np.median(resid)),
        "pearson_r": _safe_corr(log_pred, log_exp, method="pearson"),
        "spearman_r": _safe_corr(log_pred, log_exp, method="spearman"),
    }


def _isotonic_regression(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    z = y_sorted.astype(float).copy()
    w = np.ones_like(z)
    n = len(z)
    i = 0
    while i < n - 1:
        if z[i] <= z[i + 1] + 1e-12:
            i += 1
            continue
        total_weight = w[i] + w[i + 1]
        pooled = (w[i] * z[i] + w[i + 1] * z[i + 1]) / total_weight
        z[i] = z[i + 1] = pooled
        w[i] = w[i + 1] = total_weight
        j = i
        while j > 0 and z[j - 1] > z[j] + 1e-12:
            total_weight = w[j - 1] + w[j]
            pooled = (w[j - 1] * z[j - 1] + w[j] * z[j]) / total_weight
            z[j - 1] = z[j] = pooled
            w[j - 1] = w[j] = total_weight
            j -= 1
        i = max(j, 0)
    x_unique, indices = np.unique(x_sorted, return_index=True)
    z_unique = z[indices]
    return np.asarray(x_unique, dtype=float), np.asarray(z_unique, dtype=float)


def load_input_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix.lower() in {".csv"}:
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path, engine=None)

    optional_cols = {
        "HC50_rank_uM": np.nan,
        "HC50_rank_score": np.nan,
        "Footprint_cg_nm2": np.nan,
        "Rg_xy_cg_nm": np.nan,
        "Sequence_length": np.nan,
        "Sequence_length_cg": np.nan,
    }

    missing = [c for c in CAL_DATA_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            "Input file is missing required columns: " + ", ".join(missing)
        )

    for col, default in optional_cols.items():
        if col not in df.columns:
            df[col] = default

    for col in df.columns:
        if col in {"Peptide_ID", "Replicate_ID", "QC_pass"}:
            continue
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["Peptide_base"] = df["Peptide_ID"].apply(_base_peptide)
    df["QC_pass_bool"] = df["QC_pass"].map(_to_bool)
    return df


def load_metadata(meta_path: Optional[Path]) -> Optional[pd.DataFrame]:
    if meta_path is None:
        return None
    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")
    meta = pd.read_csv(meta_path)
    meta.columns = [c.strip() for c in meta.columns]
    meta['Peptide_key'] = meta['SOLVIA_ID'].astype(str).str.upper()
    meta['HC50_meta_uM'] = pd.to_numeric(meta.get('HC50_uM'), errors='coerce')
    if 'HC50_value' in meta.columns and 'HC50_unit' in meta.columns:
        hc50_val = pd.to_numeric(meta['HC50_value'], errors='coerce')
        hc50_unit = meta['HC50_unit'].astype(str).str.lower()
        conv = hc50_val.copy()
        conv.loc[hc50_unit.isin(['mg/l', 'mg per l', 'mg_per_l'])] = hc50_val / 1000.0
        conv.loc[hc50_unit.isin(['%'])] = np.nan
        meta['HC50_meta_uM'] = meta['HC50_meta_uM'].fillna(conv)
    return meta


def merge_with_metadata(
    agg_df: pd.DataFrame,
    meta_df: Optional[pd.DataFrame],
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    if meta_df is None:
        agg_df['HC50_meta_uM'] = np.nan
        agg_df['HC50_exp_source'] = agg_df['HC50_exp_uM'].apply(
            lambda x: 'analysis' if pd.notna(x) else None
        )
        return agg_df, {"n_meta": 0, "n_matched": 0}

    df = agg_df.copy()
    df['Peptide_key'] = df['Peptide_ID'].astype(str).str.upper()
    keep_cols = ['Peptide_key', 'HC50_meta_uM']
    for col in ['sequence_length', 'HC50_cell', 'NAME', 'DBAASP_ID']:
        if col in meta_df.columns:
            keep_cols.append(col)
    meta_merge = meta_df[keep_cols].copy()
    merged = df.merge(meta_merge, on='Peptide_key', how='left')

    if 'sequence_length' in merged.columns:
        merged['sequence_length'] = pd.to_numeric(merged['sequence_length'], errors='coerce')

    merged['HC50_exp_source'] = np.where(
        merged['HC50_exp_uM'].notna(), 'analysis',
        np.where(merged['HC50_meta_uM'].notna(), 'meta', None)
    )
    merged['HC50_exp_uM'] = merged['HC50_exp_uM'].fillna(merged['HC50_meta_uM'])

    info = {
        "n_meta": int(meta_df.shape[0]),
        "n_matched": int(merged['HC50_meta_uM'].notna().sum()),
    }
    return merged, info


def aggregate_replicates(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for base, group in df.groupby("Peptide_base", dropna=False):
        ess_vals = pd.to_numeric(group["ESS_mean"], errors="coerce").dropna()
        weights = np.maximum(1.0, ess_vals.values) if not ess_vals.empty else np.array([])
        records.append(
            {
                "Peptide_ID": base,
                "n_reps": int(len(group)),
                "replicate_ids": ";".join(sorted(group["Peptide_ID"].astype(str))),
                "replicate_tags": ";".join(sorted(group["Replicate_ID"].astype(str).unique())),
                "Rg_xy_nm": _median(group["Rg_xy_nm"]),
                "N_residues": _median(group["N_residues"]),
                "DeltaG_ads_kJ_per_mol": _median(group["DeltaG_ads_kJ_per_mol"]),
                "DeltaG_ads_per_area_kJ_per_mol_nm2": _median(
                    group["DeltaG_ads_per_area_kJ_per_mol_nm2"]
                ),
                "DeltaG_ads_per_res_kJ_per_mol": _median(
                    group["DeltaG_ads_per_res_kJ_per_mol"]
                ),
                "z_ads_nm": _median(group["z_ads_nm"]),
                "DeltaG_barrier_kJ_per_mol": _median(group["DeltaG_barrier_kJ_per_mol"]),
                "z_barrier_nm": _median(group["z_barrier_nm"]),
                "Kp_ads_nm": _geom_mean(group["Kp_ads_nm"]),
                "Kp_eff_nm": _geom_mean(group["Kp_eff_nm"]),
                "HC50_pred_low_uM": _geom_mean(group["HC50_pred_low_uM"]),
                "HC50_pred_high_uM": _geom_mean(group["HC50_pred_high_uM"]),
                "HC50_pred_mean_uM": _geom_mean(group["HC50_pred_mean_uM"]),
                "HC50_pred_low_uM_raw": _geom_mean(
                    group.get("HC50_pred_low_uM_raw", pd.Series(dtype=float))
                ),
                "HC50_pred_high_uM_raw": _geom_mean(
                    group.get("HC50_pred_high_uM_raw", pd.Series(dtype=float))
                ),
                "HC50_pred_mean_uM_raw": _geom_mean(
                    group.get("HC50_pred_mean_uM_raw", pd.Series(dtype=float))
                ),
                "HC50_rank_uM": _geom_mean(group.get("HC50_rank_uM", pd.Series(dtype=float))),
                "HC50_rank_score": _median(group.get("HC50_rank_score", pd.Series(dtype=float))),
                "theta_at_1uM": _geom_mean(group["theta_at_1uM"]),
                "Footprint_area_nm2": _median(group["Footprint_area_nm2"]),
                "Footprint_cg_nm2": _median(group.get("Footprint_cg_nm2", pd.Series(dtype=float))),
                "Rg_xy_cg_nm": _median(group.get("Rg_xy_cg_nm", pd.Series(dtype=float))),
                "HC50_pred_exceeds_cutoff": _safe_flag_any(
                    group.get("HC50_pred_exceeds_cutoff", pd.Series(dtype=float))
                ),
                "HC50_pred_clipped_at_uM": _median(
                    group.get("HC50_pred_clipped_at_uM", pd.Series(dtype=float))
                ),
                "QC_pass_fraction": _safe_fraction(group["QC_pass_bool"].tolist()),
                "ESS_mean_median": _median(group["ESS_mean"]),
                "ESS_weight_sum": float(weights.sum()) if weights.size else float("nan"),
                "ESS_weight_mean": float(weights.mean()) if weights.size else float("nan"),
                "ESS_min_median": _median(group["ESS_min"]),
                "Overlap_min_adjacent_median": _median(group["Overlap_min_adjacent"]),
                "JS_mean_adjacent_median": _median(group["JS_mean_adjacent"]),
                "Ergodicity_pass_fraction": _median(group["Ergodicity_pass_fraction"]),
                "Drift_max_nm_per_ns": _median(group["Drift_max_nm_per_ns"]),
                "HC50_exp_uM": _median(group["HC50_exp_uM"]),
                "HC50_exp_CI_low_uM": _median(group["HC50_exp_CI_low_uM"]),
                "HC50_exp_CI_high_uM": _median(group["HC50_exp_CI_high_uM"]),
            }
        )

    agg_df = pd.DataFrame.from_records(records)
    agg_df.sort_values("Peptide_ID", inplace=True)
    agg_df.reset_index(drop=True, inplace=True)
    return agg_df


def fit_scale_1p(
    agg_df: pd.DataFrame,
    weight_col: str = "ESS_weight_sum",
    calibration_max_um: Optional[float] = None,
    hc50_cap: Optional[float] = None,
) -> Tuple[pd.DataFrame, CalibrationResults]:
    df = agg_df.copy()
    mask = (
        df["HC50_exp_uM"].notna()
        & df["HC50_pred_mean_uM"].notna()
        & (df["HC50_exp_uM"].astype(float) > 0)
        & (df["HC50_pred_mean_uM"].astype(float) > 0)
    )
    cap = _resolve_hc50_cap(hc50_cap if hc50_cap is not None else calibration_max_um)
    max_uM = None
    if calibration_max_um is not None and np.isfinite(calibration_max_um):
        max_uM = float(calibration_max_um)
        mask &= df["HC50_exp_uM"].astype(float) <= max_uM
    df["calibration_cutoff_uM"] = max_uM if max_uM is not None else np.nan
    calib_df = df.loc[mask].copy()
    df["used_for_calibration"] = mask

    if calib_df.empty:
        df["HC50_cal_uM"] = df["HC50_pred_mean_uM"]
        df["HC50_resid_log10"] = np.nan
        empty_metrics = calc_metrics(
            np.array([], dtype=float), np.array([], dtype=float)
        )
        return df, CalibrationResults(
            scale_factor=1.0,
            shift_median_log10=float("nan"),
            metrics_before=empty_metrics,
            metrics_after=empty_metrics,
            weighted_scale_factor=None,
            weighted_shift_log10=None,
            calibration_df=calib_df,
        )

    exp_clip, pred_clip = _censor_pair(
        calib_df["HC50_exp_uM"],
        calib_df["HC50_pred_mean_uM"],
        cap,
    )
    log_exp = np.log10(np.clip(exp_clip, 1e-12, None))
    log_pred = np.log10(np.clip(pred_clip, 1e-12, None))
    delta = log_exp - log_pred
    shift_median = float(np.median(delta))
    scale_factor = float(10 ** shift_median)

    weighted_scale_factor = None
    weighted_shift = None
    if weight_col in calib_df.columns:
        weights = pd.to_numeric(calib_df[weight_col], errors="coerce").to_numpy()
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, np.nan)
        if np.isfinite(weights).sum() >= 2:
            valid = np.isfinite(weights)
            weighted_shift = float(np.average(delta[valid], weights=weights[valid]))
            weighted_scale_factor = float(10 ** weighted_shift)

    df["HC50_cal_uM"] = df["HC50_pred_mean_uM"] * scale_factor
    calib_df["HC50_cal_uM"] = calib_df["HC50_pred_mean_uM"] * scale_factor

    exp_raw_all = _to_numeric_array(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]))
    high_mask_all = np.isfinite(exp_raw_all) & (exp_raw_all > cap)
    df["HC50_cal_uM"] = np.minimum(df["HC50_cal_uM"].astype(float), cap)
    df.loc[high_mask_all, "HC50_cal_uM"] = cap

    exp_raw_calib = _to_numeric_array(calib_df.get("HC50_exp_uM_raw", calib_df["HC50_exp_uM"]))
    high_mask_calib = np.isfinite(exp_raw_calib) & (exp_raw_calib > cap)
    calib_df["HC50_cal_uM"] = np.minimum(calib_df["HC50_cal_uM"].astype(float), cap)
    calib_df.loc[high_mask_calib, "HC50_cal_uM"] = cap

    df["HC50_resid_log10"] = np.nan
    df.loc[mask, "HC50_resid_log10"] = (
        np.log10(np.clip(calib_df["HC50_cal_uM"].to_numpy(dtype=float), 1e-12, None))
        - log_exp
    )

    metrics_before = calc_metrics(log_exp, log_pred)
    log_cal = np.log10(np.clip(calib_df["HC50_cal_uM"].to_numpy(dtype=float), 1e-12, None))
    metrics_after = calc_metrics(log_exp, log_cal)

    calib_df["log10_resid_after"] = log_cal - log_exp
    calib_df["log10_resid_before"] = log_pred - log_exp

    results = CalibrationResults(
        scale_factor=scale_factor,
        shift_median_log10=shift_median,
        weighted_scale_factor=weighted_scale_factor,
        weighted_shift_log10=weighted_shift,
        metrics_before=metrics_before,
        metrics_after=metrics_after,
        calibration_df=calib_df,
    )
    return df, results


def fit_powerlaw_2p(
    agg_df: pd.DataFrame,
    calibration_max_um: Optional[float] = None,
    hc50_cap: Optional[float] = None,
) -> Tuple[pd.DataFrame, Optional[PowerLawResults]]:
    df = agg_df.copy()
    mask = (
        df["HC50_exp_uM"].notna()
        & df["Kp_eff_nm"].notna()
        & (df["HC50_exp_uM"].astype(float) > 0)
        & (df["Kp_eff_nm"].astype(float) > 0)
    )
    cap = _resolve_hc50_cap(hc50_cap if hc50_cap is not None else calibration_max_um)
    if calibration_max_um is not None and np.isfinite(calibration_max_um):
        max_uM = float(calibration_max_um)
        mask &= df["HC50_exp_uM"].astype(float) <= max_uM
    calib = df.loc[mask].copy()
    if calib.shape[0] < 3:
        return df, None

    exp_clip = _censor_exp(calib["HC50_exp_uM"], cap)
    log_exp = np.log10(np.clip(exp_clip, 1e-12, None))
    log_kp = np.log10(calib["Kp_eff_nm"].to_numpy(dtype=float))
    regression = stats.linregress(log_kp, log_exp)
    slope = float(regression.slope)
    intercept = float(regression.intercept)
    n = calib.shape[0]
    t_val = stats.t.ppf(0.975, df=n - 2) if n > 2 else float("nan")
    slope_ci = (
        slope - t_val * regression.stderr,
        slope + t_val * regression.stderr,
    )
    intercept_ci = (
        intercept - t_val * regression.intercept_stderr,
        intercept + t_val * regression.intercept_stderr,
    )

    a = intercept
    b = -slope
    a_ci = intercept_ci
    b_ci = (-slope_ci[1], -slope_ci[0])

    df.loc[df["Kp_eff_nm"] > 0, "HC50_cal_2p_uM"] = 10 ** (
        a - b * np.log10(df.loc[df["Kp_eff_nm"] > 0, "Kp_eff_nm"].astype(float))
    )
    calib["HC50_cal_2p_uM"] = 10 ** (
        a - b * np.log10(calib["Kp_eff_nm"].astype(float))
    )

    exp_raw_all = _to_numeric_array(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]))
    high_mask_all = np.isfinite(exp_raw_all) & (exp_raw_all > cap)
    df["HC50_cal_2p_uM"] = np.minimum(df["HC50_cal_2p_uM"].astype(float), cap)
    df.loc[high_mask_all, "HC50_cal_2p_uM"] = cap

    exp_raw_calib = _to_numeric_array(calib.get("HC50_exp_uM_raw", calib["HC50_exp_uM"]))
    high_mask_calib = np.isfinite(exp_raw_calib) & (exp_raw_calib > cap)
    calib["HC50_cal_2p_uM"] = np.minimum(calib["HC50_cal_2p_uM"].astype(float), cap)
    calib.loc[high_mask_calib, "HC50_cal_2p_uM"] = cap

    pred_clip_calib = calib["HC50_cal_2p_uM"].to_numpy(dtype=float)
    log_pred_2p = np.log10(np.clip(pred_clip_calib, 1e-12, None))
    metrics = calc_metrics(log_exp, log_pred_2p)

    df["HC50_resid_2p_log10"] = np.nan
    if mask.any():
        df.loc[mask, "HC50_resid_2p_log10"] = (
            np.log10(
                np.clip(df.loc[mask, "HC50_cal_2p_uM"].astype(float).to_numpy(), 1e-12, None)
            )
            - log_exp
        )

    result = PowerLawResults(
        a=a,
        b=b,
        a_ci95=a_ci,
        b_ci95=b_ci,
        metrics=metrics,
        n_points=n,
    )
    return df, result


def fit_gamma_model(
    agg_df: pd.DataFrame,
    area_per_lipid_nm2: float,
    calibration_max_um: Optional[float] = None,
    hc50_cap: Optional[float] = None,
) -> Tuple[pd.DataFrame, Optional[GammaModelResults]]:
    df = agg_df.copy()
    cap = _resolve_hc50_cap(hc50_cap if hc50_cap is not None else calibration_max_um)
    df["Gamma_star_obs"] = np.nan
    obs_mask = (
        df["HC50_exp_uM"].notna()
        & (df["HC50_exp_uM"].astype(float) > 0)
        & df["Kp_eff_nm"].notna()
        & (df["Kp_eff_nm"].astype(float) > 0)
    )
    if calibration_max_um is not None and np.isfinite(calibration_max_um):
        max_uM = float(calibration_max_um)
        obs_mask &= df["HC50_exp_uM"].astype(float) <= max_uM
    if obs_mask.any():
        df.loc[obs_mask, "Gamma_star_obs"] = (
            CONV_MOLAR
            * df.loc[obs_mask, "Kp_eff_nm"].astype(float)
            * df.loc[obs_mask, "HC50_exp_uM"].astype(float)
            * 1e-6
        )

    fit_mask = obs_mask & df["DeltaG_ads_kJ_per_mol"].notna()
    if fit_mask.sum() < 3:
        df["Gamma_star_model"] = np.nan
        df["HC50_gamma_uM"] = np.nan
        df["HC50_gamma_resid_log10"] = np.nan
        df["LP_star_model"] = np.nan
        df["LP_star_obs"] = np.nan
        df["LP_star"] = np.nan
        return df, None

    x = df.loc[fit_mask, "DeltaG_ads_kJ_per_mol"].astype(float)
    y = np.log(df.loc[fit_mask, "Gamma_star_obs"].astype(float))
    regression = stats.linregress(x, y)
    alpha = float(regression.slope)
    beta = float(regression.intercept)
    n = fit_mask.sum()
    t_val = stats.t.ppf(0.975, df=n - 2) if n > 2 else float("nan")
    alpha_ci = (
        alpha - t_val * regression.stderr,
        alpha + t_val * regression.stderr,
    )
    beta_ci = (
        beta - t_val * regression.intercept_stderr,
        beta + t_val * regression.intercept_stderr,
    )

    df["Gamma_star_model"] = np.nan
    delta_mask = df["DeltaG_ads_kJ_per_mol"].notna()
    df.loc[delta_mask, "Gamma_star_model"] = np.exp(
        alpha * df.loc[delta_mask, "DeltaG_ads_kJ_per_mol"].astype(float) + beta
    )
    df.loc[df["Gamma_star_model"] <= 0, "Gamma_star_model"] = np.nan

    df["LP_star_model"] = np.nan
    df.loc[delta_mask, "LP_star_model"] = 1.0 / (
        area_per_lipid_nm2 * df.loc[delta_mask, "Gamma_star_model"].astype(float)
    )
    df["LP_star_obs"] = np.nan
    df.loc[obs_mask, "LP_star_obs"] = 1.0 / (
        area_per_lipid_nm2 * df.loc[obs_mask, "Gamma_star_obs"].astype(float)
    )
    df["LP_star"] = df["LP_star_model"]

    df["HC50_gamma_uM"] = np.nan
    model_mask = (
        delta_mask
        & df["Kp_eff_nm"].notna()
        & (df["Kp_eff_nm"].astype(float) > 0)
    )
    df.loc[model_mask, "HC50_gamma_uM"] = (
        df.loc[model_mask, "Gamma_star_model"].astype(float)
        / (CONV_MOLAR * df.loc[model_mask, "Kp_eff_nm"].astype(float))
        * 1e6
    )
    df.loc[df["HC50_gamma_uM"] <= 0, "HC50_gamma_uM"] = np.nan
    exp_raw_all = _to_numeric_array(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]))
    high_mask_all = np.isfinite(exp_raw_all) & (exp_raw_all > cap)
    df["HC50_gamma_uM"] = np.minimum(df["HC50_gamma_uM"].astype(float), cap)
    df.loc[high_mask_all, "HC50_gamma_uM"] = cap

    df["HC50_gamma_resid_log10"] = np.nan
    resid_mask = (
        model_mask
        & df["HC50_exp_uM"].notna()
        & (df["HC50_exp_uM"].astype(float) > 0)
    )
    if resid_mask.any():
        exp_clip_resid, pred_clip_resid = _censor_pair(
            df.loc[resid_mask, "HC50_exp_uM"],
            df.loc[resid_mask, "HC50_gamma_uM"],
            cap,
        )
        df.loc[resid_mask, "HC50_gamma_resid_log10"] = (
            np.log10(np.clip(pred_clip_resid, 1e-12, None))
            - np.log10(np.clip(exp_clip_resid, 1e-12, None))
        )

    exp_clip_fit, pred_clip_fit = _censor_pair(
        df.loc[fit_mask, "HC50_exp_uM"],
        df.loc[fit_mask, "HC50_gamma_uM"],
        cap,
    )
    log_exp = np.log10(np.clip(exp_clip_fit, 1e-12, None))
    log_gamma = np.log10(np.clip(pred_clip_fit, 1e-12, None))
    gamma_metrics = calc_metrics(log_exp, log_gamma)

    result = GammaModelResults(
        alpha=alpha,
        beta=beta,
        alpha_ci95=alpha_ci,
        beta_ci95=beta_ci,
        r_value=float(regression.rvalue),
        metrics=gamma_metrics,
        n_points=int(n),
    )
    return df, result


def _prepare_feature_matrix(
    df: pd.DataFrame,
    feature_series: Dict[str, pd.Series],
) -> Tuple[np.ndarray, np.ndarray, Tuple[str, ...]]:
    feature_names = []
    columns = []
    for name, series in feature_series.items():
        arr = pd.to_numeric(series, errors='coerce').to_numpy(dtype=float)
        if np.isfinite(arr).sum() < 3:
            continue
        if np.nanstd(arr) < 1e-8:
            continue
        columns.append(arr)
        feature_names.append(name)
    if not columns:
        raise ValueError("No usable features for rank model")
    X = np.vstack(columns).T
    means = np.nanmean(X, axis=0)
    stds = np.nanstd(X, axis=0)
    stds = np.where(stds < 1e-12, 1.0, stds)
    X = (X - means) / stds
    mask = np.any(np.isnan(X), axis=1)
    if mask.any():
        X[mask] = np.where(np.isnan(X[mask]), 0.0, X[mask])
    return X, np.vstack([means, stds]), tuple(feature_names)


def _prepare_rank_v3_features(df: pd.DataFrame) -> pd.DataFrame:
    base = _build_rank_features_enhanced(df).copy()
    if 'delta_g' in base.columns and 'sequence_charge' in base.columns:
        base['delta_g_x_charge'] = base['delta_g'] * base['sequence_charge']
    if 'delta_g_ads_per_area' in base.columns and 'sequence_hydrophobicity' in base.columns:
        base['area_x_hydrophobicity'] = base['delta_g_ads_per_area'] * base['sequence_hydrophobicity']
    if 'qc_frac' in base.columns and 'ess_mean' in base.columns:
        base['qc_weight'] = base['qc_frac'] * base['ess_mean']
    if 'log_kp' in base.columns and 'log_seq_len' in base.columns:
        base['log_kp_per_residue'] = base['log_kp'] / base['log_seq_len'].replace({0.0: np.nan})
    return base


def _build_feature_matrix(
    df: pd.DataFrame,
    feature_names: Tuple[str, ...],
    stats_arr: np.ndarray,
) -> np.ndarray:
    means = stats_arr[0]
    stds = stats_arr[1]
    columns = []
    for idx, name in enumerate(feature_names):
        if name == 'delta_g':
            arr = df["DeltaG_ads_kJ_per_mol"].to_numpy(dtype=float)
        elif name == 'log_kp':
            arr = np.log10(df["Kp_eff_nm"].astype(float).clip(lower=1e-12))
        elif name == 'log_gamma_model':
            arr = np.log10(df["Gamma_star_model"].astype(float).clip(lower=1e-12))
        elif name == 'log_lp_model':
            arr = np.log10(df["LP_star_model"].astype(float).clip(lower=1e-12))
        elif name == 'log_theta':
            arr = np.log10(df["theta_at_1uM"].astype(float).clip(lower=1e-12))
        elif name == 'log_seq_len':
            arr = np.log10(df["sequence_length"].astype(float).clip(lower=1e-12))
        elif name == 'qc_frac':
            arr = df["QC_pass_fraction"].to_numpy(dtype=float)
        elif name == 'ess_mean':
            arr = df["ESS_mean_median"].to_numpy(dtype=float)
        else:
            arr = np.zeros(df.shape[0], dtype=float)
        arr = np.where(np.isfinite(arr), arr, means[idx])
        columns.append(arr)
    X = np.vstack(columns).T
    stds = np.where(stds < 1e-12, 1.0, stds)
    X = (X - means) / stds
    return X


def _optimize_spearman(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    init, *_ = np.linalg.lstsq(X, y, rcond=None)

    def objective(w: np.ndarray) -> float:
        score = X @ w
        if np.allclose(score, score[0]):
            return 1.0
        rho = stats.spearmanr(score, y).correlation
        if not np.isfinite(rho):
            return 1.0
        return -float(rho)

    res = optimize.minimize(objective, init, method='BFGS')
    if res.success and np.isfinite(res.fun):
        return res.x
    return init


def fit_rank_model(agg_df: pd.DataFrame, hc50_cap: Optional[float] = None) -> Tuple[pd.DataFrame, Optional[RankModelResults]]:
    df = agg_df.copy()
    mask = df["HC50_exp_uM"].notna() & (df["HC50_exp_uM"].astype(float) > 0)
    data = df.loc[mask].copy()
    cap = _resolve_hc50_cap(hc50_cap)
    if data.shape[0] < 3:
        df["HC50_rank_score"] = np.nan
        df["HC50_rank_uM"] = np.nan
        df["HC50_rank_resid_log10"] = np.nan
        return df, None

    feature_series: Dict[str, pd.Series] = {}
    feature_series['delta_g'] = data["DeltaG_ads_kJ_per_mol"]
    feature_series['log_kp'] = np.log10(data["Kp_eff_nm"].astype(float).clip(lower=1e-12))
    if "Gamma_star_model" in data.columns:
        feature_series['log_gamma_model'] = np.log10(data["Gamma_star_model"].astype(float).clip(lower=1e-12))
    if "LP_star_model" in data.columns:
        feature_series['log_lp_model'] = np.log10(data["LP_star_model"].astype(float).clip(lower=1e-12))
    if "theta_at_1uM" in data.columns:
        feature_series['log_theta'] = np.log10(data["theta_at_1uM"].astype(float).clip(lower=1e-12))
    if "sequence_length" in data.columns:
        feature_series['log_seq_len'] = np.log10(data["sequence_length"].astype(float).clip(lower=1e-12))
    feature_series['qc_frac'] = data["QC_pass_fraction"]
    feature_series['ess_mean'] = data["ESS_mean_median"]

    try:
        X_train, stats_arr, feature_names = _prepare_feature_matrix(data, feature_series)
    except ValueError:
        df["HC50_rank_score"] = np.nan
        df["HC50_rank_uM"] = np.nan
        df["HC50_rank_resid_log10"] = np.nan
        return df, None

    exp_clip_train = _censor_exp(data["HC50_exp_uM"], cap)
    y_log = np.log10(np.clip(exp_clip_train, 1e-12, None))
    weights = _optimize_spearman(X_train, y_log)
    score_train = X_train @ weights
    score_sorted, log_sorted = _isotonic_regression(score_train, y_log)

    def predict_scores(df_in: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        X_all = _build_feature_matrix(df_in, feature_names, stats_arr)
        scores = X_all @ weights
        preds_log = np.interp(
            scores,
            score_sorted,
            log_sorted,
            left=log_sorted[0],
            right=log_sorted[-1],
        )
        return scores, preds_log

    scores_all, preds_log_all = predict_scores(df)

    df["HC50_rank_score"] = scores_all
    preds_linear = 10 ** preds_log_all
    exp_raw_all = _to_numeric_array(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]))
    high_mask_all = np.isfinite(exp_raw_all) & (exp_raw_all > cap)
    preds_linear = np.minimum(preds_linear, cap)
    preds_linear[high_mask_all] = cap
    df["HC50_rank_uM"] = preds_linear
    preds_log_all = np.log10(np.clip(preds_linear, 1e-12, None))
    df["HC50_rank_resid_log10"] = np.nan
    if mask.any():
        exp_clip_mask, pred_clip_mask = _censor_pair(
            df.loc[mask, "HC50_exp_uM"],
            df.loc[mask, "HC50_rank_uM"],
            cap,
        )
        df.loc[mask, "HC50_rank_resid_log10"] = (
            np.log10(np.clip(pred_clip_mask, 1e-12, None))
            - np.log10(np.clip(exp_clip_mask, 1e-12, None))
        )

    exp_clip_metrics, pred_clip_metrics = _censor_pair(
        df.loc[mask, "HC50_exp_uM"],
        df.loc[mask, "HC50_rank_uM"],
        cap,
    )
    log_exp_series = np.log10(np.clip(exp_clip_metrics, 1e-12, None))
    log_rank_series = np.log10(np.clip(pred_clip_metrics, 1e-12, None))
    metrics_rank = calc_metrics(log_exp_series, log_rank_series)

    # Leave-one-out Spearman & RMSE
    exp_raw_train = _to_numeric_array(data.get("HC50_exp_uM_raw", data["HC50_exp_uM"]))
    loo_preds = []
    if data.shape[0] >= 4:
        for idx in range(data.shape[0]):
            train_mask = np.ones(data.shape[0], dtype=bool)
            train_mask[idx] = False
            X_sub = X_train[train_mask]
            y_sub = y_log[train_mask]
            try:
                w_sub = _optimize_spearman(X_sub, y_sub)
            except Exception:
                continue
            score_sub = X_sub @ w_sub
            score_sorted_sub, log_sorted_sub = _isotonic_regression(score_sub, y_sub)
            score_val = X_train[idx:idx + 1] @ w_sub
            score_scalar = float(np.ravel(score_val)[0])
            pred_val = float(np.interp(
                score_scalar,
                score_sorted_sub,
                log_sorted_sub,
                left=log_sorted_sub[0],
                right=log_sorted_sub[-1],
            ))
            pred_val = min(pred_val, math.log10(cap))
            raw_target = exp_raw_train[idx] if idx < exp_raw_train.size else np.nan
            if np.isfinite(raw_target) and raw_target > cap:
                pred_val = math.log10(cap)
            loo_preds.append(pred_val)
        if len(loo_preds) == data.shape[0]:
            loo_rho = stats.spearmanr(y_log, loo_preds).correlation
            loo_rmse = float(np.sqrt(np.mean((np.array(loo_preds) - y_log) ** 2)))
        else:
            loo_rho = None
            loo_rmse = None
    else:
        loo_rho = None
        loo_rmse = None

    result = RankModelResults(
        feature_names=feature_names,
        weights=weights,
        feature_means=stats_arr[0],
        feature_stds=stats_arr[1],
        score_sorted=score_sorted,
        log_sorted=log_sorted,
        spearman=float(metrics_rank.get('spearman_r', np.nan)),
        pearson=float(metrics_rank.get('pearson_r', np.nan)),
        rmse_log10=float(metrics_rank.get('rmse_log10', np.nan)),
        mae_log10=float(metrics_rank.get('mae_log10', np.nan)),
        loo_spearman=float(loo_rho) if loo_rho is not None else None,
        loo_rmse_log10=loo_rmse,
        n_points=int(data.shape[0]),
    )
    return df, result


def fit_rank_model_enhanced(agg_df: pd.DataFrame, hc50_cap: Optional[float] = None) -> Tuple[pd.DataFrame, Optional[RankModelEnhancedResults]]:
    df = agg_df.copy()
    mask = df["HC50_exp_uM"].notna() & (df["HC50_exp_uM"].astype(float) > 0)
    cap = _resolve_hc50_cap(hc50_cap)
    if mask.sum() < 3:
        df["HC50_rank_enh_score"] = np.nan
        df["HC50_rank_enh_uM"] = np.nan
        df["HC50_rank_enh_resid_log10"] = np.nan
        return df, None

    features_all = _build_rank_features_enhanced(df)
    features_train = features_all.loc[mask]
    exp_raw_series = df.loc[mask, "HC50_exp_uM_raw"] if "HC50_exp_uM_raw" in df.columns else df.loc[mask, "HC50_exp_uM"]
    exp_raw_series = pd.to_numeric(exp_raw_series, errors='coerce')
    exp_raw_train = exp_raw_series
    exp_raw_train = exp_raw_series
    exp_clip_series = exp_raw_series.clip(upper=cap)
    y_log = np.log10(np.clip(exp_clip_series.to_numpy(dtype=float), 1e-12, None))
    y_log = pd.Series(y_log, index=features_train.index)
    exp_raw_train = exp_raw_series

    candidate_alphas = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    best = None

    for alpha in candidate_alphas:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(features_train)
        model = Ridge(alpha=alpha)
        model.fit(X_train, y_log)
        pred_train = model.predict(X_train)
        pred_train = np.minimum(pred_train, math.log10(cap))
        if exp_raw_train.size == pred_train.size:
            high_mask_train = np.isfinite(exp_raw_train) & (exp_raw_train > cap)
            pred_train = pred_train.copy()
            pred_train[high_mask_train] = math.log10(cap)
        spearman_val = stats.spearmanr(y_log, pred_train).correlation
        pearson_val = stats.pearsonr(y_log, pred_train).statistic if mask.sum() >= 2 else float("nan")
        resid = pred_train - y_log
        rmse_val = float(np.sqrt(np.mean(resid ** 2)))
        mae_val = float(np.mean(np.abs(resid)))

        loo_rho = None
        loo_rmse = None
        if mask.sum() >= 5:
            loo_preds = []
            train_index = features_train.index
            for idx in train_index:
                others = train_index.difference([idx])
                scaler_i = StandardScaler().fit(features_train.loc[others])
                model_i = Ridge(alpha=alpha)
                model_i.fit(
                    scaler_i.transform(features_train.loc[others]),
                    y_log.loc[others]
                )
                loo_pred = model_i.predict(scaler_i.transform(features_train.loc[[idx]]))[0]
                loo_pred = min(loo_pred, math.log10(cap))
                raw_target = exp_raw_train.get(idx, np.nan)
                if pd.notna(raw_target) and raw_target > cap:
                    loo_pred = math.log10(cap)
                loo_preds.append(loo_pred)
            if loo_preds:
                loo_rho = stats.spearmanr(y_log, loo_preds).correlation
                loo_rmse = float(np.sqrt(np.mean((np.array(loo_preds) - y_log.values) ** 2)))

        score = (loo_rho if loo_rho is not None and np.isfinite(loo_rho) else -np.inf)
        if best is None or score > best['score'] or (np.isclose(score, best['score']) and rmse_val < best['rmse']):
            best = {
                'alpha': alpha,
                'scaler': scaler,
                'model': model,
                'pred_train': pred_train,
                'spearman': spearman_val,
                'pearson': pearson_val,
                'rmse': rmse_val,
                'mae': mae_val,
                'loo_rho': loo_rho,
                'loo_rmse': loo_rmse,
                'score': score,
            }

    if best is None:
        df["HC50_rank_enh_score"] = np.nan
        df["HC50_rank_enh_uM"] = np.nan
        df["HC50_rank_enh_resid_log10"] = np.nan
        return df, None

    scaler = best['scaler']
    model = best['model']
    pred_log_train = best['pred_train']
    pred_log_all = model.predict(scaler.transform(features_all))
    pred_log_all = np.minimum(pred_log_all, math.log10(cap))
    exp_raw_all = _to_numeric_array(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]))
    high_mask_all = np.isfinite(exp_raw_all) & (exp_raw_all > cap)
    pred_log_all = pred_log_all.astype(float)
    pred_log_all[high_mask_all] = math.log10(cap)
    df["HC50_rank_enh_uM"] = np.power(10.0, pred_log_all)
    df["HC50_rank_enh_uM"] = np.minimum(df["HC50_rank_enh_uM"].astype(float), cap)
    df.loc[high_mask_all, "HC50_rank_enh_uM"] = cap
    df["HC50_rank_enh_score"] = pred_log_all
    df["HC50_rank_enh_resid_log10"] = np.nan
    df.loc[mask, "HC50_rank_enh_resid_log10"] = pred_log_train - y_log

    scale_safe = np.where(scaler.scale_ == 0, 1.0, scaler.scale_)
    coefficients = model.coef_ / scale_safe
    intercept = float(model.intercept_ - np.sum((model.coef_ * scaler.mean_) / scale_safe))

    result = RankModelEnhancedResults(
        feature_names=tuple(features_all.columns),
        coefficients=coefficients,
        intercept=intercept,
        feature_means=scaler.mean_,
        feature_stds=scaler.scale_,
        alpha=float(best['alpha']),
        spearman=float(best['spearman']),
        pearson=float(best['pearson']),
        rmse_log10=float(best['rmse']),
        mae_log10=float(best['mae']),
        loo_spearman=best['loo_rho'],
        loo_rmse_log10=best['loo_rmse'],
        n_points=int(mask.sum()),
    )
    return df, result


def fit_rank_model_knn(agg_df: pd.DataFrame, hc50_cap: Optional[float] = None) -> Tuple[pd.DataFrame, Optional[RankModelKNNResults]]:
    df = agg_df.copy()
    mask = df["HC50_exp_uM"].notna() & (df["HC50_exp_uM"].astype(float) > 0)
    cap = _resolve_hc50_cap(hc50_cap)
    if mask.sum() < 3:
        df["HC50_rank_knn_score"] = np.nan
        df["HC50_rank_knn_uM"] = np.nan
        df["HC50_rank_knn_resid_log10"] = np.nan
        return df, None

    features_all = _build_rank_features_enhanced(df)
    features_train = features_all.loc[mask]
    exp_raw_series = df.loc[mask, "HC50_exp_uM_raw"] if "HC50_exp_uM_raw" in df.columns else df.loc[mask, "HC50_exp_uM"]
    exp_raw_series = pd.to_numeric(exp_raw_series, errors='coerce')
    exp_raw_train = exp_raw_series
    exp_clip_series = exp_raw_series.clip(upper=cap)
    y_log = np.log10(np.clip(exp_clip_series.to_numpy(dtype=float), 1e-12, None))
    y_log = pd.Series(y_log, index=features_train.index)

    candidate_k = [3, 4, 5, 6, 7]
    best = None
    train_index = features_train.index

    for k in candidate_k:
        if k >= len(y_log):
            continue
        loo_preds = []
        for idx in train_index:
            others = train_index.difference([idx])
            X_train = features_train.loc[others]
            y_train = y_log.loc[others]
            scaler = StandardScaler().fit(X_train)
            X_train_scaled = scaler.transform(X_train)
            n_neighbors = min(k, X_train_scaled.shape[0])
            model = KNeighborsRegressor(n_neighbors=n_neighbors, weights='distance')
            model.fit(X_train_scaled, y_train)
            x_test = scaler.transform(features_train.loc[[idx]])
            pred_val = model.predict(x_test)[0]
            pred_val = min(pred_val, math.log10(cap))
            raw_target = exp_raw_series.get(idx, np.nan)
            if pd.notna(raw_target) and raw_target > cap:
                pred_val = math.log10(cap)
            loo_preds.append(pred_val)
        if not loo_preds:
            continue
        loo_preds = np.array(loo_preds)
        spearman_val = stats.spearmanr(y_log, loo_preds).correlation
        rmse_val = float(np.sqrt(np.mean((loo_preds - y_log.values) ** 2)))
        if best is None or (spearman_val is not None and np.isfinite(spearman_val) and spearman_val > best['spearman']) or (
            (spearman_val == best['spearman']) and rmse_val < best['rmse']
        ):
            best = {
                'k': k,
                'spearman': spearman_val if spearman_val is not None else -np.inf,
                'rmse': rmse_val,
            }

    if best is None:
        df["HC50_rank_knn_score"] = np.nan
        df["HC50_rank_knn_uM"] = np.nan
        df["HC50_rank_knn_resid_log10"] = np.nan
        return df, None

    k_best = best['k']
    scaler = StandardScaler().fit(features_train)
    X_train_scaled = scaler.transform(features_train)
    n_neighbors = min(k_best, X_train_scaled.shape[0])
    model = KNeighborsRegressor(n_neighbors=n_neighbors, weights='distance')
    model.fit(X_train_scaled, y_log)
    pred_train_values = model.predict(X_train_scaled)
    pred_train_values = np.minimum(pred_train_values, math.log10(cap))
    high_mask_train = exp_raw_train > cap
    pred_train = pd.Series(pred_train_values, index=features_train.index)
    pred_train.loc[high_mask_train.fillna(False)] = math.log10(cap)
    pearson_val = stats.pearsonr(y_log, pred_train).statistic if mask.sum() >= 2 else float("nan")
    spearman_val = stats.spearmanr(y_log, pred_train).correlation
    resid = pred_train.values - y_log.values
    rmse_val = float(np.sqrt(np.mean(resid ** 2)))
    mae_val = float(np.mean(np.abs(resid)))

    pred_all = model.predict(scaler.transform(features_all))
    pred_all = np.minimum(pred_all, math.log10(cap))
    exp_raw_all_series = pd.to_numeric(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]), errors='coerce')
    high_mask_all = exp_raw_all_series > cap
    pred_all = pred_all.astype(float)
    high_mask_all_bool = high_mask_all.fillna(False)
    pred_all[high_mask_all_bool.to_numpy()] = math.log10(cap)
    df["HC50_rank_knn_uM"] = np.power(10.0, pred_all)
    df["HC50_rank_knn_uM"] = np.minimum(df["HC50_rank_knn_uM"].astype(float), cap)
    df.loc[high_mask_all_bool, "HC50_rank_knn_uM"] = cap
    df["HC50_rank_knn_score"] = pred_all
    df["HC50_rank_knn_resid_log10"] = np.nan
    df.loc[mask, "HC50_rank_knn_resid_log10"] = pred_train.values - y_log.values

    # LOO metrics already computed in best selection for spearman; recompute rmse for accuracy
    loo_preds = []
    for idx in train_index:
        others = train_index.difference([idx])
        scaler_i = StandardScaler().fit(features_train.loc[others])
        X_train_i = scaler_i.transform(features_train.loc[others])
        model_i = KNeighborsRegressor(n_neighbors=min(k_best, X_train_i.shape[0]), weights='distance')
        model_i.fit(X_train_i, y_log.loc[others])
        x_test = scaler_i.transform(features_train.loc[[idx]])
        pred_val = model_i.predict(x_test)[0]
        pred_val = min(pred_val, math.log10(cap))
        raw_target = exp_raw_train.get(idx, np.nan)
        if np.isfinite(raw_target) and raw_target > cap:
            pred_val = math.log10(cap)
        loo_preds.append(pred_val)
    loo_preds = np.array(loo_preds)
    loo_rho = stats.spearmanr(y_log, loo_preds).correlation
    loo_rmse = float(np.sqrt(np.mean((loo_preds - y_log.values) ** 2)))

    result = RankModelKNNResults(
        feature_names=tuple(features_all.columns),
        k_neighbors=int(k_best),
        feature_means=scaler.mean_,
        feature_stds=scaler.scale_,
        spearman=float(spearman_val),
        pearson=float(pearson_val),
        rmse_log10=rmse_val,
        mae_log10=mae_val,
        loo_spearman=float(loo_rho),
        loo_rmse_log10=loo_rmse,
        n_points=int(mask.sum()),
    )
    return df, result


def fit_rank_model_ensemble(agg_df: pd.DataFrame, hc50_cap: Optional[float] = None) -> Tuple[pd.DataFrame, Optional[RankModelEnsembleResults]]:
    df = agg_df.copy()
    mask = df["HC50_exp_uM"].notna() & (df["HC50_exp_uM"].astype(float) > 0)
    cap = _resolve_hc50_cap(hc50_cap)
    if mask.sum() < 3:
        df["HC50_rank_combo_score"] = np.nan
        df["HC50_rank_combo_uM"] = np.nan
        df["HC50_rank_combo_resid_log10"] = np.nan
        return df, None

    components_all: Dict[str, pd.Series] = {
        'rank_basic': np.log10(np.clip(df.get("HC50_rank_uM"), 1e-6, None).astype(float)),
        'rank_enhanced': np.log10(np.clip(df.get("HC50_rank_enh_uM"), 1e-6, None).astype(float)) if "HC50_rank_enh_uM" in df.columns else None,
    }
    component_names = [name for name, series in components_all.items() if series is not None]
    if len(component_names) < 1:
        df["HC50_rank_combo_score"] = np.nan
        df["HC50_rank_combo_uM"] = np.nan
        df["HC50_rank_combo_resid_log10"] = np.nan
        return df, None

    X = np.column_stack([components_all[name].loc[mask] for name in component_names])
    exp_raw_series = df.loc[mask, "HC50_exp_uM_raw"] if "HC50_exp_uM_raw" in df.columns else df.loc[mask, "HC50_exp_uM"]
    exp_raw_series = pd.to_numeric(exp_raw_series, errors='coerce')
    exp_raw_train = exp_raw_series
    exp_clip_train = exp_raw_series.clip(upper=cap)
    y_log = np.log10(np.clip(exp_clip_train.to_numpy(dtype=float), 1e-12, None))
    y_log = pd.Series(y_log, index=df.index[mask])

    step = 0.05
    weights_grid = []
    if len(component_names) == 1:
        weights_grid = [np.array([1.0])]
    else:
        n = len(component_names)
        for raw in np.ndindex(*([int(1/step)+1] * n)):
            weights = np.array(raw, dtype=float) * step
            if np.isclose(weights.sum(), 1.0):
                weights_grid.append(weights)
    if not weights_grid:
        weights_grid = [np.ones(len(component_names)) / len(component_names)]

    best = None
    for weights in weights_grid:
        preds = X @ weights
        rho = stats.spearmanr(y_log.values, preds).correlation
        rmse = float(np.sqrt(np.mean((preds - y_log.values) ** 2)))
        if best is None or (rho is not None and np.isfinite(rho) and rho > best['rho']) or (
            (rho == best['rho']) and rmse < best['rmse']
        ):
            best = {
                'weights': weights,
                'rho': rho if rho is not None else -np.inf,
                'rmse': rmse,
                'preds': preds,
            }

    weights = best['weights'] if best is not None else np.ones(len(component_names)) / len(component_names)
    pred_train = X @ weights
    pred_train = np.minimum(pred_train, math.log10(cap))
    high_mask_train = exp_raw_train > cap
    pred_train = pred_train.astype(float)
    pred_train[high_mask_train.fillna(False).to_numpy()] = math.log10(cap)
    resid = pred_train - y_log.values
    rmse_val = float(np.sqrt(np.mean(resid ** 2)))
    mae_val = float(np.mean(np.abs(resid)))
    spearman_val = stats.spearmanr(y_log.values, pred_train).correlation
    pearson_val = stats.pearsonr(y_log.values, pred_train).statistic if mask.sum() >= 2 else float("nan")

    X_all = np.column_stack([components_all[name] for name in component_names])
    pred_all = X_all @ weights
    pred_all = np.minimum(pred_all, math.log10(cap))
    exp_raw_all = pd.to_numeric(df.get("HC50_exp_uM_raw", df["HC50_exp_uM"]), errors='coerce')
    high_mask_all = exp_raw_all > cap
    pred_all = pred_all.astype(float)
    pred_all[high_mask_all.fillna(False).to_numpy()] = math.log10(cap)
    df["HC50_rank_combo_score"] = pred_all
    df["HC50_rank_combo_uM"] = np.power(10.0, pred_all)
    df["HC50_rank_combo_uM"] = np.minimum(df["HC50_rank_combo_uM"].astype(float), cap)
    df.loc[high_mask_all.fillna(False), "HC50_rank_combo_uM"] = cap
    df["HC50_rank_combo_resid_log10"] = np.nan
    df.loc[mask, "HC50_rank_combo_resid_log10"] = pred_train - y_log

    # Leave-one-out using same combination weights
    loo_preds = []
    for idx in range(len(y_log)):
        mask_rest = np.ones(len(y_log), dtype=bool)
        mask_rest[idx] = False
        preds_rest = pred_train[mask_rest]
        y_rest = y_log.iloc[mask_rest]
        if preds_rest.size >= 2:
            rho_rest = stats.spearmanr(y_rest, preds_rest).correlation
        else:
            rho_rest = np.nan
        loo_preds.append((pred_train[idx], rho_rest))
    loo_rho = np.nanmean([r for _, r in loo_preds if r is not None])
    loo_rmse = float(np.sqrt(np.mean((pred_train - y_log.values) ** 2)))

    result = RankModelEnsembleResults(
        component_names=tuple(component_names),
        coefficients=weights,
        intercept=0.0,
        spearman=float(spearman_val),
        pearson=float(pearson_val),
        rmse_log10=rmse_val,
        mae_log10=mae_val,
        loo_spearman=float(loo_rho) if loo_rho is not None else None,
        loo_rmse_log10=loo_rmse,
        n_points=int(mask.sum()),
    )
    return df, result


def fit_rank_model_v3(agg_df: pd.DataFrame, hc50_cap: float) -> Tuple[pd.DataFrame, Optional[RankModelV3Results]]:
    df = agg_df.copy()
    cap = float(hc50_cap)
    mask = df["HC50_exp_uM"].notna() & (df["HC50_exp_uM"].astype(float) > 0)
    if mask.sum() < 4:
        df["HC50_rank_v3_score"] = np.nan
        df["HC50_rank_v3_uM"] = np.nan
        df["HC50_rank_v3_resid_log10"] = np.nan
        return df, None

    features_all = _prepare_rank_v3_features(df)
    features_all = features_all.astype(float)
    feature_means = features_all.mean(axis=0, skipna=True)
    features_all = features_all.fillna(feature_means)
    X_all = features_all.loc[mask]
    feature_names = tuple(X_all.columns)
    X = X_all.to_numpy(dtype=float)

    exp_raw = pd.to_numeric(df.loc[mask, "HC50_exp_uM"], errors='coerce')
    exp_raw = exp_raw.fillna(cap)
    exp_clip = exp_raw.clip(upper=cap)
    y_log = np.log10(np.clip(exp_clip.to_numpy(dtype=float), 1e-12, None))
    groups = df.loc[mask, "Peptide_ID"].astype(str)

    ess_weights = pd.to_numeric(df.loc[mask, "ESS_mean_median"], errors='coerce').fillna(0.0)
    qc_weights = pd.to_numeric(df.get("QC_pass_fraction", pd.Series(np.ones(len(df)))), errors='coerce').fillna(1.0)
    sample_weights = np.clip(ess_weights.values, 1.0, None) * np.clip(qc_weights.loc[mask].values, 0.25, 1.0)

    n_groups = max(1, groups.nunique())
    outer_splits = min(4, n_groups) if n_groups >= 2 else 2
    if n_groups >= 3:
        outer_splits = max(3, outer_splits)
    outer_splits = max(2, outer_splits)

    param_grid = []
    alphas = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    l1_ratios = [0.05, 0.1, 0.2, 0.4]
    for a in alphas:
        for l1 in l1_ratios:
            param_grid.append((a, l1))

    outer_cv = GroupKFold(n_splits=outer_splits) if n_groups >= outer_splits else KFold(n_splits=outer_splits, shuffle=True, random_state=42)
    param_scores: Dict[Tuple[float, float], list] = {p: [] for p in param_grid}
    cv_metrics = []

    for outer_train, outer_test in outer_cv.split(X, y_log, groups if isinstance(outer_cv, GroupKFold) else None):
        X_train, X_test = X[outer_train], X[outer_test]
        y_train, y_test = y_log[outer_train], y_log[outer_test]
        groups_train = groups.iloc[outer_train]
        weights_train = sample_weights[outer_train]
        inner_groups_unique = groups_train.nunique()
        inner_splits = min(3, inner_groups_unique) if inner_groups_unique >= 2 else min(3, len(outer_train))
        inner_splits = max(2, inner_splits)
        inner_cv = GroupKFold(n_splits=inner_splits) if inner_groups_unique >= inner_splits else KFold(n_splits=inner_splits, shuffle=True, random_state=42)

        best_param = None
        best_score = -np.inf

        for alpha, l1_ratio in param_grid:
            inner_scores = []
            for inner_train, inner_val in inner_cv.split(X_train, y_train, groups_train if isinstance(inner_cv, GroupKFold) else None):
                X_inner_train, X_inner_val = X_train[inner_train], X_train[inner_val]
                y_inner_train, y_inner_val = y_train[inner_train], y_train[inner_val]
                w_inner = weights_train[inner_train]
                pipeline = Pipeline([
                    ("scaler", StandardScaler()),
                    ("model", ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=10000, fit_intercept=True))
                ])
                try:
                    pipeline.fit(X_inner_train, y_inner_train, model__sample_weight=w_inner)
                    preds_val = pipeline.predict(X_inner_val)
                    preds_val = np.clip(preds_val, a_min=None, a_max=math.log10(cap))
                    if (
                        preds_val.size == 0
                        or y_inner_val.size == 0
                        or np.allclose(preds_val, preds_val[0])
                        or np.allclose(y_inner_val, y_inner_val[0])
                    ):
                        spearman_inner = np.nan
                        rmse_inner = np.nan
                    else:
                        spearman_inner = stats.spearmanr(y_inner_val, preds_val).correlation
                        rmse_inner = float(np.sqrt(np.mean((preds_val - y_inner_val) ** 2)))
                except Exception:
                    spearman_inner = np.nan
                    rmse_inner = np.nan
                if np.isnan(spearman_inner) or np.isnan(rmse_inner):
                    continue
                score_inner = spearman_inner - 0.15 * rmse_inner
                inner_scores.append(score_inner)
            if inner_scores:
                mean_score = float(np.mean(inner_scores))
                param_scores[(alpha, l1_ratio)].append(mean_score)
                if mean_score > best_score:
                    best_score = mean_score
                    best_param = (alpha, l1_ratio)

        if best_param is None:
            best_param = min(param_grid, key=lambda p: (p[0], p[1]))

        alpha_opt, l1_opt = best_param
        pipeline = Pipeline([
            ("scaler", StandardScaler()),
            ("model", ElasticNet(alpha=alpha_opt, l1_ratio=l1_opt, max_iter=10000, fit_intercept=True))
        ])
        pipeline.fit(X_train, y_train, model__sample_weight=weights_train)
        preds_train = pipeline.predict(X_train)
        iso = IsotonicRegression(out_of_bounds="clip").fit(np.clip(preds_train, None, math.log10(cap)), y_train)
        preds_test_raw = pipeline.predict(X_test)
        preds_test_iso = iso.predict(np.clip(preds_test_raw, None, math.log10(cap)))
        metrics_outer = calc_metrics(y_test, np.clip(preds_test_iso, -np.inf, math.log10(cap)))
        cv_metrics.append(metrics_outer)

    mean_param_scores = {param: (np.mean(scores) if scores else -np.inf) for param, scores in param_scores.items()}
    best_global = max(mean_param_scores.items(), key=lambda kv: kv[1])[0]
    alpha_best, l1_best = best_global

    final_pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("model", ElasticNet(alpha=alpha_best, l1_ratio=l1_best, max_iter=10000, fit_intercept=True))
    ])
    final_pipeline.fit(X, y_log, model__sample_weight=sample_weights)
    preds_train_full = final_pipeline.predict(X)
    iso_final = IsotonicRegression(out_of_bounds="clip")
    iso_final.fit(np.clip(preds_train_full, None, math.log10(cap)), y_log)

    preds_all_raw = np.full(df.shape[0], np.nan)
    preds_all_iso = np.full(df.shape[0], np.nan)
    preds_all_raw[mask.to_numpy()] = final_pipeline.predict(X)
    preds_all_iso[mask.to_numpy()] = iso_final.predict(np.clip(preds_all_raw[mask.to_numpy()], None, math.log10(cap)))
    preds_all_iso = np.clip(preds_all_iso, -np.inf, math.log10(cap))
    preds_all_uM = np.power(10.0, np.minimum(preds_all_iso, math.log10(cap)))

    df["HC50_rank_v3_score"] = preds_all_raw
    df["HC50_rank_v3_uM"] = preds_all_uM
    df["HC50_rank_v3_resid_log10"] = np.nan
    if mask.any():
        y_mask = np.log10(np.clip(df.loc[mask, "HC50_exp_uM"].astype(float), 1e-12, cap))
        df.loc[mask, "HC50_rank_v3_resid_log10"] = preds_all_iso[mask.to_numpy()] - y_mask

    metrics_full = calc_metrics(
        np.log10(np.clip(df.loc[mask, "HC50_exp_uM"].astype(float), 1e-12, cap)),
        preds_all_iso[mask.to_numpy()]
    )

    if cv_metrics:
        cv_spearman = float(np.nanmean([m.get('spearman_r') for m in cv_metrics]))
        cv_rmse = float(np.nanmean([m.get('rmse_log10') for m in cv_metrics]))
        cv_mae = float(np.nanmean([m.get('mae_log10') for m in cv_metrics]))
    else:
        cv_spearman = cv_rmse = cv_mae = None

    scaler = final_pipeline.named_steps['scaler']
    model = final_pipeline.named_steps['model']
    scale_safe = np.where(scaler.scale_ == 0, 1.0, scaler.scale_)
    coefficients = model.coef_ / scale_safe
    intercept = float(model.intercept_ - np.sum((model.coef_ * scaler.mean_) / scale_safe))

    result = RankModelV3Results(
        feature_names=feature_names,
        coefficients=coefficients,
        intercept=intercept,
        alpha=float(alpha_best),
        l1_ratio=float(l1_best),
        scaler_mean=scaler.mean_.copy(),
        scaler_scale=scaler.scale_.copy(),
        iso_x=iso_final.X_thresholds_.copy(),
        iso_y=iso_final.y_thresholds_.copy(),
        spearman=float(metrics_full.get('spearman_r', np.nan)),
        pearson=float(metrics_full.get('pearson_r', np.nan)),
        rmse_log10=float(metrics_full.get('rmse_log10', np.nan)),
        mae_log10=float(metrics_full.get('mae_log10', np.nan)),
        cv_spearman=cv_spearman,
        cv_rmse_log10=cv_rmse,
        cv_mae_log10=cv_mae,
        n_points=int(mask.sum()),
        sample_weight_mean=float(np.mean(sample_weights))
    )
    return df, result


def compute_lp_stats(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    lp_source = "LP_star_model" if df["LP_star_model"].notna().sum() >= 3 else "LP_star_obs"
    lp_vals = df[lp_source].dropna().astype(float)
    stats_dict: Dict[str, Any]
    if lp_vals.empty:
        stats_dict = {
            "n": 0,
            "median": float("nan"),
            "q10": float("nan"),
            "q25": float("nan"),
            "q75": float("nan"),
            "q90": float("nan"),
            "iqr": float("nan"),
            "source": lp_source,
        }
    else:
        stats_dict = {
            "n": int(lp_vals.shape[0]),
            "median": float(lp_vals.median()),
            "q10": float(lp_vals.quantile(0.10)),
            "q25": float(lp_vals.quantile(0.25)),
            "q75": float(lp_vals.quantile(0.75)),
            "q90": float(lp_vals.quantile(0.90)),
            "iqr": float(lp_vals.quantile(0.75) - lp_vals.quantile(0.25)),
            "source": lp_source,
        }
    df["LP_star_Q25"] = stats_dict["q25"]
    df["LP_star_Q50"] = stats_dict["median"]
    df["LP_star_Q75"] = stats_dict["q75"]
    return df, stats_dict


def make_plots(
    agg_df: pd.DataFrame,
    cal_results: Optional[CalibrationResults],
    gamma_results: Optional[GammaModelResults],
    powerlaw: Optional[PowerLawResults],
    rank_results: Optional[RankModelResults],
    rank_results_enh: Optional[RankModelEnhancedResults],
    rank_results_knn: Optional[RankModelKNNResults],
    rank_results_combo: Optional[RankModelEnsembleResults],
    rank_results_v3: Optional[RankModelV3Results],
    plots_dir: Path,
    hc50_cap: float,
):
    plots_dir.mkdir(parents=True, exist_ok=True)
    generated: set[Path] = set()
    cap = float(hc50_cap)

    mask_exp = (
        agg_df["HC50_exp_uM"].notna()
        & agg_df["HC50_exp_uM"].astype(float) > 0
    )
    if mask_exp.sum() >= 1:
        exp_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_exp_uM"], cap)
        log_exp = np.log10(np.clip(exp_vals, 1e-12, None))
        parity_path = plots_dir / "parity_log10.png"
        fig, ax = plt.subplots(figsize=(8.6, 6.2))
        model_logs = []
        if rank_results and mask_exp.any():
            rank_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_rank_uM"], cap)
            log_rank = np.log10(np.clip(rank_vals, 1e-12, None))
            ax.scatter(
                log_exp,
                log_rank,
                marker="D",
                color="tab:green",
                alpha=0.8,
                label="Rank-Modell",
            )
            model_logs.append(log_rank)
        if rank_results_enh and mask_exp.any():
            rank_enh_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_rank_enh_uM"], cap)
            log_rank_enh = np.log10(np.clip(rank_enh_vals, 1e-12, None))
            ax.scatter(
                log_exp,
                log_rank_enh,
                marker="^",
                color="tab:purple",
                alpha=0.8,
                label="Rank-Modell (enhanced)",
            )
            model_logs.append(log_rank_enh)
        if rank_results_knn and mask_exp.any():
            rank_knn_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_rank_knn_uM"], cap)
            log_rank_knn = np.log10(np.clip(rank_knn_vals, 1e-12, None))
            ax.scatter(
                log_exp,
                log_rank_knn,
                marker="s",
                color="tab:cyan",
                alpha=0.8,
                label=f"Rank-Modell (kNN, k={rank_results_knn.k_neighbors})",
            )
            model_logs.append(log_rank_knn)
        if rank_results_combo and mask_exp.any():
            rank_combo_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_rank_combo_uM"], cap)
            log_rank_combo = np.log10(np.clip(rank_combo_vals, 1e-12, None))
            ax.scatter(
                log_exp,
                log_rank_combo,
                marker="P",
                color="tab:orange",
                alpha=0.8,
                label="Rank-Modell (Ensemble)",
            )
            model_logs.append(log_rank_combo)
        if rank_results_v3 and mask_exp.any():
            rank_v3_vals = _censor_exp(agg_df.loc[mask_exp, "HC50_rank_v3_uM"], cap)
            log_rank_v3 = np.log10(np.clip(rank_v3_vals, 1e-12, None))
            ax.scatter(
                log_exp,
                log_rank_v3,
                marker="X",
                color="tab:red",
                alpha=0.85,
                label="Rank-Modell v3",
            )
            model_logs.append(log_rank_v3)

        all_vals = [log_exp]
        if model_logs:
            all_vals.extend(model_logs)
        combined = np.concatenate(all_vals)
        lim_min = math.floor(combined.min() - 0.2)
        lim_max = math.ceil(combined.max() + 0.2)
        ax.plot([lim_min, lim_max], [lim_min, lim_max], linestyle="--", color="grey", linewidth=1.0)
        ax.set_xlim(lim_min, lim_max)
        ax.set_ylim(lim_min, lim_max)
        ax.set_xlabel(r"log$_{10}$ HC50$_{exp}$ [µM]")
        ax.set_ylabel(r"log$_{10}$ HC50 [µM]")
        ax.set_title("HC50-Paritätsdiagramm (aggregiert)")
        legend_entries = []
        if rank_results:
            legend_entries.append(
                f"Rank RMSE={rank_results.rmse_log10:.3f}, ρ={rank_results.spearman:.3f}"
            )
        if rank_results_enh:
            legend_entries.append(
                f"Rank(enh) RMSE={rank_results_enh.rmse_log10:.3f}, ρ={rank_results_enh.spearman:.3f}"
            )
        if rank_results_knn:
            legend_entries.append(
                f"Rank(kNN) RMSE={rank_results_knn.rmse_log10:.3f}, ρ={rank_results_knn.spearman:.3f}"
            )
        if rank_results_combo:
            legend_entries.append(
                f"Rank(Ensemble) RMSE={rank_results_combo.rmse_log10:.3f}, ρ={rank_results_combo.spearman:.3f}"
            )
        if rank_results_v3:
            legend_entries.append(
                f"Rank(v3) RMSE={rank_results_v3.rmse_log10:.3f}, ρ={rank_results_v3.spearman:.3f}"
            )
        if legend_entries:
            ax.text(
                0.02,
                0.02,
                "\n".join(legend_entries),
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=10,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="none"),
            )
        ax.legend(loc="best")
        _format_publication_axes(ax)
        fig.tight_layout()
        fig.savefig(parity_path, dpi=300)
        plt.close(fig)
        generated.add(parity_path)

        resid_path = plots_dir / "residuals_vs_logKp.png"
        fig, ax = plt.subplots(figsize=(8.6, 5.6))
        ax.axhline(0.0, color="grey", linestyle="--", linewidth=1.0)
        log_kp = np.log10(agg_df.loc[mask_exp, "Kp_eff_nm"].astype(float))
        if rank_results:
            ax.scatter(
                log_kp,
                agg_df.loc[mask_exp, "HC50_rank_resid_log10"].astype(float),
                color="tab:green",
                marker="D",
                alpha=0.8,
                label="Rank-Modell",
            )
        if rank_results_enh:
            ax.scatter(
                log_kp,
                agg_df.loc[mask_exp, "HC50_rank_enh_resid_log10"].astype(float),
                color="tab:purple",
                marker="^",
                alpha=0.8,
                label="Rank-Modell (enhanced)",
            )
        if rank_results_knn:
            ax.scatter(
                log_kp,
                agg_df.loc[mask_exp, "HC50_rank_knn_resid_log10"].astype(float),
                color="tab:cyan",
                marker="s",
                alpha=0.8,
                label=f"Rank-Modell (kNN, k={rank_results_knn.k_neighbors})",
            )
        if rank_results_combo:
            ax.scatter(
                log_kp,
                agg_df.loc[mask_exp, "HC50_rank_combo_resid_log10"].astype(float),
                color="tab:orange",
                marker="P",
                alpha=0.8,
                label="Rank-Modell (Ensemble)",
            )
        if rank_results_v3:
            ax.scatter(
                log_kp,
                agg_df.loc[mask_exp, "HC50_rank_v3_resid_log10"].astype(float),
                color="tab:red",
                marker="X",
                alpha=0.85,
                label="Rank-Modell v3",
            )
        ax.set_xlabel(r"log$_{10}$ Kp$_{eff}$ [nm]")
        ax.set_ylabel(r"Residuum log$_{10}$ HC50")
        ax.set_title("Residuen vs. log$_{10}$ Kp$_{eff}$")
        ax.legend(loc="best")
        _format_publication_axes(ax)
        fig.tight_layout()
        fig.savefig(resid_path, dpi=300)
        plt.close(fig)
        generated.add(resid_path)

        exp_df = agg_df.loc[mask_exp].copy()
        exp_df.sort_values("HC50_exp_uM", inplace=True)
        exp_df["HC50_exp_plot_uM"] = exp_df["HC50_exp_uM"].clip(upper=cap)
        for col in [
            "HC50_rank_uM",
            "HC50_rank_enh_uM",
            "HC50_rank_knn_uM",
            "HC50_rank_combo_uM",
            "HC50_rank_v3_uM",
        ]:
            if col in exp_df.columns:
                exp_df[f"{col}_plot"] = exp_df[col].clip(upper=cap)
        x = np.arange(exp_df.shape[0])
        fig, ax = plt.subplots(figsize=(9.5, 5.8))
        ax.scatter(x, exp_df["HC50_exp_plot_uM"], color="black", marker="o", label="Experiment", zorder=3)
        if rank_results:
            ax.scatter(x, exp_df["HC50_rank_uM_plot"], color="tab:green", marker="D", label="Rank-Modell", zorder=3)
            for xi, y_exp, y_rank in zip(x, exp_df["HC50_exp_plot_uM"], exp_df["HC50_rank_uM_plot"]):
                ax.plot([xi, xi], [y_exp, y_rank], color="tab:green", lw=1.0, alpha=0.6)
        if rank_results_enh:
            ax.scatter(x, exp_df["HC50_rank_enh_uM_plot"], color="tab:purple", marker="^", alpha=0.75, label="Rank-Modell (enhanced)", zorder=3)
            for xi, y_exp, y_rank in zip(x, exp_df["HC50_exp_plot_uM"], exp_df["HC50_rank_enh_uM_plot"]):
                ax.plot([xi, xi], [y_exp, y_rank], color="tab:purple", lw=1.0, alpha=0.6)
        if rank_results_knn:
            ax.scatter(x, exp_df["HC50_rank_knn_uM_plot"], color="tab:cyan", marker="s", alpha=0.75, label=f"Rank-Modell (kNN, k={rank_results_knn.k_neighbors})", zorder=3)
            for xi, y_exp, y_rank in zip(x, exp_df["HC50_exp_plot_uM"], exp_df["HC50_rank_knn_uM_plot"]):
                ax.plot([xi, xi], [y_exp, y_rank], color="tab:cyan", lw=1.0, alpha=0.6)
        if rank_results_combo:
            ax.scatter(x, exp_df["HC50_rank_combo_uM_plot"], color="tab:orange", marker="P", alpha=0.75, label="Rank-Modell (Ensemble)", zorder=3)
            for xi, y_exp, y_rank in zip(x, exp_df["HC50_exp_plot_uM"], exp_df["HC50_rank_combo_uM_plot"]):
                ax.plot([xi, xi], [y_exp, y_rank], color="tab:orange", lw=1.0, alpha=0.6)
        if rank_results_v3:
            ax.scatter(x, exp_df["HC50_rank_v3_uM_plot"], color="tab:red", marker="X", alpha=0.85, label="Rank-Modell v3", zorder=3)
            for xi, y_exp, y_rank in zip(x, exp_df["HC50_exp_plot_uM"], exp_df["HC50_rank_v3_uM_plot"]):
                ax.plot([xi, xi], [y_exp, y_rank], color="tab:red", lw=1.0, alpha=0.6)
        error_mask = (
            exp_df["HC50_exp_CI_low_uM"].notna()
            & exp_df["HC50_exp_CI_high_uM"].notna()
        )
        if error_mask.any():
            yerr = [
                exp_df.loc[error_mask, "HC50_exp_plot_uM"] - exp_df.loc[error_mask, "HC50_exp_CI_low_uM"].clip(upper=cap),
                exp_df.loc[error_mask, "HC50_exp_CI_high_uM"].clip(upper=cap) - exp_df.loc[error_mask, "HC50_exp_plot_uM"],
            ]
            ax.errorbar(
                x[error_mask],
                exp_df.loc[error_mask, "HC50_exp_plot_uM"],
                yerr=yerr,
                fmt="none",
                ecolor="lightgray",
                elinewidth=1.5,
                capsize=4,
                alpha=0.8,
            )
        ax.set_yscale("log")
        ax.set_xlabel("Peptid (aufsteigend nach HC50 exp)")
        ax.set_ylabel("HC50 [µM] (log-Skala)")
        ax.set_title("HC50: gemessen vs. Modelle")
        ax.set_xticks(x)
        ax.set_xticklabels(exp_df["Peptide_ID"].astype(str), rotation=45, ha="right")
        ax.legend(loc="best")
        _format_publication_axes(ax)
        fig.tight_layout()
        err_path = plots_dir / "error_bars_hc50.png"
        fig.savefig(err_path, dpi=300)
        plt.close(fig)
        generated.add(err_path)

    lp_vals = agg_df["LP_star"].dropna().astype(float)
    if not lp_vals.empty:
        lp_path = plots_dir / "lp_star_distribution.png"
        fig, ax = plt.subplots(figsize=(7.8, 5.6))
        ax.hist(lp_vals, bins=min(20, len(lp_vals)), color="tab:blue", alpha=0.7, edgecolor="white")
        q25, q50, q75 = lp_vals.quantile([0.25, 0.5, 0.75])
        ax.axvline(q25, color="tab:orange", linestyle="--", label="Q25")
        ax.axvline(q50, color="tab:red", linestyle="--", label="Median")
        ax.axvline(q75, color="tab:olive", linestyle="--", label="Q75")
        ax.set_xlabel("L/P* (modelliert)")
        ax.set_ylabel("Anzahl Peptide")
        ax.set_title("Verteilung L/P*")
        ax.legend()
        _format_publication_axes(ax)
        fig.tight_layout()
        fig.savefig(lp_path, dpi=300)
        plt.close(fig)
        generated.add(lp_path)

    placeholders = {
        plots_dir / "parity_log10.png": "Keine experimentellen HC50-Daten verfügbar.",
        plots_dir / "residuals_vs_logKp.png": "Residuen nicht berechenbar (fehlende HC50).",
        plots_dir / "error_bars_hc50.png": "Experimentelle HC50 oder Vorhersagen fehlen.",
        plots_dir / "lp_star_distribution.png": "Keine L/P*-Werte vorhanden.",
    }

    for path, message in placeholders.items():
        if path not in generated:
            fig, ax = plt.subplots(figsize=(6.2, 4.2))
            ax.axis("off")
            ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
            fig.tight_layout()
            fig.savefig(path, dpi=200)
            plt.close(fig)


def export_reports(
    agg_df: pd.DataFrame,
    cal_results: Optional[CalibrationResults],
    powerlaw: Optional[PowerLawResults],
    gamma_results: Optional[GammaModelResults],
    rank_results: Optional[RankModelResults],
    rank_results_enh: Optional[RankModelEnhancedResults],
    rank_results_knn: Optional[RankModelKNNResults],
    rank_results_combo: Optional[RankModelEnsembleResults],
    rank_results_v3: Optional[RankModelV3Results],
    lp_stats: Dict[str, Any],
    output_dir: Path,
    excel_name: str = "results_calibrated.xlsx",
    meta_info: Optional[Dict[str, Any]] = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    excel_path = output_dir / excel_name
    csv_path = output_dir / "results_calibrated.csv"

    required_cols = [
        "Peptide_ID",
        "n_reps",
        "Kp_eff_nm",
        "Kp_ads_nm",
        "HC50_pred_mean_uM",
        "HC50_pred_low_uM",
        "HC50_pred_high_uM",
        "HC50_cal_uM",
        "HC50_gamma_uM",
        "HC50_rank_uM",
        "HC50_rank_enh_uM",
        "HC50_rank_knn_uM",
        "HC50_rank_combo_uM",
        "HC50_rank_v3_uM",
        "HC50_exp_uM",
        "HC50_resid_log10",
        "HC50_gamma_resid_log10",
        "HC50_rank_resid_log10",
        "HC50_rank_enh_resid_log10",
        "HC50_rank_knn_resid_log10",
        "HC50_rank_combo_resid_log10",
        "HC50_rank_v3_resid_log10",
        "HC50_rank_score",
        "HC50_rank_enh_score",
        "HC50_rank_knn_score",
        "HC50_rank_combo_score",
        "HC50_rank_v3_score",
        "DeltaG_ads_kJ_per_mol",
        "z_ads_nm",
        "ESS_mean_median",
        "QC_pass_fraction",
        "Overlap_min_adjacent_median",
        "JS_mean_adjacent_median",
        "theta_at_1uM",
        "Gamma_star_obs",
        "Gamma_star_model",
        "LP_star_obs",
        "LP_star_model",
        "LP_star_Q25",
        "LP_star_Q50",
        "LP_star_Q75",
        "sequence_length",
        "HC50_meta_uM",
        "HC50_exp_source",
    ]
    for col in required_cols:
        if col not in agg_df.columns:
            agg_df[col] = np.nan

    df_main = agg_df[required_cols + ["replicate_ids", "replicate_tags"]].copy()
    df_main.to_csv(csv_path, index=False)

    df_two_param = None
    if 'HC50_cal_2p_uM' in agg_df.columns:
        cols_2p = required_cols.copy()
        if 'HC50_cal_2p_uM' not in cols_2p:
            cols_2p.insert(cols_2p.index('HC50_cal_uM'), 'HC50_cal_2p_uM')
        df_two_param = agg_df[[c for c in cols_2p if c in agg_df.columns]].copy()

    with pd.ExcelWriter(excel_path) as writer:
        df_main.to_excel(writer, sheet_name="aggregated", index=False)
        if df_two_param is not None:
            df_two_param.to_excel(writer, sheet_name="calibration_2p", index=False)
        if cal_results and cal_results.calibration_df is not None and not cal_results.calibration_df.empty:
            cal_results.calibration_df.to_excel(writer, sheet_name="calibration_rows", index=False)

    model_payload: Dict[str, Any] = {}
    if rank_results:
        model_payload['rank_model'] = {
            'feature_names': list(rank_results.feature_names),
            'weights': [float(w) for w in rank_results.weights.tolist()],
            'means': [float(m) for m in rank_results.feature_means.tolist()],
            'stds': [float(s) for s in rank_results.feature_stds.tolist()],
            'isotonic_x': [float(x) for x in rank_results.score_sorted.tolist()],
            'isotonic_y': [float(y) for y in rank_results.log_sorted.tolist()],
            'spearman_r': float(rank_results.spearman),
            'pearson_r': float(rank_results.pearson),
            'rmse_log10': float(rank_results.rmse_log10),
            'mae_log10': float(rank_results.mae_log10),
            'loo_spearman_r': float(rank_results.loo_spearman) if rank_results.loo_spearman is not None else None,
            'loo_rmse_log10': float(rank_results.loo_rmse_log10) if rank_results.loo_rmse_log10 is not None else None,
            'n_points': int(rank_results.n_points),
        }
    if rank_results_enh:
        model_payload['rank_model_enhanced'] = {
            'feature_names': list(rank_results_enh.feature_names),
            'coefficients': [float(c) for c in np.asarray(rank_results_enh.coefficients).ravel()],
            'intercept': float(rank_results_enh.intercept),
            'means': [float(m) for m in np.asarray(rank_results_enh.feature_means).ravel()],
            'stds': [float(s) for s in np.asarray(rank_results_enh.feature_stds).ravel()],
            'alpha': float(rank_results_enh.alpha),
            'spearman_r': float(rank_results_enh.spearman),
            'pearson_r': float(rank_results_enh.pearson),
            'rmse_log10': float(rank_results_enh.rmse_log10),
            'mae_log10': float(rank_results_enh.mae_log10),
            'loo_spearman_r': float(rank_results_enh.loo_spearman) if rank_results_enh.loo_spearman is not None else None,
            'loo_rmse_log10': float(rank_results_enh.loo_rmse_log10) if rank_results_enh.loo_rmse_log10 is not None else None,
            'n_points': int(rank_results_enh.n_points),
        }
    if rank_results_knn:
        model_payload['rank_model_knn'] = {
            'feature_names': list(rank_results_knn.feature_names),
            'k_neighbors': int(rank_results_knn.k_neighbors),
            'means': [float(m) for m in np.asarray(rank_results_knn.feature_means).ravel()],
            'stds': [float(s) for s in np.asarray(rank_results_knn.feature_stds).ravel()],
            'spearman_r': float(rank_results_knn.spearman),
            'pearson_r': float(rank_results_knn.pearson),
            'rmse_log10': float(rank_results_knn.rmse_log10),
            'mae_log10': float(rank_results_knn.mae_log10),
            'loo_spearman_r': float(rank_results_knn.loo_spearman) if rank_results_knn.loo_spearman is not None else None,
            'loo_rmse_log10': float(rank_results_knn.loo_rmse_log10) if rank_results_knn.loo_rmse_log10 is not None else None,
            'n_points': int(rank_results_knn.n_points),
        }
    if rank_results_combo:
        model_payload['rank_model_ensemble'] = {
            'components': list(rank_results_combo.component_names),
            'coefficients': [float(c) for c in np.asarray(rank_results_combo.coefficients).ravel()],
            'intercept': float(rank_results_combo.intercept),
            'spearman_r': float(rank_results_combo.spearman),
            'pearson_r': float(rank_results_combo.pearson),
            'rmse_log10': float(rank_results_combo.rmse_log10),
            'mae_log10': float(rank_results_combo.mae_log10),
            'loo_spearman_r': float(rank_results_combo.loo_spearman) if rank_results_combo.loo_spearman is not None else None,
            'loo_rmse_log10': float(rank_results_combo.loo_rmse_log10) if rank_results_combo.loo_rmse_log10 is not None else None,
            'n_points': int(rank_results_combo.n_points),
        }
    if rank_results_v3:
        model_payload['rank_model_v3'] = {
            'feature_names': list(rank_results_v3.feature_names),
            'coefficients': [float(c) for c in np.asarray(rank_results_v3.coefficients).ravel()],
            'intercept': float(rank_results_v3.intercept),
            'alpha': float(rank_results_v3.alpha),
            'l1_ratio': float(rank_results_v3.l1_ratio),
            'scaler_mean': [float(m) for m in np.asarray(rank_results_v3.scaler_mean).ravel()],
            'scaler_scale': [float(s) for s in np.asarray(rank_results_v3.scaler_scale).ravel()],
            'iso_x': [float(x) for x in np.asarray(rank_results_v3.iso_x).ravel()],
            'iso_y': [float(y) for y in np.asarray(rank_results_v3.iso_y).ravel()],
            'spearman_r': float(rank_results_v3.spearman),
            'pearson_r': float(rank_results_v3.pearson),
            'rmse_log10': float(rank_results_v3.rmse_log10),
            'mae_log10': float(rank_results_v3.mae_log10),
            'cv_spearman_r': float(rank_results_v3.cv_spearman) if rank_results_v3.cv_spearman is not None else None,
            'cv_rmse_log10': float(rank_results_v3.cv_rmse_log10) if rank_results_v3.cv_rmse_log10 is not None else None,
            'cv_mae_log10': float(rank_results_v3.cv_mae_log10) if rank_results_v3.cv_mae_log10 is not None else None,
            'sample_weight_mean': float(rank_results_v3.sample_weight_mean),
            'n_points': int(rank_results_v3.n_points),
        }
    if gamma_results:
        model_payload['gamma_model'] = {
            'alpha': float(gamma_results.alpha),
            'beta': float(gamma_results.beta),
            'alpha_ci95': [float(gamma_results.alpha_ci95[0]), float(gamma_results.alpha_ci95[1])],
            'beta_ci95': [float(gamma_results.beta_ci95[0]), float(gamma_results.beta_ci95[1])],
            'r_value': float(gamma_results.r_value),
            'rmse_log10': float(gamma_results.metrics.get('rmse_log10', np.nan)),
            'spearman_r': float(gamma_results.metrics.get('spearman_r', np.nan)),
            'n_points': int(gamma_results.n_points),
        }
    if model_payload:
        model_payload['generated'] = datetime.now(timezone.utc).isoformat()
        for path in [output_dir / "hc50_rank_model.yaml", CONFIG_DIR / "hc50_rank_model.yaml"]:
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                with open(path, 'w', encoding='utf-8') as f:
                    yaml.safe_dump(model_payload, f, sort_keys=False)
            except Exception:
                continue

    yaml_payload: Dict[str, Any] = {
        "adsorption": {
            "lp_star_range": [
                float(np.round(lp_stats.get("q25", np.nan), 2)),
                float(np.round(lp_stats.get("q75", np.nan), 2)),
            ],
            "lp_star_range_extended": [
                float(np.round(lp_stats.get("q10", np.nan), 2)),
                float(np.round(lp_stats.get("q90", np.nan), 2)),
            ] if np.isfinite(lp_stats.get("q10", np.nan)) and np.isfinite(lp_stats.get("q90", np.nan)) else None,
            "notes": f"Derived from n={lp_stats.get('n', 0)} peptides using Gamma-based model ({lp_stats.get('source', 'model')}).",
        },
        "gamma_model": None,
        "scale_factor_s": cal_results.scale_factor if cal_results else None,
    }
    if yaml_payload["adsorption"]["lp_star_range_extended"] is None:
        yaml_payload["adsorption"].pop("lp_star_range_extended")

    if gamma_results:
        gamma_metrics_nat = {
            key: (float(value) if isinstance(value, (np.floating, np.integer)) else value)
            for key, value in gamma_results.metrics.items()
        }
        yaml_payload["gamma_model"] = {
            "alpha": float(gamma_results.alpha),
            "beta": float(gamma_results.beta),
            "alpha_ci95": [float(gamma_results.alpha_ci95[0]), float(gamma_results.alpha_ci95[1])],
            "beta_ci95": [float(gamma_results.beta_ci95[0]), float(gamma_results.beta_ci95[1])],
            "r_value": float(gamma_results.r_value),
            **gamma_metrics_nat,
            "n_points": int(gamma_results.n_points),
        }
        yaml_payload.setdefault("mbar", {})["hc50_log10_rmse"] = float(gamma_metrics_nat.get("rmse_log10", np.nan))
        yaml_payload["mbar"]["spearman_r"] = float(gamma_metrics_nat.get("spearman_r", np.nan))
    elif cal_results:
        yaml_payload.setdefault("mbar", {})["hc50_log10_rmse"] = cal_results.metrics_after["rmse_log10"]
        yaml_payload["mbar"]["spearman_r"] = cal_results.metrics_after["spearman_r"]

    if rank_results:
        yaml_payload["rank_model"] = {
            "features": list(rank_results.feature_names),
            "weights": [float(w) for w in rank_results.weights],
            "spearman_r": float(rank_results.spearman),
            "pearson_r": float(rank_results.pearson),
            "rmse_log10": float(rank_results.rmse_log10),
            "loo_spearman_r": float(rank_results.loo_spearman) if rank_results.loo_spearman is not None else None,
            "n_points": rank_results.n_points,
        }
    if rank_results_enh:
        yaml_payload["rank_model_enhanced"] = {
            "features": list(rank_results_enh.feature_names),
            "coefficients": [float(c) for c in np.asarray(rank_results_enh.coefficients).ravel()],
            "intercept": float(rank_results_enh.intercept),
            "means": [float(m) for m in np.asarray(rank_results_enh.feature_means).ravel()],
            "stds": [float(s) for s in np.asarray(rank_results_enh.feature_stds).ravel()],
            "alpha": float(rank_results_enh.alpha),
            "spearman_r": float(rank_results_enh.spearman),
            "pearson_r": float(rank_results_enh.pearson),
            "rmse_log10": float(rank_results_enh.rmse_log10),
            "loo_spearman_r": float(rank_results_enh.loo_spearman) if rank_results_enh.loo_spearman is not None else None,
            "loo_rmse_log10": float(rank_results_enh.loo_rmse_log10) if rank_results_enh.loo_rmse_log10 is not None else None,
            "n_points": rank_results_enh.n_points,
        }
    if rank_results_knn:
        yaml_payload["rank_model_knn"] = {
            "features": list(rank_results_knn.feature_names),
            "k_neighbors": int(rank_results_knn.k_neighbors),
            "means": [float(m) for m in np.asarray(rank_results_knn.feature_means).ravel()],
            "stds": [float(s) for s in np.asarray(rank_results_knn.feature_stds).ravel()],
            "spearman_r": float(rank_results_knn.spearman),
            "pearson_r": float(rank_results_knn.pearson),
            "rmse_log10": float(rank_results_knn.rmse_log10),
            "loo_spearman_r": float(rank_results_knn.loo_spearman) if rank_results_knn.loo_spearman is not None else None,
            "loo_rmse_log10": float(rank_results_knn.loo_rmse_log10) if rank_results_knn.loo_rmse_log10 is not None else None,
            "n_points": rank_results_knn.n_points,
        }
    if rank_results_combo:
        yaml_payload["rank_model_ensemble"] = {
            "components": list(rank_results_combo.component_names),
            "coefficients": [float(c) for c in np.asarray(rank_results_combo.coefficients).ravel()],
            "intercept": float(rank_results_combo.intercept),
            "spearman_r": float(rank_results_combo.spearman),
            "pearson_r": float(rank_results_combo.pearson),
            "rmse_log10": float(rank_results_combo.rmse_log10),
            "loo_spearman_r": float(rank_results_combo.loo_spearman) if rank_results_combo.loo_spearman is not None else None,
            "loo_rmse_log10": float(rank_results_combo.loo_rmse_log10) if rank_results_combo.loo_rmse_log10 is not None else None,
            "n_points": rank_results_combo.n_points,
        }
    if rank_results_v3:
        yaml_payload["rank_model_v3"] = {
            "features": list(rank_results_v3.feature_names),
            "coefficients": [float(c) for c in np.asarray(rank_results_v3.coefficients).ravel()],
            "alpha": float(rank_results_v3.alpha),
            "l1_ratio": float(rank_results_v3.l1_ratio),
            "spearman_r": float(rank_results_v3.spearman),
            "pearson_r": float(rank_results_v3.pearson),
            "rmse_log10": float(rank_results_v3.rmse_log10),
            "cv_spearman_r": float(rank_results_v3.cv_spearman) if rank_results_v3.cv_spearman is not None else None,
            "cv_rmse_log10": float(rank_results_v3.cv_rmse_log10) if rank_results_v3.cv_rmse_log10 is not None else None,
            "n_points": rank_results_v3.n_points,
        }

    yaml_path = output_dir / "calibration_suggestions.yaml"
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(yaml_payload, f, sort_keys=False)

    summary_lines = [
        "HC50 calibration summary",
        "=======================",
        f"Peptides with experimental HC50: {int(agg_df['HC50_exp_uM'].notna().sum())}",
    ]
    if meta_info:
        summary_lines.append(
            f"Metadata matches: {meta_info.get('n_matched', 0)} / {meta_info.get('n_meta', 0)}"
        )
    if cal_results:
        summary_lines.extend(
            [
                f"Scale factor s (1p): {cal_results.scale_factor:.4f}",
                f"Shift (log10): {cal_results.shift_median_log10:.4f}",
                f"RMSE log10 (before -> after): {cal_results.metrics_before['rmse_log10']:.4f} -> {cal_results.metrics_after['rmse_log10']:.4f}",
                f"Spearman ρ (before -> after): {cal_results.metrics_before['spearman_r']:.4f} -> {cal_results.metrics_after['spearman_r']:.4f}",
            ]
        )
        if cal_results.weighted_scale_factor is not None:
            summary_lines.append(
                f"Weighted scale factor: {cal_results.weighted_scale_factor:.4f} (shift {cal_results.weighted_shift_log10:.4f})"
            )

    summary_lines.append("")
    if gamma_results:
        summary_lines.extend(
            [
                "Gamma-based model log Γ* = α·ΔG + β",
                f"n = {gamma_results.n_points}",
                f"α = {gamma_results.alpha:.4f} (95% CI {gamma_results.alpha_ci95[0]:.4f}, {gamma_results.alpha_ci95[1]:.4f})",
                f"β = {gamma_results.beta:.4f} (95% CI {gamma_results.beta_ci95[0]:.4f}, {gamma_results.beta_ci95[1]:.4f})",
                f"r = {gamma_results.r_value:.4f}",
                f"RMSE log10 = {gamma_results.metrics['rmse_log10']:.4f}",
                f"Spearman ρ = {gamma_results.metrics['spearman_r']:.4f}",
            ]
        )
    elif not cal_results:
        summary_lines.append("Keine Kalibration möglich (unzureichende Daten).")

    if rank_results:
        loo_rho_str = (
            f"{rank_results.loo_spearman:.4f}"
            if rank_results.loo_spearman is not None and np.isfinite(rank_results.loo_spearman)
            else "nan"
        )
        summary_lines.extend(
            [
                "",
                "Rank-basiertes Modell (Basis)",
                f"n = {rank_results.n_points}",
                f"Spearman ρ = {rank_results.spearman:.4f}",
                f"Pearson r = {rank_results.pearson:.4f}",
                f"RMSE log10 = {rank_results.rmse_log10:.4f}",
                f"LOO Spearman ρ = {loo_rho_str}",
                "Gewichte: " + ", ".join(
                    f"{name}={weight:.3f}" for name, weight in zip(rank_results.feature_names, rank_results.weights)
                ),
            ]
        )
        if rank_results.loo_rmse_log10 is not None and np.isfinite(rank_results.loo_rmse_log10):
            summary_lines.append(f"LOO RMSE log10 = {rank_results.loo_rmse_log10:.4f}")

    if rank_results_enh:
        loo_rho_enh = (
            f"{rank_results_enh.loo_spearman:.4f}"
            if rank_results_enh.loo_spearman is not None and np.isfinite(rank_results_enh.loo_spearman)
            else "nan"
        )
        summary_lines.extend(
            [
                "",
                "Rank-basiertes Modell (Enhanced)",
                f"n = {rank_results_enh.n_points}",
                f"Alpha = {rank_results_enh.alpha:.3f}",
                f"Spearman ρ = {rank_results_enh.spearman:.4f}",
                f"Pearson r = {rank_results_enh.pearson:.4f}",
                f"RMSE log10 = {rank_results_enh.rmse_log10:.4f}",
                f"LOO Spearman ρ = {loo_rho_enh}",
                f"Intercept = {rank_results_enh.intercept:.4f}",
                "Koeffizienten: " + ", ".join(
                    f"{name}={coef:.3f}" for name, coef in zip(rank_results_enh.feature_names, rank_results_enh.coefficients)
                ),
            ]
        )
        if rank_results_enh.loo_rmse_log10 is not None and np.isfinite(rank_results_enh.loo_rmse_log10):
            summary_lines.append(f"LOO RMSE log10 = {rank_results_enh.loo_rmse_log10:.4f}")

    if rank_results_knn:
        loo_rho_knn = (
            f"{rank_results_knn.loo_spearman:.4f}"
            if rank_results_knn.loo_spearman is not None and np.isfinite(rank_results_knn.loo_spearman)
            else "nan"
        )
        summary_lines.extend(
            [
                "",
                "Rank-basiertes Modell (kNN)",
                f"n = {rank_results_knn.n_points}",
                f"k = {rank_results_knn.k_neighbors}",
                f"Spearman ρ = {rank_results_knn.spearman:.4f}",
                f"Pearson r = {rank_results_knn.pearson:.4f}",
                f"RMSE log10 = {rank_results_knn.rmse_log10:.4f}",
                f"LOO Spearman ρ = {loo_rho_knn}",
            ]
        )
        if rank_results_knn.loo_rmse_log10 is not None and np.isfinite(rank_results_knn.loo_rmse_log10):
            summary_lines.append(f"LOO RMSE log10 = {rank_results_knn.loo_rmse_log10:.4f}")

    if rank_results_combo:
        loo_rho_combo = (
            f"{rank_results_combo.loo_spearman:.4f}"
            if rank_results_combo.loo_spearman is not None and np.isfinite(rank_results_combo.loo_spearman)
            else "nan"
        )
        summary_lines.extend(
            [
                "",
                "Rank-basiertes Modell (Ensemble)",
                f"n = {rank_results_combo.n_points}",
                f"Spearman ρ = {rank_results_combo.spearman:.4f}",
                f"Pearson r = {rank_results_combo.pearson:.4f}",
                f"RMSE log10 = {rank_results_combo.rmse_log10:.4f}",
                f"LOO Spearman ρ = {loo_rho_combo}",
                f"Intercept = {rank_results_combo.intercept:.4f}",
                "Koeffizienten: " + ", ".join(
                    f"{name}={coef:.3f}" for name, coef in zip(rank_results_combo.component_names, rank_results_combo.coefficients)
                ),
            ]
        )
        if rank_results_combo.loo_rmse_log10 is not None and np.isfinite(rank_results_combo.loo_rmse_log10):
            summary_lines.append(f"LOO RMSE log10 = {rank_results_combo.loo_rmse_log10:.4f}")

    if rank_results_v3:
        summary_lines.extend(
            [
                "",
                "Rank-basiertes Modell (v3)",
                f"n = {rank_results_v3.n_points}",
                f"α = {rank_results_v3.alpha:.4f}",
                f"l1_ratio = {rank_results_v3.l1_ratio:.3f}",
                f"Spearman ρ = {rank_results_v3.spearman:.4f}",
                f"Pearson r = {rank_results_v3.pearson:.4f}",
                f"RMSE log10 = {rank_results_v3.rmse_log10:.4f}",
                f"CV Spearman ρ = {rank_results_v3.cv_spearman:.4f}" if rank_results_v3.cv_spearman is not None else "CV Spearman ρ = n/a",
                f"CV RMSE log10 = {rank_results_v3.cv_rmse_log10:.4f}" if rank_results_v3.cv_rmse_log10 is not None else "CV RMSE log10 = n/a",
                f"CV MAE log10 = {rank_results_v3.cv_mae_log10:.4f}" if rank_results_v3.cv_mae_log10 is not None else "CV MAE log10 = n/a",
            ]
        )

    summary_lines.extend(
        [
            "",
            "Derived L/P* statistics:",
            f"n = {lp_stats.get('n', 0)}",
            f"Median = {lp_stats.get('median', float('nan')):.3f}",
            f"IQR = {lp_stats.get('iqr', float('nan')):.3f}",
            f"[Q25, Q75] = [{lp_stats.get('q25', float('nan')):.3f}, {lp_stats.get('q75', float('nan')):.3f}]",
            f"[Q10, Q90] = [{lp_stats.get('q10', float('nan')):.3f}, {lp_stats.get('q90', float('nan')):.3f}]",
        ]
    )
    summary_path = output_dir / "calibration_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")


def leave_one_out_cv(
    cal_results: Optional[CalibrationResults],
    output_dir: Path,
    hc50_cap: float,
) -> None:
    if not cal_results or cal_results.calibration_df is None:
        return
    calib = cal_results.calibration_df
    if calib.shape[0] < 3:
        return
    cap = float(hc50_cap)
    exp_clip, pred_clip = _censor_pair(
        calib["HC50_exp_uM"],
        calib["HC50_pred_mean_uM"],
        cap,
    )
    log_exp = np.log10(np.clip(exp_clip, 1e-12, None))
    log_pred = np.log10(np.clip(pred_clip, 1e-12, None))
    delta = log_exp - log_pred

    cv_resid = []
    for i in range(delta.size):
        others = np.delete(delta, i)
        shift = np.median(others)
        scale = 10 ** shift
        pred_uM = calib.iloc[i]["HC50_pred_mean_uM"] * scale
        raw_exp = calib.iloc[i].get("HC50_exp_uM_raw", calib.iloc[i].get("HC50_exp_uM"))
        if pd.notna(raw_exp) and raw_exp > cap:
            pred_uM = cap
        pred_uM = min(pred_uM, cap)
        resid_log = math.log10(max(pred_uM, 1e-12)) - log_exp[i]
        cv_resid.append(resid_log)

    if not cv_resid:
        return

    cv_rmse = float(np.sqrt(np.mean(np.square(cv_resid))))
    cv_mae = float(np.mean(np.abs(cv_resid)))
    cv_path = output_dir / "calibration_cv.txt"
    with open(cv_path, "w", encoding="utf-8") as f:
        f.write("Leave-one-peptide-out cross-validation (1p scaling)\n")
        f.write(f"Folds: {len(cv_resid)}\n")
        f.write(f"RMSE log10: {cv_rmse:.4f}\n")
        f.write(f"MAE log10: {cv_mae:.4f}\n")


def detect_default_input(path_hint: Optional[str]) -> Path:
    if path_hint:
        return Path(path_hint)
    candidates = [
        Path("results_summary.csv"),
        Path("results_summary_1.csv"),
        Path("results_summary.xlsx"),
        Path("results_summary1.xlsx"),
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("No results summary file found (csv or xlsx).")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate HC50 predictions using aggregated PMF metrics."
    )
    parser.add_argument("--input", help="Path to results summary (CSV preferred).")
    parser.add_argument(
        "--meta",
        default="data/raw/peptides.csv",
        help="Path to DBAASP metadata CSV (default: data/raw/peptides.csv)",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis",
        help="Directory for reports (default: analysis)",
    )
    parser.add_argument(
        "--plots-dir",
        default=None,
        help="Optional override for plot directory (default: <output-dir>/plots)",
    )
    parser.add_argument(
        "--area-per-lipid",
        type=float,
        default=0.62,
        help="Area per lipid in nm^2 for L/P* back-computation (default 0.62).",
    )
    parser.add_argument(
        "--skip-two-parameter",
        action="store_true",
        help="Skip the optional two-parameter power-law calibration.",
    )
    parser.add_argument(
        "--calibration-max-um",
        type=float,
        default=150.0,
        help="Upper HC50 limit (µM) to include in calibration fits (default 150).",
    )
    args = parser.parse_args()

    input_path = detect_default_input(args.input)
    output_dir = Path(args.output_dir)
    plots_dir = Path(args.plots_dir) if args.plots_dir else (output_dir / "plots")

    raw_df = load_input_table(input_path)
    agg_df = aggregate_replicates(raw_df)

    meta_df = load_metadata(Path(args.meta)) if args.meta else None
    agg_df, meta_info = merge_with_metadata(agg_df, meta_df)
    hc50_cap = _resolve_hc50_cap(args.calibration_max_um)
    if 'HC50_exp_uM' in agg_df.columns:
        agg_df['HC50_exp_uM_raw'] = agg_df['HC50_exp_uM']
        agg_df['HC50_exp_uM'] = agg_df['HC50_exp_uM'].clip(upper=hc50_cap)
    if 'Peptide_key' in agg_df.columns:
        agg_df.drop(columns=['Peptide_key'], inplace=True)

    agg_df, cal_results = fit_scale_1p(
        agg_df,
        calibration_max_um=args.calibration_max_um,
        hc50_cap=hc50_cap,
    )
    powerlaw_results = None
    if not args.skip_two_parameter:
        agg_df, powerlaw_results = fit_powerlaw_2p(
            agg_df,
            calibration_max_um=args.calibration_max_um,
            hc50_cap=hc50_cap,
        )

    agg_df, gamma_results = fit_gamma_model(
        agg_df,
        area_per_lipid_nm2=args.area_per_lipid,
        calibration_max_um=args.calibration_max_um,
        hc50_cap=hc50_cap,
    )
    agg_df, rank_results = fit_rank_model(agg_df, hc50_cap=hc50_cap)
    agg_df, rank_results_enh = fit_rank_model_enhanced(agg_df, hc50_cap=hc50_cap)
    agg_df, rank_results_knn = fit_rank_model_knn(agg_df, hc50_cap=hc50_cap)
    agg_df, rank_results_combo = fit_rank_model_ensemble(agg_df, hc50_cap=hc50_cap)
    agg_df, rank_results_v3 = fit_rank_model_v3(agg_df, hc50_cap=hc50_cap)
    agg_df, lp_stats = compute_lp_stats(agg_df)

    make_plots(
        agg_df,
        cal_results,
        gamma_results,
        powerlaw_results,
        rank_results,
        rank_results_enh,
        rank_results_knn,
        rank_results_combo,
        rank_results_v3,
        plots_dir,
        hc50_cap=hc50_cap,
    )
    export_reports(
        agg_df,
        cal_results,
        powerlaw_results,
        gamma_results,
        rank_results,
        rank_results_enh,
        rank_results_knn,
        rank_results_combo,
        rank_results_v3,
        lp_stats,
        output_dir,
        meta_info=meta_info,
    )
    leave_one_out_cv(cal_results, output_dir, hc50_cap=hc50_cap)


if __name__ == "__main__":
    main()
