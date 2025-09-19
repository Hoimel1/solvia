from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any

# Physical constants
K_BOLTZ_KJ_PER_MOL_K = 0.008314462618  # kJ/mol/K
NA_PER_L_NM3 = 0.602214                # 1 M -> nm^-3


@dataclass
class AdsorptionParams:
    temperature_K: float = 310.0
    energy_band_kj: float = 3.0
    area_per_lipid_nm2: float = 0.62
    bilayer_factor: float = 2.0
    lp_star_range: Tuple[float, float] = (154.0, 515.0)
    z_ads_max_nm: float = 2.0
    z_lo_hi_nm: Tuple[float, float] = (-1.5, 1.5)
    min_bin_count: int = 25
    exp_clip_kj: float = 30.0  # numerical stabilization


def _beta(T: float) -> float:
    return 1.0 / (K_BOLTZ_KJ_PER_MOL_K * T)


def find_adsorption_minimum(z_nm: np.ndarray,
                            W_kj: np.ndarray,
                            counts: Optional[np.ndarray],
                            z_max_nm: float) -> Optional[int]:
    """Index of PMF minimum on water side in [0, z_max_nm]."""
    mask = (z_nm >= 0.0) & (z_nm <= z_max_nm)
    if counts is not None:
        mask &= (np.asarray(counts) >= 1)
    if not np.any(mask):
        return None
    idx_rel = int(np.argmin(W_kj[mask]))
    return int(np.nonzero(mask)[0][idx_rel])


def _window_mask(z_nm: np.ndarray,
                 W_kj: np.ndarray,
                 min_idx: int,
                 params: AdsorptionParams,
                 counts: Optional[np.ndarray]) -> np.ndarray:
    z_lo, z_hi = params.z_lo_hi_nm
    band = params.energy_band_kj
    z_min = float(z_nm[min_idx])
    W_min = float(W_kj[min_idx])

    # water-side only, clamp to membrane bounds and adsorption search range
    mask = (z_nm >= 0.0) & (z_nm <= params.z_ads_max_nm) & (z_nm >= z_lo) & (z_nm <= z_hi)
    # energy band around adsorption minimum
    mask &= (W_kj <= (W_min + band))
    # minimum support if available
    if counts is not None:
        mask &= (np.asarray(counts) >= params.min_bin_count)
    return mask


def compute_kp_ads_nm(z_nm: np.ndarray,
                      W_kj: np.ndarray,
                      counts: Optional[np.ndarray],
                      params: AdsorptionParams) -> Tuple[float, Dict[str, Any]]:
    """Compute Kp_ads in units of nm (length integral) on the water side.

    Integrand: exp(-beta*W) - 1 over the adsorption energy band around the
    water-side minimum. Returns (kp_nm, meta).
    """
    z_nm = np.asarray(z_nm, dtype=float)
    W_kj = np.asarray(W_kj, dtype=float)
    counts = None if counts is None else np.asarray(counts)

    beta = _beta(params.temperature_K)
    min_idx = find_adsorption_minimum(z_nm, W_kj, counts, params.z_ads_max_nm)
    if min_idx is None:
        return 0.0, {"ok": False, "reason": "no_adsorption_minimum"}

    mask = _window_mask(z_nm, W_kj, min_idx, params, counts)
    if not np.any(mask):
        return 0.0, {"ok": False, "reason": "empty_band", "z_min_nm": float(z_nm[min_idx]), "W_min_kj_per_mol": float(W_kj[min_idx])}

    z_win = z_nm[mask]
    W_win = W_kj[mask].copy()

    # Stabilize exponentials: clip to avoid extremely negative W in flat regions
    W_win = np.maximum(W_win, -params.exp_clip_kj)
    integrand = np.exp(-beta * W_win) - 1.0
    kp_nm = float(np.trapz(integrand, z_win))

    meta = {
        "ok": True,
        "z_min_nm": float(z_nm[min_idx]),
        "W_min_kj_per_mol": float(W_kj[min_idx]),
        "n_bins": int(np.sum(mask)),
    }
    return max(kp_nm, 0.0), meta


def compute_hc50_uM_from_kp(kp_nm: float, params: AdsorptionParams) -> Dict[str, Any]:
    """Compute HC50 band (µM) from Kp (nm) and project constants.

    Returns dict with ok flag, kp_eff_nm, and hc50_uM_range=(low, high).
    """
    bilayer = float(params.bilayer_factor)
    kp_eff = kp_nm * bilayer
    if not np.isfinite(kp_eff) or kp_eff <= 0.0:
        return {"ok": False, "reason": "kp_eff_nonpositive", "kp_eff_nm": float(kp_eff)}

    Lmin, Lmax = float(params.lp_star_range[0]), float(params.lp_star_range[1])
    # Gamma* = 1 / (A_lipid * L/P*)
    gamma_star_min = 1.0 / (params.area_per_lipid_nm2 * Lmin)  # larger -> larger HC50
    gamma_star_max = 1.0 / (params.area_per_lipid_nm2 * Lmax)

    # c = Gamma*/(NA * Kp_eff) [M]; Convert to µM
    c_high_uM = 1e6 * gamma_star_min / (NA_PER_L_NM3 * kp_eff)
    c_low_uM = 1e6 * gamma_star_max / (NA_PER_L_NM3 * kp_eff)

    return {
        "ok": True,
        "kp_nm": float(kp_nm),
        "kp_eff_nm": float(kp_eff),
        "hc50_uM_range": (float(c_low_uM), float(c_high_uM)),
        "params": {
            "area_per_lipid_nm2": float(params.area_per_lipid_nm2),
            "bilayer_factor": float(params.bilayer_factor),
            "lp_star_range": (float(Lmin), float(Lmax)),
        },
    }

