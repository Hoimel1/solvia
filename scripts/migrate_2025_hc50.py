#!/usr/bin/env python3
from __future__ import annotations

"""
Migration script: enrich existing pmf_analysis_results.yaml with adsorption.kp_nm,
adsorption.kp_eff_nm, adsorption.hc50_uM_range (if missing) using analysis.hc50.

Usage:
  # Run from repo root; scans recursively under base
  python scripts/migrate_2025_hc50.py /path/to/base --write [--verbose]
"""

import argparse
import sys
from pathlib import Path
from typing import Optional
import yaml
import numpy as np

# Ensure repository root is on sys.path so that "analysis.hc50" can be imported
REPO_ROOT = Path(__file__).parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from analysis.hc50 import AdsorptionParams, compute_kp_ads_nm, compute_hc50_uM_from_kp
except Exception:
    # Retry once after adjusting sys.path explicitly
    try:
        if str(REPO_ROOT) not in sys.path:
            sys.path.insert(0, str(REPO_ROOT))
        from analysis.hc50 import AdsorptionParams, compute_kp_ads_nm, compute_hc50_uM_from_kp
    except Exception:
        AdsorptionParams = None


def load_config(base: Path) -> dict:
    cfg_path = base / "config" / "pmf_standard_config.yaml"
    if not cfg_path.exists():
        # try repo root relative to this script
        cfg_path = Path(__file__).parent.parent / "config" / "pmf_standard_config.yaml"
    if cfg_path.exists():
        with open(cfg_path, 'r') as f:
            return yaml.safe_load(f) or {}
    return {}


def migrate_file(path: Path, cfg: dict, write: bool) -> tuple[bool, str]:
    with open(path, 'r') as f:
        data = yaml.safe_load(f) or {}

    ads = data.get('adsorption') or {}
    pmf = (data.get('pmf_data') or {})
    z = np.asarray(pmf.get('z') or [])
    W = np.asarray(pmf.get('pmf') or [])
    counts = np.asarray(pmf.get('bin_counts') or []) if pmf.get('bin_counts') is not None else None

    # skip if already present
    already = ('kp_nm' in ads and 'kp_eff_nm' in ads and 'hc50_uM_range' in ads)
    if already:
        return False, 'already_up_to_date'
    if z.size == 0 or W.size == 0 or z.shape != W.shape:
        return False, 'missing_or_mismatched_pmf_arrays'

    # Build params from config
    ana = ((cfg.get('pmf') or {}).get('analysis') or {})
    ads_cfg = (ana.get('adsorption') or {})
    feat_cfg = (ana.get('feature_params') or {})
    sim_cfg = (ana.get('simulation') or {})
    ap = AdsorptionParams(
        temperature_K=float(sim_cfg.get('temperature', 310.0)),
        energy_band_kj=float(ads_cfg.get('energy_band_kj', 3.0)),
        area_per_lipid_nm2=float(ads_cfg.get('area_per_lipid_nm2', 0.62)),
        bilayer_factor=float(ads_cfg.get('bilayer_factor', 2.0)),
        lp_star_range=tuple((ads_cfg.get('lp_star_range') or [154.0, 515.0])[:2]),
        z_ads_max_nm=float(feat_cfg.get('ads_max_z', 2.0)),
        z_lo_hi_nm=(float((ads_cfg.get('membrane') or {}).get('z_lo', -1.5)),
                    float((ads_cfg.get('membrane') or {}).get('z_hi', 1.5))),
        min_bin_count=int((ana.get('plot') or {}).get('pmf_min_bin_count', 25)),
    )

    kp_nm, meta = compute_kp_ads_nm(z, W, counts, ap)
    hc = compute_hc50_uM_from_kp(kp_nm, ap)
    # update adsorption
    ads.setdefault('params', {
        'area_per_lipid_nm2': ap.area_per_lipid_nm2,
        'bilayer_factor': ap.bilayer_factor,
        'lp_star_range': [float(ap.lp_star_range[0]), float(ap.lp_star_range[1])],
    })
    ads['kp_nm'] = float(kp_nm)
    if hc.get('ok'):
        ads['kp_eff_nm'] = float(hc['kp_eff_nm'])
        ads['hc50_uM_range'] = [float(hc['hc50_uM_range'][0]), float(hc['hc50_uM_range'][1])]
    data['adsorption'] = ads

    if write:
        with open(path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False)
    return True, 'updated'


def main():
    p = argparse.ArgumentParser(description="Migrate pmf_analysis_results.yaml to include HC50/Kp adsorption fields.")
    p.add_argument('base', help='Base directory to scan recursively')
    p.add_argument('--write', action='store_true', help='Write changes back to files')
    p.add_argument('--verbose', action='store_true', help='Print per-file status (up-to-date/skipped)')
    args = p.parse_args()

    base = Path(args.base)
    if not base.exists():
        print(f"Base path not found: {base}")
        sys.exit(1)
    if AdsorptionParams is None:
        print("analysis.hc50 not importable. Ensure repository is intact.")
        sys.exit(1)

    cfg = load_config(Path.cwd())
    total = 0
    updated = 0
    skipped = 0
    for yml in base.rglob('pmf_analysis_results.yaml'):
        try:
            total += 1
            changed, reason = migrate_file(yml, cfg, write=args.write)
            if changed:
                updated += 1
                print(("Updated" if args.write else "Needs update") + f": {yml}")
            else:
                skipped += 1
                if args.verbose:
                    msg = "Up-to-date" if reason == 'already_up_to_date' else f"Skipped ({reason})"
                    print(f"{msg}: {yml}")
        except Exception as e:
            print(f"Skip {yml}: {e}")
    print(f"Done. scanned={total}, updated={updated}, skipped={skipped}.")


if __name__ == '__main__':
    main()
