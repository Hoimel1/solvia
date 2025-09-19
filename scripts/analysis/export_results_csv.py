#!/usr/bin/env python3
from __future__ import annotations

import sys
import csv
import argparse
from pathlib import Path
from typing import List, Dict, Any
import yaml


def _safe_get(d: dict, path: List[str], default=None):
    cur: Any = d
    try:
        for k in path:
            if cur is None:
                return default
            cur = cur.get(k)
        return cur if cur is not None else default
    except Exception:
        return default


def _infer_peptide_and_replicate(yaml_path: Path) -> tuple[str, str]:
    # Heuristic: peptide = run dir name prefix (e.g., solvia_XX_run_1 → solvia_XX)
    # replicate: from parent dir name if matches replicate_*, else "aggregate" or "unknown"
    try:
        parts = list(yaml_path.parts)
        # replicate is parent dir of results when under pmf/replicate_X
        rep = "unknown"
        if len(parts) >= 2 and parts[-2].startswith('replicate_'):
            rep = parts[-2].split('_', 1)[-1]
        elif yaml_path.name == 'pmf_analysis_aggregate.yaml':
            rep = 'aggregate'
        # peptide: search upward for simulations/<run> or top-level run dir
        pep = yaml_path.parent
        while pep.parent != pep and pep.name not in ('simulations', ''):
            if pep.parent.name == 'simulations':
                run_name = pep.name
                break
            pep = pep.parent
        else:
            run_name = yaml_path.parent.name
        # Extract peptide id prefix before _run
        pfx = run_name.split('_run')[0]
        return pfx, rep
    except Exception:
        return "unknown", "unknown"


def collect_rows(base: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    # Include both replicate and aggregate results
    patterns = ['**/pmf_analysis_results.yaml', '**/pmf_analysis_aggregate.yaml']
    files: List[Path] = []
    for pat in patterns:
        files.extend(base.glob(pat))
    for yml in sorted(files):
        try:
            with open(yml, 'r') as f:
                data = yaml.safe_load(f) or {}
        except Exception:
            continue
        peptide, replicate = _infer_peptide_and_replicate(yml)
        features = data.get('features', {}) or {}
        ads = data.get('adsorption', {}) or {}
        qc = data.get('quality_metrics', {}) or {}
        label = _safe_get(data, ['classification', 'label']) or _safe_get(data, ['classification', 'severity'])
        qc_pass = _safe_get(qc, ['qc_status', 'passed'])
        hc = ads.get('hc50_uM_range')
        lo_uM, hi_uM = (hc[0], hc[1]) if (isinstance(hc, (list, tuple)) and len(hc) == 2) else (None, None)
        row = {
            'peptide': peptide,
            'replicate': replicate,
            'kp_nm': ads.get('kp_nm'),
            'kp_eff_nm': ads.get('kp_eff_nm') or ads.get('kp_ads_nm'),
            'hc50_low_uM': lo_uM,
            'hc50_high_uM': hi_uM,
            'delta_g_ads': features.get('delta_g_ads'),
            'z_ads': features.get('z_ads'),
            'label': label,
            'qc_pass': qc_pass,
            'path': str(yml)
        }
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description='Export PMF analysis results to CSV (HC50/Kp).')
    ap.add_argument('base', help='Base directory to scan recursively')
    ap.add_argument('-o', '--out', help='Output CSV file (default: stdout)')
    args = ap.parse_args()

    base = Path(args.base)
    if not base.exists():
        print(f"Base path does not exist: {base}", file=sys.stderr)
        sys.exit(1)
    rows = collect_rows(base)
    cols = ['peptide', 'replicate', 'kp_nm', 'kp_eff_nm', 'hc50_low_uM', 'hc50_high_uM', 'delta_g_ads', 'z_ads', 'label', 'qc_pass', 'path']
    if args.out:
        with open(args.out, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()

