

#!/usr/bin/env python3
"""
Create PMF meta file (pmf_meta.yaml) from an existing PMF summary.json.

Usage:
  python3 scripts/analysis/create_pmf_meta.py RUN_DIR [--tag TAG] [--meta-out PATH]

RUN_DIR must contain pmf/summary.json (written by 08_run_pmf.py).
"""

import os, sys, json, yaml, argparse

def main():
    p = argparse.ArgumentParser(description="Create pmf_meta.yaml from pmf/summary.json")
    p.add_argument("run_dir", help="Simulation run directory (with pmf/summary.json)")
    p.add_argument("--tag", default=None, help="Replicate tag (default from summary.json)")
    p.add_argument("--meta-out", default=None, help="Output YAML path (default run_dir/pmf_meta.yaml)")
    args = p.parse_args()

    run_dir = args.run_dir
    summ_path = os.path.join(run_dir, "pmf", "summary.json")
    if not os.path.exists(summ_path):
        sys.exit(f"Error: {summ_path} not found.")

    with open(summ_path) as f:
        summ = json.load(f)

    tag = args.tag or summ.get("tag", "rep1")
    T_K = 310.0  # default, adjust if needed
    k = float(summ.get("spring_k", 700.0))

    rep = {"name": str(tag), "windows": []}
    for w in summ.get("windows", []):
        if not w.get("ok", False):
            continue
        wdir = os.path.join(run_dir, w["dir"])
        xvg = os.path.join(wdir, "umbrella_pullx.xvg")
        rep["windows"].append({
            "center_nm": float(w.get("center_nm", 0.0)),
            "k_kj_mol_nm2": k,
            "xvg": xvg,
            "skip_ns": 1.0,
            "stride": 1,
        })

    rep["windows"].sort(key=lambda a: a["center_nm"], reverse=True)

    meta = {
        "temperature_K": T_K,
        "bin_width_nm": 0.06,
        "pad_nm": 0.05,
        "drop_empty_bins": True,
        "decorrelate_for_fit": False,
        "bootstrap": 300,
        "z_bulk_min_nm": 2.4,
        "ranges": {
            "surf_min": [-0.4, 0.8],
            "head_min": [-1.4, -0.2],
        },
        "replicates": [rep],
    }

    out = args.meta_out or os.path.join(run_dir, "pmf_meta.yaml")
    with open(out, "w") as f:
        yaml.safe_dump(meta, f, sort_keys=False)
    print(f"Wrote {out} with {len(rep['windows'])} windows.")

if __name__ == "__main__":
    main()