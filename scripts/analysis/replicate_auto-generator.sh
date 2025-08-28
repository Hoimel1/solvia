RUN_DIR="simulations/solvia_12_run_1"
python3 - <<'PY'
import os, json, yaml, glob, numpy as np, math

run_dir = os.environ.get("RUN_DIR", "simulations/solvia_12_run_1")
files = sorted(glob.glob(f"{run_dir}/pmf/summary*.json")) or [f"{run_dir}/pmf/summary.json"]

def read_z(xvg, skip_ns=1.0):
    t, z = [], []
    with open(xvg) as f:
        for line in f:
            if not line or line[0] in ('@', '#'): continue
            parts = line.split()
            if len(parts) < 2: continue
            t.append(float(parts[0]))  # ps
            z.append(float(parts[1]))  # nm
    t = np.asarray(t); z = np.asarray(z)
    if t.size == 0: return z
    mask = t >= (skip_ns*1000.0)  # ps
    return z[mask]

# sammle alle z
Z = []
reps = []
for sf in files:
    s = json.load(open(sf))
    tag = s.get("tag") or os.path.splitext(os.path.basename(sf))[0]
    k  = float(s.get("spring_k", 700.0))
    rep = {"name": str(tag), "windows": [], "k": k}
    for w in s.get("windows", []):
        if not w.get("ok", True): continue
        wdir = os.path.join(run_dir, w["dir"])
        xvg  = os.path.join(wdir, "umbrella_pullx.xvg")
        rep["windows"].append({"center_nm": float(w["center_nm"]), "xvg": xvg, "skip_ns": 1.0, "stride": 1})
        if os.path.exists(xvg):
            Z.append(read_z(xvg, 1.0))
    reps.append(rep)

# --- Overlap precheck (a priori) ---
KB_KJ_MOL_K = 0.008314462618  # kJ/mol/K
T_K = float(os.environ.get("PMF_T_K", 310.0))
overlap_target = float(os.environ.get("PMF_OVERLAP_TARGET", 0.10))

def Phi(x):  # standard normal CDF
    return 0.5*(1.0 + math.erf(x/math.sqrt(2.0)))

def expected_overlap(dz, sigma):
    # Two-sided overlap for neighboring harmonic umbrellas with equal k
    return 2.0*(1.0 - Phi(dz/(2.0*sigma)))

overlap_pairs_below = []
min_overlap = 1.0
suggest_midpoints = []
for rep in reps:
    k = float(rep["k"])
    sigma = math.sqrt(KB_KJ_MOL_K * T_K / k)
    centers = sorted([w["center_nm"] for w in rep["windows"]], reverse=True)
    for i in range(len(centers)-1):
        z_hi, z_lo = centers[i], centers[i+1]
        dz = abs(z_hi - z_lo)
        O = expected_overlap(dz, sigma)
        min_overlap = min(min_overlap, O)
        if O < overlap_target:
            midpoint = round((z_hi + z_lo) / 2.0, 3)
            overlap_pairs_below.append({
                "rep": rep["name"],
                "z_hi_nm": z_hi,
                "z_lo_nm": z_lo,
                "dz_nm": round(dz, 3),
                "O": round(O, 3),
                "suggest_midpoint_nm": midpoint,
            })
            suggest_midpoints.append(midpoint)

# If there are theoretical low-overlap pairs, do not hide with overly wide bins
cap_binwidth_if_low_overlap = len(overlap_pairs_below) > 0

Z = np.concatenate([z for z in Z if z.size], axis=0) if Z else np.array([])
# Default-Fallbacks
bin_width = 0.06
pad = 0.05

if Z.size:
    zmin, zmax = float(Z.min()), float(Z.max())
    # Start mit 0.06 nm; nur wenn keine theoretischen Overlap-Lücken vorliegen,
    # darf bis 0.10 verbreitert werden. Sonst bei 0.08 kappen, um echte Gaps nicht zu maskieren.
    width_candidates = (0.06, 0.07, 0.08, 0.10)
    if cap_binwidth_if_low_overlap:
        width_candidates = (0.06, 0.07, 0.08)
    for bw in width_candidates:
        edges = np.arange(zmin - pad, zmax + pad + 1e-9, bw)
        counts, _ = np.histogram(Z, bins=edges)
        # „interne“ Bins = ohne allererste/allerletzte
        internal = counts[1:-1] if counts.size > 2 else counts
        if internal.size == 0 or (internal > 0).all():
            bin_width = bw
            break

meta = {
  "temperature_K": 310,
  "bin_width_nm": round(bin_width, 3),
  "pad_nm": round(pad, 3),
  "drop_empty_bins": True,          # Sicherheit: äußere leere Bins entfernen
  "decorrelate_for_fit": False,     # nicht vorab ausdünnen (vergrößert Leer-Bin-Risiko)
  "bootstrap": 1000,
  "z_bulk_min_nm": 2.4,
  "ranges": {"surf_min": [-0.4, 0.8], "head_min": [-1.4, -0.2]},
  "overlap_precheck": {
    "target": overlap_target,
    "min_expected_overlap": round(min_overlap, 3),
    "pairs_below_target": overlap_pairs_below,
    "suggest_intermediate_centers_nm": sorted(set(suggest_midpoints)),
  },
  "replicates": []
}
for rep in reps:
    repd = {"name": rep["name"], "windows": []}
    for w in sorted(rep["windows"], key=lambda x: x["center_nm"], reverse=True):
        repd["windows"].append({
          "center_nm": float(w["center_nm"]),
          "k_kj_mol_nm2": float(rep["k"]),
          "xvg": w["xvg"],
          "skip_ns": w["skip_ns"],
          "stride": w["stride"]
        })
    meta["replicates"].append(repd)

out = f"{run_dir}/pmf_meta.yaml"
with open(out, "w") as f: yaml.safe_dump(meta, f, sort_keys=False)
print(f"{out} geschrieben (bin_width_nm={meta['bin_width_nm']}, pad_nm={meta['pad_nm']})")
if overlap_pairs_below:
    print(f"Overlap-Precheck: min O={min_overlap:.3f} < target {overlap_target:.2f}; vorgeschlagene Zwischenfenster: {sorted(set(suggest_midpoints))}")
else:
    print(f"Overlap-Precheck: min O={min_overlap:.3f} ≥ target {overlap_target:.2f}; keine Zwischenfenster nötig.")


# --- Auto-Insert: write per-replicate lists of missing midpoint windows ---
# Only create additional windows that are truly missing (no xvg yet).

def _norm(x: float) -> float:
    v = float(x)
    return 0.0 if abs(v) < 5e-4 else round(v, 3)

def _wdir_for(z: float) -> str:
    return os.path.join(run_dir, "pmf", "windows", f"z_{_norm(z):+0.2f}")

for rep in reps:
    # Existing centers of this replicate
    existing = sorted({_norm(w["center_nm"]) for w in rep["windows"]}, reverse=True)
    # Suggested midpoints *for this replicate*
    mids = sorted({_norm(p["suggest_midpoint_nm"]) for p in overlap_pairs_below if p["rep"] == rep["name"]}, reverse=True)

    # Filter: only those midpoint windows that do not already exist on disk
    missing = []
    for z in mids:
        wdir = _wdir_for(z)
        xvg = os.path.join(wdir, "umbrella_pullx.xvg")
        if not os.path.exists(xvg):
            missing.append(z)

    if not missing:
        continue

    # Write YAML override and TXT helper
    ypath = os.path.join(run_dir, "pmf", f"auto_windows_{rep['name']}.yaml")
    tpath = os.path.join(run_dir, "pmf", f"auto_windows_{rep['name']}.txt")
    os.makedirs(os.path.dirname(ypath), exist_ok=True)
    with open(ypath, "w") as yf:
        yaml.safe_dump({"pmf": {"window_centers": missing}}, yf, sort_keys=False)
    with open(tpath, "w") as tf:
        tf.write(",".join(f"{z:.2f}" for z in missing))

    print(f"Auto-Insert: wrote {ypath} and {tpath} with {len(missing)} midpoint window(s) for replicate {rep['name']}.")
    print("Run additional windows with:")
    print(f"  python3 scripts/universal/08_run_pmf.py {run_dir} --tag {rep['name']} --window-centers-file {ypath} --time-per-window 15 --k 700")

PY