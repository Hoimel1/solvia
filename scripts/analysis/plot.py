import json
import matplotlib.pyplot as plt
import sys, os

# Expect JSON path as first argument
if len(sys.argv) < 2:
    print("Usage: python plot.py <pmf_rep.json> [output.png]")
    sys.exit(1)
json_path = sys.argv[1]
with open(json_path) as f:
    data = json.load(f)

outpath = os.path.splitext(json_path)[0] + ".png"
if len(sys.argv) > 2:
    outpath = sys.argv[2]

z = data["z_centers_nm"]
G = data["G_kJmol"]
Gse = data.get("G_se_kJmol", None)

plt.figure(figsize=(8,5))
plt.errorbar(z, G, yerr=Gse if Gse else None, fmt='-o', markersize=3, linewidth=1)
plt.axhline(0, color='k', linestyle='--')
plt.xlabel("z (nm)")
plt.ylabel("Free Energy G(z) [kJ/mol]")
plt.title(f"PMF Replicate {data['replicate']} (shift {data.get('shift_applied_nm',0):+.2f} nm)")
plt.grid(True)
plt.tight_layout()
plt.savefig(outpath, dpi=300)
print(f"Plot saved to {os.path.abspath(outpath)}")