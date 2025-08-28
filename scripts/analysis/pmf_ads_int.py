# save as scripts/analysis/pmf_ads_int.py
import json, sys, math
j = json.load(open(sys.argv[1]))
z   = j["z_centers_nm"]
G   = j["G_kJmol"]
kT  = float(j["kT_kJmol"])
def mean_exp_minus_beta(G, z, zmin, zmax):
    # simple trapez-Integration über alle z-Bins im Bereich
    sel = [(zi,Gi) for zi,Gi in zip(z,G) if zmin <= zi <= zmax]
    if len(sel)<2: raise SystemExit(f"Not enough points in [{zmin},{zmax}]")
    sel.sort()
    S = 0.0
    for (z0,G0),(z1,G1) in zip(sel, sel[1:]):
        f0 = math.exp(-(G0)/kT); f1 = math.exp(-(G1)/kT)
        S += 0.5*(f0+f1)*(z1-z0)
    return S
surf = mean_exp_minus_beta(G, z, 0.6, 1.2)
bulk = mean_exp_minus_beta(G, z, 2.4, 2.8)
dG_ads_int = -kT*math.log(surf/bulk)
print(f"ΔG_ads,int = {dG_ads_int:.2f} kJ/mol  ({dG_ads_int/kT:.2f} kT)")