import numpy as np
import sys
from pathlib import Path

# Ensure repository root is importable
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.hc50 import AdsorptionParams, compute_kp_ads_nm, compute_hc50_uM_from_kp


def make_pmf(z):
    # synthetic adsorption minimum at z=0.6 nm, depth ~ -8 kJ/mol,
    # parabolic well on water side; core side unfavorable
    W = 0.5 * (z - 0.6) ** 2 / 0.1 - 8.0
    W[z < 0] = 5.0 + z[z < 0]
    return W


def test_kp_and_hc50_monotonic():
    z = np.linspace(-1.5, 2.0, 351)
    W = make_pmf(z)
    counts = np.full_like(z, 50, dtype=int)
    p = AdsorptionParams()
    kp, meta = compute_kp_ads_nm(z, W, counts, p)
    assert meta["ok"] and kp > 0.0
    hc = compute_hc50_uM_from_kp(kp, p)
    assert hc["ok"]
    lo, hi = hc["hc50_uM_range"]
    assert 0 < lo < hi
