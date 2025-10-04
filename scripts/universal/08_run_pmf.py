#!/usr/bin/env python3
"""
SOLVIA Enhanced PMF Runner with Local Patch Reference
Implements all improvements from aenderungen.md for robust PMF calculations
"""

import os
import sys
import yaml
import math
import json
import time
import argparse
import subprocess
import numpy as np
from pathlib import Path
from datetime import datetime
import hashlib
import shutil
import re
import platform
import logging

# Physical constants

KB_KJ_MOL_K = 0.008314462618  # kJ/mol/K

# Prefer SciPy's PPF; fallback to Acklam approximation
try:
    from scipy.stats import norm as _scipy_norm  # type: ignore
    _ppf = _scipy_norm.ppf
except Exception:
    def _ppf(p: float) -> float:
        """Acklam approximation to inverse normal CDF for p in (0,1)."""
        if not (0.0 < p < 1.0):
            raise ValueError("p in (0,1)")
        a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
        b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
             6.680131188771972e+01, -1.328068155288572e+01]
        c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
        d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00]
        plow = 0.02425; phigh = 1 - plow
        import math as _m
        if p < plow:
            q = _m.sqrt(-2.0 * _m.log(p))
            return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                   ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        if p > phigh:
            q = _m.sqrt(-2.0 * _m.log(1.0 - p))
            return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                     ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        q = p - 0.5; r = q*q
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0)

# --- math/helpers for QC and spacing ---
def _dz_max_for_overlap(T_K: float, k_kj_mol_nm2: float, overlap: float) -> float:
    """Given target histogram overlap O and harmonic k, return max allowed Δz."""
    sigma = math.sqrt(KB_KJ_MOL_K * T_K / k_kj_mol_nm2)
    # 1D Gaussian: Δz_max ≈ 2 * sigma * z_alpha with alpha = 1 - O/2
    z = _ppf(1.0 - overlap / 2.0)
    return 2.0 * sigma * z

def _round3(x: float) -> float:
    return float(f"{x:.3f}")

def _dedup_sorted_centers(centers):
    """Round to 3 decimals, drop duplicates, return sorted(desc)."""
    s = sorted({_round3(c) for c in centers}, reverse=True)
    return s

def _js_divergence_from_series(z1: np.ndarray, z2: np.ndarray, bins: int = 100) -> float:
    """Compute Jensen-Shannon divergence between two 1D series via histograms."""
    z1 = np.asarray(z1, dtype=float)
    z2 = np.asarray(z2, dtype=float)
    if z1.size == 0 or z2.size == 0:
        return float('nan')
    zmin = float(min(np.min(z1), np.min(z2)))
    zmax = float(max(np.max(z1), np.max(z2)))
    if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax <= zmin:
        return float('nan')
    hist1, edges = np.histogram(z1, bins=bins, range=(zmin, zmax), density=True)
    hist2, _     = np.histogram(z2, bins=edges, density=True)
    p = np.clip(hist1, 1e-12, None); p = p / p.sum()
    q = np.clip(hist2, 1e-12, None); q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: float((a * np.log(a / b)).sum())
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

def _integrated_autocorr_time(x: np.ndarray, max_lag: int = None) -> float:
    """Compute simple integrated autocorrelation time tau_int."""
    x = np.asarray(x)
    x = x - x.mean()
    n = len(x)
    if n < 2:
        return 0.0
    max_lag = max_lag or min(1000, n // 2)
    var = np.var(x)
    if var <= 0:
        return 0.0
    acsum = 0.0
    for lag in range(1, max_lag+1):
        c = np.dot(x[:-lag], x[lag:]) / (n - lag)
        rho = c / var
        if rho <= 0:
            break
        acsum += rho
    tau_int = 1.0 + 2.0 * acsum
    return tau_int

def _robust_ess(x: np.ndarray, W: int | None = None) -> float:
    """Robuste ESS-Berechnung per FFT-ACF und automatischem Windowing.

    Orientierung an Wolff/Stan: summiere ACF bis Lag ~ 6*tau_int iterativ;
    begrenze auf n/4. Gibt ESS = n / (2*tau_int) konservativ zurück.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4:
        return float(max(1, n))
    x = x - np.mean(x)
    # FFT-basierte ACF
    m = 1 << (2*n-1).bit_length()
    fx = np.fft.rfft(x, m)
    ac = np.fft.irfft(fx * np.conj(fx), m)[:n].real
    ac /= ac[0] if ac[0] != 0 else 1.0
    max_lag = min(n//4, 10000)
    tau_sum = 0.0
    window = 0
    for lag in range(1, max_lag+1):
        rho = float(ac[lag])
        tau_sum += rho
        tau_int = 1.0 + 2.0 * tau_sum
        if lag >= 6.0 * tau_int:
            window = lag
            break
        if W is not None and lag >= int(W):
            window = lag
            break
    if window == 0:
        window = max_lag
    tau_int = 1.0 + 2.0 * float(np.sum(ac[1:window+1]))
    tau_int = max(1.0, tau_int)
    ess = float(n / (2.0 * tau_int))
    return max(1.0, ess)

def _ess_pymbar(x: np.ndarray) -> float:
    """Prefer pymbar statistical inefficiency if available; fallback to robust ESS."""
    try:
        from pymbar.timeseries import statistical_inefficiency  # type: ignore
        x = np.asarray(x, dtype=float)
        if x.size < 2:
            return float(max(1, x.size))
        g = float(statistical_inefficiency(x))
        n = len(x)
        if not np.isfinite(g) or g <= 0:
            return _robust_ess(x)
        return max(1.0, n / max(1.0, g))
    except Exception:
        return _robust_ess(x)

class PMFRunner:
    """Enhanced PMF runner with local patch reference and QC gates"""
    
    def __init__(self, run_dir, config=None):
        self.run_dir = Path(run_dir)
        self.config = config or self.load_pmf_config()
        self.pmf_config = self.config.get('pmf', {})
        self.umbrella_config = self.pmf_config.get('umbrella', {})
        self.qc_config = self.pmf_config.get('qc', {})
        self.temperature_K = (self.config.get('simulation', {}) or {}).get('temperature', 310.0)
        self.min_overlap = float(self.qc_config.get('min_neighbor_overlap', 0.10))
        self.target_overlap = float(self.qc_config.get('target_overlap', 0.20))
        self.auto_extend = bool((self.qc_config.get('auto_extend', True)))
        self.extend_ns = float(self.qc_config.get('extend_ns', 10.0))
        self.max_extend_ns = float(self.qc_config.get('max_extend_ns', 100.0))
        # QC: discard initial transient fraction when computing ESS/half-time
        self.qc_discard_frac = float(self.qc_config.get('discard_fraction', 0.1))
        # QC: optional stride for ESS to mitigate oversampling bias (supports int or "auto")
        self.qc_ess_stride = self.qc_config.get('ess_stride', 1)
        self.strict_qc = bool(self.qc_config.get('strict_mode', False))
        self.strict_failure = False
        # Ensure same reference group across all windows to avoid analysis bias
        self.consistent_reference = bool(self.umbrella_config.get('consistent_reference', True))
        self.common_index_file = None
        self.common_ref_group = None
        self.logger = logging.getLogger("pmf")
        self.logger.setLevel(logging.INFO)
        # Handlers are attached lazily in run_pmf_calculation when pmf_dir is known
        # Round history for QC provenance
        self.round_history = []
        # Cache for expensive QC computations keyed by pullx path metadata
        self._window_qc_cache: dict[tuple, dict] = {}

    # --- Region helpers (bulk/approach/interface) ---
    def _region_for_center(self, z: float) -> str:
        if z >= 2.0:
            return "bulk"
        if z >= 0.6:
            return "approach"
        return "interface"

    # --- PBC helpers (XY) ---
    @staticmethod
    def _pbc_wrap(d: float, L: float) -> float:
        try:
            import numpy as _np
            return d - float(_np.rint(d / L)) * L if L and L > 0 else d
        except Exception:
            return d

    def _pbc_d2_xy(self, x: float, y: float, x0: float, y0: float, Lx: float | None, Ly: float | None) -> float:
        dx = self._pbc_wrap(x - x0, Lx) if Lx else (x - x0)
        dy = self._pbc_wrap(y - y0, Ly) if Ly else (y - y0)
        return float(dx*dx + dy*dy)

    @staticmethod
    def _box_lengths_from_gro(path: Path) -> tuple[float | None, float | None]:
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
            if not lines:
                return (None, None)
            parts = [float(x) for x in lines[-1].strip().split()]
            if len(parts) >= 2:
                return (parts[0], parts[1])
        except Exception:
            pass
        return (None, None)

    def _region_cfg(self, z: float) -> dict:
        """Return region-specific QC thresholds with YAML overrides (deep-merge)."""
        defaults = {
            "bulk":      {"min_ess_frames":150, "min_neighbor_overlap":0.12, "half_energy_tol_kj":3.0, "half_z_tol_sigma":2.5, "min_time_ns":20.0},
            "approach":  {"min_ess_frames":200, "min_neighbor_overlap":0.15, "half_energy_tol_kj":2.0, "half_z_tol_sigma":2.0, "min_time_ns":30.0},
            "interface": {"min_ess_frames":300, "min_neighbor_overlap":0.20, "half_energy_tol_kj":1.0, "half_z_tol_sigma":1.5, "min_time_ns":40.0},
        }
        rt = (self.qc_config.get("region_thresholds") or {})
        import copy
        cfg = copy.deepcopy(defaults)
        for k, v in rt.items():
            if k in cfg and isinstance(v, dict):
                cfg[k].update(v)
        return cfg[self._region_for_center(float(z))]

    def _setup_logging(self, pmf_dir: Path):
        """Attach file + console handlers once we know the output directory."""
        # Avoid duplicate handlers on resume calls
        if any(isinstance(h, logging.FileHandler) for h in self.logger.handlers):
            return
        pmf_dir.mkdir(parents=True, exist_ok=True)
        log_path = pmf_dir / "run.log"
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
        fh = logging.FileHandler(log_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(fmt)
        self.logger.addHandler(fh)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(fmt)
        self.logger.addHandler(ch)
        self.logger.info(f"Logging to {log_path}")
        
    def load_pmf_config(self):
        """Load PMF-specific configuration"""
        pmf_config_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
        if pmf_config_path.exists():
            with open(pmf_config_path, 'r') as f:
                return yaml.safe_load(f)
        else:
            # Fallback to standard config
            std_config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
            if not std_config_path.exists():
                raise FileNotFoundError(f"No PMF config found at {pmf_config_path} or {std_config_path}")
            with open(std_config_path, 'r') as f:
                return yaml.safe_load(f)
    
    def read_gro_atoms(self, gro_path):
        """Parse GRO file and return atom information"""
        atoms = []
        with open(gro_path, "r") as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            return atoms
        
        for line in lines[2:-1]:
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            try:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                atoms.append((resname, atomname, x, y, z))
            except ValueError:
                continue
        
        return atoms
    
    def find_peptide_indices(self, gro_atoms):
        """Find peptide atom indices (1-based)"""
        # Standard amino acids
        AA_RESNAMES = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        
        peptide_idx = []
        for i, (resname, _, _, _, _) in enumerate(gro_atoms):
            if resname.upper() in AA_RESNAMES:
                peptide_idx.append(i + 1)  # 1-based
        
        return peptide_idx

    

    def _compute_pbcatom(self, ndx_path: Path, ref_group_name: str, gro_atoms) -> int:
        """Pick a centrally located atom from the reference group for PBC handling.

        We choose the atom closest (in 3D) to the group's COM.
        Returns absolute atom index (1-based), or 0 on failure (COM).
        """
        try:
            ids = self._parse_ndx_group(ndx_path, ref_group_name)
            if not ids:
                return 0
            # Build array of coords
            coords = []
            for i in ids:
                _, _, x, y, z = gro_atoms[i-1]
                coords.append((x, y, z))
            import numpy as _np
            xyz = _np.array(coords, dtype=float)
            com = xyz.mean(axis=0)
            d2 = ((xyz - com)**2).sum(axis=1)
            j = int(d2.argmin())
            return int(ids[j])
        except Exception:
            return 0

    def _get_po4_atoms(self, gro_atoms):
        """Return list of tuples (index1, x, y, z) for phosphate atoms (PO4/P)."""
        atoms = []
        for i, (_res, atom, x, y, z) in enumerate(gro_atoms):
            if atom in ("PO4", "P"):
                atoms.append((i+1, x, y, z))
        return atoms

    def load_outer_po4_indices(self):
        """Load OuterPO4 atom indices from index_leaflets.ndx if present."""
        ndx_path = self.run_dir / "index_leaflets.ndx"
        indices = []
        if not ndx_path.exists():
            return indices
        group = None
        with open(ndx_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('[') and line.endswith(']'):
                    group = line.strip('[]').strip()
                    continue
                if group == 'OuterPO4':
                    parts = re.split(r"\s+", line)
                    for p in parts:
                        if p.isdigit():
                            indices.append(int(p))
        return indices

    def load_inner_po4_indices(self):
        """Load InnerPO4 atom indices from index_leaflets.ndx if present."""
        ndx_path = self.run_dir / "index_leaflets.ndx"
        indices = []
        if not ndx_path.exists():
            return indices
        group = None
        with open(ndx_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('[') and line.endswith(']'):
                    group = line.strip('[]').strip()
                    continue
                if group == 'InnerPO4':
                    parts = re.split(r"\s+", line)
                    for p in parts:
                        if p.isdigit():
                            indices.append(int(p))
        return indices

    def _split_leaflets(self, gro_atoms):
        """Return (outer_ids, inner_ids) for PO4/P atoms.

        Prefer index_leaflets.ndx; fallback by splitting PO4/P by global z-median.
        """
        outer = set(self.load_outer_po4_indices())
        inner = set(self.load_inner_po4_indices())
        if outer or inner:
            return outer, inner
        # Fallback: z-median split
        po4 = self._get_po4_atoms(gro_atoms)
        if not po4:
            return set(), set()
        zs = [z for (_i, _x, _y, z) in po4]
        zmed = float(np.median(zs))
        outer_ids = {i for (i, _x, _y, z) in po4 if z >= zmed}
        inner_ids = {i for (i, _x, _y, z) in po4 if z < zmed}
        return outer_ids, inner_ids

    def generate_leaflet_index(self, gro_path: Path, out_path: Path = None) -> Path:
        """Create index_leaflets.ndx with OuterPO4/InnerPO4 by global z-median split.

        - Uses PO4/P atoms from the provided GRO file.
        - Writes groups [ OuterPO4 ] and [ InnerPO4 ] in self.run_dir by default.
        """
        out_path = out_path or (self.run_dir / "index_leaflets.ndx")
        atoms = self.read_gro_atoms(gro_path)
        po4 = self._get_po4_atoms(atoms)
        if not po4:
            raise RuntimeError("No PO4/P atoms found to split into leaflets.")
        # Global z-median split to define leaflets consistently
        zs = [z for (_i, _x, _y, z) in po4]
        zmed = float(np.median(zs))
        outer = [i for (i, _x, _y, z) in po4 if z >= zmed]
        inner = [i for (i, _x, _y, z) in po4 if z < zmed]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, 'w') as f:
            f.write("[ OuterPO4 ]\n")
            for k, idx in enumerate(outer):
                if k and (k % 15 == 0):
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n\n[ InnerPO4 ]\n")
            for k, idx in enumerate(inner):
                if k and (k % 15 == 0):
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n")
        self.logger.info(
            f"Generated leaflet index with {len(outer)} OuterPO4 and {len(inner)} InnerPO4 atoms: "
            f"{out_path.relative_to(self.run_dir)}"
        )
        return out_path

    def create_midplane_patch_reference_group(self, gro_atoms, peptide_indices, patch_radius=2.5, min_per_leaflet=10, box_xy: tuple | None = None):
        """Create a local midplane reference group by balancing Outer/Inner PO4 around peptide.

        - Select PO4/P within XY radius of peptide COM.
        - Split by leaflet and take equal counts from each side (trim to min count).
        - If not enough atoms, gradually increase radius up to +2.0 nm.
        Returns list of atom indices.
        """
        # Peptide COM in XY
        pep_xyz = np.array([[x, y, z] for idx, (_res, _a, x, y, z) in enumerate(gro_atoms) if (idx+1) in peptide_indices])
        if pep_xyz.size == 0:
            raise ValueError("No peptide atoms found")
        pep_xy = pep_xyz[:, :2].mean(axis=0)
        Lx = Ly = None
        try:
            if box_xy:
                Lx, Ly = box_xy
        except Exception:
            Lx = Ly = None

        # Split leaflets
        outer_ids, inner_ids = self._split_leaflets(gro_atoms)
        po4 = self._get_po4_atoms(gro_atoms)
        # Build candidates by leaflet
        outer_cand = []
        inner_cand = []
        for i, x, y, z in po4:
            dx = x - pep_xy[0]
            dy = y - pep_xy[1]
            # radius filtered below
            if outer_ids and i in outer_ids:
                outer_cand.append((i, x, y))
            elif inner_ids and i in inner_ids:
                inner_cand.append((i, x, y))
            else:
                # If leaflet sets are empty (no index file), split by z relative to peptide mid z
                pass
        # If index file missing, infer leaflet by z relative to peptide mean z
        if not outer_ids and not inner_ids:
            pep_z = float(pep_xyz[:, 2].mean())
            for i, x, y, z in po4:
                if z >= pep_z:
                    outer_cand.append((i, x, y))
                else:
                    inner_cand.append((i, x, y))

        # Radius growth to satisfy minimum per leaflet
        R0 = patch_radius
        R = R0
        def in_radius(lst, R):
            return [(i, x, y) for (i, x, y) in lst if self._pbc_d2_xy(x, y, pep_xy[0], pep_xy[1], Lx, Ly) <= R*R]
        for _ in range(10):
            o_sel = in_radius(outer_cand, R)
            i_sel = in_radius(inner_cand, R)
            if len(o_sel) >= min_per_leaflet and len(i_sel) >= min_per_leaflet:
                break
            R += 0.5
            if R > R0 + 2.0:
                break
        # Balance counts
        m = min(len(o_sel), len(i_sel))
        if m == 0:
            # Fallback: take whatever available, may bias
            combined = [i for (i, *_xy) in (o_sel + i_sel)]
            return combined
        # Take nearest m from each side by XY distance
        def take_nearest(sel, m):
            sel = sorted(sel, key=lambda t: self._pbc_d2_xy(t[1], t[2], pep_xy[0], pep_xy[1], Lx, Ly))
            return [i for (i, _x, _y) in sel[:m]]
        o_ids = take_nearest(o_sel, m)
        i_ids = take_nearest(i_sel, m)
        return o_ids + i_ids

    def create_patch_reference_group(self, gro_atoms, peptide_indices, patch_radius=2.0, box_xy: tuple | None = None):
        """
        Create local patch reference from phosphates near peptide using OuterPO4 indices
        if available; fallback to all PO4/P atoms.
        """
        # peptide COM (XYZ), project to XY
        pep_xyz = np.array([ [x,y,z] for idx, (_,_,x,y,z) in enumerate(gro_atoms) if (idx+1) in peptide_indices ])
        if pep_xyz.size == 0:
            raise ValueError("No peptide atoms found")
        pep_xy = pep_xyz[:, :2].mean(axis=0)
        Lx = Ly = None
        try:
            if box_xy:
                Lx, Ly = box_xy
        except Exception:
            Lx = Ly = None

        # candidate phosphate indices
        outer_indices = set(self.load_outer_po4_indices())
        candidates = []
        for i, (res, atom, x, y, z) in enumerate(gro_atoms):
            is_phos = atom in ("PO4", "P")
            if not is_phos:
                continue
            if outer_indices and (i+1) not in outer_indices:
                continue
            candidates.append( (i+1, x, y) )

        # distance filter in XY
        patch = []
        R = patch_radius
        while True:
            patch = [idx for idx,x,y in candidates if self._pbc_d2_xy(x, y, pep_xy[0], pep_xy[1], Lx, Ly) <= R*R]
            if len(patch) >= 10 or R >= patch_radius + 2.0:
                break
            R += 0.5
        if len(patch) < 10:
            print(f"Warning: only {len(patch)} phosphates in patch (R={R:.1f} nm)")
        return patch
    
    def create_dynamic_index_file(self, gro_path, output_ndx, window_center=None, ref_mode_override=None):
        """Create index file with (by default) a local patch reference.

        For PMF, we prefer to create this once from a single reference structure
        and reuse it for all windows (see ensure_common_index) so that the
        biased coordinate is defined against the same atoms across windows.
        """
        gro_atoms = self.read_gro_atoms(gro_path)
        Lx, Ly = self._box_lengths_from_gro(Path(gro_path))

        # Find peptide
        peptide_indices = self.find_peptide_indices(gro_atoms)
        if not peptide_indices:
            raise ValueError("No peptide found in structure")

        # Get reference mode (default to global for stability)
        ref_mode = ref_mode_override or self.umbrella_config.get('ref_mode', 'global')

        if ref_mode == 'midplane_local':
            patch_radius = float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.5)))
            reference_indices = self.create_midplane_patch_reference_group(
                gro_atoms, peptide_indices, patch_radius, box_xy=(Lx, Ly)
            )
            ref_group_name = "LocalMidplane"
        elif ref_mode == 'hybrid':
            # Try midplane_local first, validate; fall back to local patch, then global PO4
            ref_group_name = "LocalMidplane"
            patch_radius = float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.5)))
            try:
                reference_indices = self.create_midplane_patch_reference_group(
                    gro_atoms, peptide_indices, patch_radius, box_xy=(Lx, Ly)
                )
            except Exception:
                reference_indices = []
            if not reference_indices or not self._validate_reference_group(gro_atoms, reference_indices, 0.95, ref_group_name, box_xy=(Lx, Ly)):
                ref_group_name = "LocalPatch"
                patch_radius = float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.0)))
                try:
                    reference_indices = self.create_patch_reference_group(
                        gro_atoms, peptide_indices, patch_radius, box_xy=(Lx, Ly)
                    )
                except Exception:
                    reference_indices = []
            if not reference_indices or not self._validate_reference_group(gro_atoms, reference_indices, 0.95, ref_group_name, box_xy=(Lx, Ly)):
                ref_group_name = "GlobalMidplane"
                reference_indices = [i+1 for i, (_r, a, *_xyz) in enumerate(gro_atoms) if a in ("PO4","P")]
        elif ref_mode == 'patch':
            # LOCAL PATCH REFERENCE (key improvement)
            patch_radius = float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.0)))
            reference_indices = self.create_patch_reference_group(
                gro_atoms, peptide_indices, patch_radius, box_xy=(Lx, Ly)
            )
            ref_group_name = "LocalPatch"
        elif ref_mode in ('global', 'po4', 'po4_all', 'upperpo4'):
            # Global midplane reference: all PO4/P atoms across both leaflets
            reference_indices = [i+1 for i, (_r, a, *_xyz) in enumerate(gro_atoms) if a in ("PO4","P")]
            ref_group_name = "GlobalMidplane"
        else:
            reference_indices = [i+1 for i, (_r, a, *_xyz) in enumerate(gro_atoms) if a in ("PO4","P")]
            ref_group_name = "GlobalMidplane"

        # Write index file
        with open(output_ndx, 'w') as f:
            # System group (all atoms)
            f.write("[ System ]\n")
            num_atoms = len(gro_atoms)
            for i in range(1, num_atoms + 1):
                if (i - 1) > 0 and (i - 1) % 15 == 0:
                    f.write("\n")
                f.write(f"{i} ")
            f.write("\n\n")

            # Peptide group
            f.write("[ Peptide ]\n")
            for i, idx in enumerate(peptide_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n\n")

            # Reference group (patch or global)
            f.write(f"[ {ref_group_name} ]\n")
            for i, idx in enumerate(reference_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n")

        print(f"Created index with {len(peptide_indices)} peptide atoms, "
              f"{len(reference_indices)} reference atoms ({ref_group_name})")

        return output_ndx, ref_group_name

    def ensure_common_index(self, gro_path, pmf_dir: Path):
        """Create a single index file to be shared by all windows.

        This avoids reference drift/bias in analysis by ensuring the pull
        coordinate is defined against the same LocalPatch (or UpperPO4) atoms
        for every window.
        """
        if self.common_index_file and Path(self.common_index_file).exists():
            return self.common_index_file, self.common_ref_group
        idx_path = pmf_dir / "common_index.ndx"
        idx_path.parent.mkdir(parents=True, exist_ok=True)
        # If a common index already exists (existing run), reuse it to keep reference consistent
        if idx_path.exists():
            self.common_index_file = str(idx_path)
            # Try to infer a valid reference group from file header with priority
            try:
                with open(idx_path, 'r') as f:
                    txt = f.read()
                groups = re.findall(r"\[\s*([^\]]+)\s*\]", txt)
                # Normalize names (strip)
                groups = [g.strip() for g in groups]
                # Choose the best available, in this priority
                for candidate in ("GlobalMidplane", "LocalMidplane", "LocalPatch", "UpperPO4"):
                    if candidate in groups:
                        self.common_ref_group = candidate
                        break
                else:
                    # Fallback if none found
                    self.common_ref_group = "LocalPatch"
                # Compute a central pbcatom for the reference group (or COM for LocalMidplane)
                gro_atoms = self.read_gro_atoms(gro_path)
                if self.common_ref_group in ('LocalMidplane','GlobalMidplane'):
                    self.common_pbcatom1_index = 0
                else:
                    self.common_pbcatom1_index = self._compute_pbcatom(idx_path, self.common_ref_group, gro_atoms)
            except Exception:
                self.common_ref_group = 'LocalPatch'
            return self.common_index_file, self.common_ref_group
        # Otherwise, create new according to config and validate it
        idx_file, ref_group = self.create_dynamic_index_file(gro_path, idx_path)
        gro_atoms = self.read_gro_atoms(gro_path)
        LxLy = self._box_lengths_from_gro(Path(gro_path))
        ref_ids = self._parse_ndx_group(idx_file, ref_group)
        if (not ref_ids) or (not self._validate_reference_group(gro_atoms, ref_ids, 0.95, ref_group, box_xy=LxLy)):
            idx_file, ref_group = self.create_dynamic_index_file(gro_path, idx_path, ref_mode_override='global')
        self.common_index_file = str(idx_file)
        self.common_ref_group = ref_group
        if self.common_ref_group in ('LocalMidplane','GlobalMidplane'):
            self.common_pbcatom1_index = 0
        else:
            self.common_pbcatom1_index = self._compute_pbcatom(idx_path, self.common_ref_group, gro_atoms)
        return self.common_index_file, self.common_ref_group

    def _parse_ndx_group(self, ndx_path: Path, group_name: str):
        group = None
        idxs = []
        with open(ndx_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('[') and line.endswith(']'):
                    group = line.strip('[]').strip()
                    continue
                if group == group_name:
                    parts = re.split(r"\s+", line)
                    for p in parts:
                        if p.isdigit():
                            idxs.append(int(p))
        return idxs

    def _preflight_reference_index(self, gro_path: Path, pmf_dir: Path):
        """Validate the presence and health of the common index and reference group.

        Ensures:
        - common_index.ndx exists
        - reference group exists and has sufficient atoms
        - pbcatom selection is resolvable (or COM for midplane groups)
        """
        if not self.consistent_reference:
            return
        if not self.common_index_file or not Path(self.common_index_file).exists():
            raise RuntimeError("Preflight: common index missing; cannot proceed")
        idx_path = Path(self.common_index_file)
        ref_group = self.common_ref_group or 'GlobalMidplane'
        ids = self._parse_ndx_group(idx_path, ref_group)
        if not ids or len(ids) < 5:
            raise RuntimeError(f"Preflight: reference group '{ref_group}' missing/too small in {idx_path}")
        # Compactness validation (XY for midplane groups)
        gro_atoms = self.read_gro_atoms(gro_path)
        LxLy = self._box_lengths_from_gro(Path(gro_path))
        if not self._validate_reference_group(gro_atoms, ids, 0.95, ref_group, box_xy=LxLy):
            raise RuntimeError(f"Preflight: reference group '{ref_group}' failed compactness validation")
        # Attempt to compute pbcatom unless using a midplane group
        if ref_group not in ('LocalMidplane','GlobalMidplane'):
            gro_atoms = self.read_gro_atoms(gro_path)
            pb = int(self._compute_pbcatom(idx_path, ref_group, gro_atoms))
            if pb <= 0:
                raise RuntimeError("Preflight: could not determine pbcatom for reference group; aborting")

    def _validate_reference_group(self, gro_atoms, ref_indices, threshold: float = 0.95, group_name: str | None = None, box_xy: tuple | None = None) -> bool:
        """Validate reference group compactness with midplane-aware logic.

        - For bimodal z distributions (e.g., LocalMidplane) or when explicitly
          named 'LocalMidplane', evaluate compactness in XY only.
        - Otherwise evaluate 3D compactness around the mean.
        Uses IQR-based outlier fraction and compares against threshold.
        """
        try:
            if not ref_indices or len(ref_indices) < 5:
                return False
            coords = np.array([gro_atoms[i-1][2:5] for i in ref_indices], dtype=float)
            z = coords[:, 2]
            z_spread = float(np.std(z))
            bimodal_like = (z_spread > 0.5) or (group_name in ('LocalMidplane','GlobalMidplane'))
            if bimodal_like:
                xy = coords[:, :2]
                # PBC-aware XY distances around COM if box provided
                Lx = Ly = None
                try:
                    if box_xy:
                        Lx, Ly = box_xy
                except Exception:
                    Lx = Ly = None
                com = xy.mean(axis=0)
                d = np.array([np.sqrt(self._pbc_d2_xy(p[0], p[1], com[0], com[1], Lx, Ly)) for p in xy], dtype=float)
            else:
                d = np.linalg.norm(coords - coords.mean(axis=0), axis=1)
            q1, q3 = np.percentile(d, [25, 75])
            iqr = q3 - q1 if q3 > q1 else (np.std(d) * 1.349)
            outliers = np.sum((d < q1 - 1.5*iqr) | (d > q3 + 1.5*iqr))
            return (outliers / len(d)) < (1.0 - float(threshold))
        except Exception:
            return True

    def robust_pbc_distance(self, pos1: np.ndarray, pos2: np.ndarray, box: np.ndarray | None, pbc_type: str = 'orthorhombic') -> float:
        """Robust PBC distance between two positions.

        - Orthorhombic: minimal image using lengths (Lx,Ly,Lz)
        - Triclinic: provide full 3x3 box matrix H; uses fractional wrapping
        - Else: try MDAnalysis; if unavailable, raise ValueError (no Euclidean fallback)
        """
        if box is None:
            raise ValueError("PBC distance requested without box.")
        p = np.asarray(pos1, dtype=float).reshape(3)
        q = np.asarray(pos2, dtype=float).reshape(3)
        B = np.asarray(box, dtype=float)
        if pbc_type == 'orthorhombic' and B.ndim == 1 and B.size >= 3:
            d = q - p
            Lx, Ly, Lz = float(B[0]), float(B[1]), float(B[2])
            import numpy as _np
            if Lx > 0: d[0] -= float(_np.rint(d[0] / Lx)) * Lx
            if Ly > 0: d[1] -= float(_np.rint(d[1] / Ly)) * Ly
            if Lz > 0: d[2] -= float(_np.rint(d[2] / Lz)) * Lz
            return float(np.linalg.norm(d))
        if B.shape == (3, 3):
            H = B
            Hinv = np.linalg.inv(H)
            dfrac = Hinv.dot((q - p).reshape(3))
            dfrac -= np.round(dfrac)
            dcart = H.dot(dfrac)
            return float(np.linalg.norm(dcart))
        # Try MDAnalysis as a last resort for unusual PBC configs
        try:
            from MDAnalysis.lib.distances import distance_array  # type: ignore
            return float(distance_array(p.reshape(1, 3), q.reshape(1, 3), box=B)[0, 0])
        except Exception:
            raise ValueError("Unsupported PBC configuration; provide orthorhombic lengths or 3x3 box matrix.")

    def _current_z_from_gro(self, gro_path: Path, ndx_path: Path, ref_group_name: str):
        """Compute z = COM(Peptide) - COM(RefGroup) with robust PBC handling along z."""
        atoms = self.read_gro_atoms(gro_path)
        pep = self.find_peptide_indices(atoms)
        ref = self._parse_ndx_group(ndx_path, ref_group_name)
        if not pep or not ref:
            return 0.0
        # Parse box length in z from last line of GRO
        try:
            with open(gro_path, 'r') as f:
                lines = f.readlines()
            parts = (lines[-1].strip().split() if lines else [])
            zbox = float(parts[2]) if len(parts) >= 3 else None
            box = np.array([float(parts[i]) for i in range(3)]) if len(parts) >= 3 else None
        except Exception:
            zbox = None; box = None
        import numpy as _np
        pep_z = _np.mean([atoms[i-1][4] for i in pep])
        ref_z = _np.mean([atoms[i-1][4] for i in ref])
        dz = float(pep_z - ref_z)
        if zbox and zbox > 0:
            import numpy as _np
            dz = dz - float(_np.rint(dz / zbox)) * zbox
        return dz

    def run_smd_ladder(self, start_structure: Path, target_z: float, output_dir: Path,
                        index_file: Path, ref_group: str, pbcatom1: int, seed: int):
        """Run short SMD to nudge structure toward target_z and return final gro."""
        smd_dir = output_dir / "smd_prep"
        smd_dir.mkdir(exist_ok=True)

        # Determine initial z
        z0 = self._current_z_from_gro(start_structure, index_file, ref_group)
        # Cap SMD time to prevent long start-up; skip if already near target
        # Allow longer SMD pre-pull by default; 1.0 ns cap for CG is reasonable
        max_smd = float(self.umbrella_config.get('smd_max_ns', 1.0))
        smd_time = min(float(self.umbrella_config.get('pre_smd_ns', 1.0)), max_smd)
        if abs(float(target_z) - z0) <= 0.05 or smd_time <= 0.0:
            return start_structure
        # Compute rate in nm/ns and convert to nm/ps for GROMACS
        rate_nm_per_ns = (float(target_z) - z0) / smd_time
        # For CG, rates up to 1.5–2.0 nm/ns are acceptable for short nudges
        max_rate_nm_per_ns = float(self.umbrella_config.get('smd_max_rate_nm_per_ns', 1.5))
        if abs(rate_nm_per_ns) > max_rate_nm_per_ns:
            rate_nm_per_ns = max_rate_nm_per_ns if rate_nm_per_ns > 0 else -max_rate_nm_per_ns
        rate = rate_nm_per_ns / 1000.0  # nm/ps expected by GROMACS
        print(f"  SMD prep: z0={z0:+.3f} → {target_z:+.3f} nm over {smd_time:.2f} ns (rate {rate_nm_per_ns:.3f} nm/ns)")

        # Write MDP
        smd_mdp = smd_dir / "smd.mdp"
        mdp_content = f"""
; SMD pulling for window preparation (conservative settings)
integrator              = sd
dt                      = 0.01
nsteps                  = {int(max(1, smd_time * 1000 / 0.01))}
nstcomm                 = 100

; Velocity init for short SMD prep (avoid 0 K start)
gen-vel                 = yes
gen-temp                = {int(self.temperature_K)}
gen-seed                = {seed}
continuation            = no
ld-seed                 = {seed}

; Output (sparse)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstenergy               = 500
nstxout-compressed      = 500

; Constraints
constraints             = none
constraint-algorithm    = lincs

; Thermostat/Barostat
tcoupl                  = no
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = {int(self.temperature_K)}
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau-p                   = 12.0
ref-p                   = 1.0 1.0
compressibility         = 4.5e-5 0

; Nonbonded / PBC
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1

; Pull code - SMD
pull                    = yes
pull-ngroups            = 2
pull-ncoords            = 1
pull-group1-name        = {ref_group}
pull-group2-name        = Peptide
pull-group1-pbcatom     = {pbcatom1}
pull-group2-pbcatom     = 0
pull-pbc-ref-prev-step-com = yes
pull-coord1-type        = umbrella
pull-coord1-geometry    = distance
pull-coord1-dim         = N N Y
pull-coord1-groups      = 1 2
pull-coord1-rate        = {rate}
pull-coord1-k           = {int(self.umbrella_config.get('smd_force_constant', 800))}
pull-coord1-init        = 0.0
pull-coord1-start       = yes
"""
        with open(smd_mdp, 'w') as f:
            f.write(mdp_content)

        # Prepare container paths
        project_root = Path(__file__).parent.parent.parent
        rel_index = os.path.relpath(index_file, project_root)
        rel_mdp = os.path.relpath(smd_mdp, project_root)
        rel_start = os.path.relpath(start_structure, project_root)
        smd_tpr = smd_dir / "smd.tpr"
        rel_tpr = os.path.relpath(smd_tpr, project_root)
        rel_smd = os.path.relpath(smd_dir / "smd", project_root)
        rel_wd = os.path.relpath(smd_dir, project_root)
        container_wd = f"/work/{rel_wd}"

        # grompp
        grompp_cmd = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "grompp",
            "-f", f"/work/{rel_mdp}",
            "-c", f"/work/{rel_start}",
            "-p", f"/work/{os.path.relpath(self.find_topology_file(), project_root)}",
            "-n", f"/work/{rel_index}",
            "-o", f"/work/{rel_tpr}",
            "-maxwarn", "2"
        ]
        r = subprocess.run(grompp_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=str(project_root))
        if r.returncode != 0:
            print("  ✗ SMD grompp failed:\n" + r.stdout)
            return start_structure

        # mdrun
        mdrun_cmd = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "mdrun",
            "-deffnm", f"/work/{rel_smd}",
            "-s", f"/work/{rel_tpr}",
            "-px", f"/work/{os.path.relpath(smd_dir / 'smd_pullx.xvg', project_root)}",
            "-pf", f"/work/{os.path.relpath(smd_dir / 'smd_pullf.xvg', project_root)}"
        ]
        r2 = subprocess.run(mdrun_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=str(project_root))
        if r2.returncode != 0:
            print("  ✗ SMD mdrun failed:\n" + r2.stdout)
            return start_structure
        final_gro = smd_dir / "smd.gro"
        return final_gro if final_gro.exists() else start_structure
    
    def calculate_window_positions(self):
        """Return umbrella window centers (nm) based on evidence/default strategy.

        Evidence one-sided (default for RBC-like): +3.0 → 0.0 nm
        - Bulk (Δz=0.2): 3.0, 2.8, 2.6, 2.4, 2.2
        - Approach (Δz=0.1): 2.0, 1.9, …, 1.2
        - Interface/core (Δz=0.05–0.1): 1.15, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80,
                                         0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.00
        """
        # If explicit list provided in config, honor it
        explicit = self.umbrella_config.get('window_centers') or self.umbrella_config.get('windows')
        if explicit:
            return _dedup_sorted_centers(explicit)

        strategy = self.umbrella_config.get('strategy', 'evidence_one_sided')

        if strategy == 'evidence_one_sided':
            windows = [
                # Bulk
                3.0, 2.8, 2.6, 2.4, 2.2,
                # Approach (0.1)
                2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2,
                # Interface/core
                1.15, 1.10, 1.05, 1.00, 0.95, 0.90, 0.85, 0.80,
                0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.00,
            ]
            return _dedup_sorted_centers(windows)

        # Fallback to legacy adaptive two-sided spacing
        z_range = self.umbrella_config.get('z_range', {})
        windows = []
        bulk_max = z_range.get('bulk_max', 2.8)
        bulk_min = z_range.get('bulk_min', 2.4)
        windows.extend([bulk_max, bulk_min])
        current = bulk_min - z_range.get('coarse_step', 0.2)
        dense_max = z_range.get('dense_max', 0.6)
        while current > dense_max:
            windows.append(current)
            current -= z_range.get('coarse_step', 0.2)
        dense_min = z_range.get('dense_min', -0.6)
        dense_step = z_range.get('dense_step', 0.15)
        current = dense_max
        while current >= dense_min:
            windows.append(current)
            current -= dense_step
        end_z = z_range.get('end_z', -2.0)
        current = dense_min - z_range.get('coarse_step', 0.2)
        while current >= end_z:
            windows.append(current)
            current -= z_range.get('coarse_step', 0.2)
        windows = _dedup_sorted_centers(windows)
        # Auto densify based on target overlap
        k_global = float(self.umbrella_config.get('force_constant', 900))
        dz_max = _dz_max_for_overlap(float(self.temperature_K), k_global, float(self.target_overlap))
        densified = [windows[0]]
        for c in windows[1:]:
            prev = densified[-1]
            gap = abs(prev - c)
            if gap > dz_max:
                n_mid = int(math.ceil(gap / dz_max)) - 1
                step = gap / (n_mid + 1)
                for j in range(1, n_mid+1):
                    densified.append(_round3(prev - math.copysign(j*step, prev - c)))
            densified.append(c)
        windows = _dedup_sorted_centers(densified)
        return windows
    
    def check_window_overlap(self, window1_pullx, window2_pullx, target_overlap=0.20):
        """Check overlap between adjacent windows from pullx using configured method.

        Methods: 'hist' (default) with FD bins; 'kde' using Gaussian KDE.
        """
        method = str(self.qc_config.get('overlap_method', 'hist')).lower()
        try:
            d1 = np.loadtxt(window1_pullx, comments=['#', '@'])
            d2 = np.loadtxt(window2_pullx, comments=['#', '@'])
        except Exception:
            return 0.0
        if d1.size == 0 or d2.size == 0:
            return 0.0
        z1 = d1[:, 1]; z2 = d2[:, 1]
        if method == 'kde':
            try:
                from scipy import stats
                kde1 = stats.gaussian_kde(z1, bw_method='scott')
                kde2 = stats.gaussian_kde(z2, bw_method='scott')
                x_min = min(z1.min(), z2.min()) - 2*np.std(z1)
                x_max = max(z1.max(), z2.max()) + 2*np.std(z2)
                xs = np.linspace(x_min, x_max, 1024)
                p1 = kde1(xs); p2 = kde2(xs)
                overlap = float(np.trapz(np.minimum(p1, p2), x=xs))
                return max(0.0, min(1.0, overlap))
            except Exception:
                pass  # fall back to histogram
        # Histogram method
        zmin = float(min(z1.min(), z2.min())); zmax = float(max(z1.max(), z2.max()))
        if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax <= zmin:
            return 0.0
        def _bins(z):
            n = len(z)
            if n < 5:
                return 10
            iqr = np.subtract(*np.percentile(z, [75, 25]))
            bw = 2 * iqr / np.cbrt(n) if iqr > 0 else (z.std() * 3.49 / np.cbrt(n))
            B = max(10, int((zmax - zmin) / max(bw, 1e-4)))
            return min(200, B)
        B = max(_bins(z1), _bins(z2))
        bins = np.linspace(zmin, zmax, B+1)
        h1, _ = np.histogram(z1, bins=bins, density=True)
        h2, _ = np.histogram(z2, bins=bins, density=True)
        bw = (zmax - zmin) / B
        overlap = float(np.minimum(h1, h2).sum() * bw)
        return max(0.0, min(1.0, overlap))

    # --- QC helpers with caching -------------------------------------------------

    def _qc_cache_key(self, pullx_path: str | Path) -> tuple:
        """Return a stat-based cache key for a pullx file."""
        try:
            st = os.stat(pullx_path)
            return (str(pullx_path), st.st_size, st.st_mtime_ns)
        except FileNotFoundError:
            return (str(pullx_path), None, None)

    def _get_window_qc(self, pullx_path: str, k_kj_mol_nm2: float, center_nm: float,
                        half_energy_tol_kj: float) -> dict:
        """Return QC metrics for a window, using a cached result when possible."""
        cache_key = (self._qc_cache_key(pullx_path),
                     round(float(k_kj_mol_nm2), 6),
                     round(float(center_nm), 6),
                     round(float(half_energy_tol_kj), 6))
        cached = self._window_qc_cache.get(cache_key)
        if cached is not None:
            return cached
        qc = self.compute_window_qc(pullx_path, k_kj_mol_nm2, center_nm, half_energy_tol_kj)
        self._window_qc_cache[cache_key] = qc
        return qc

    def _invalidate_qc_cache(self, pullx_path: str | Path) -> None:
        """Drop cached entries referencing a given pullx file (used after reruns)."""
        pullx_path = str(pullx_path)
        keys_to_drop = [key for key in self._window_qc_cache if key[0][0] == pullx_path]
        for key in keys_to_drop:
            self._window_qc_cache.pop(key, None)

    def compute_window_qc(self, pullx_path: str, k_kj_mol_nm2: float, center_nm: float,
                           half_energy_tol_kj: float = 2.0):
        """QC with telemetry: ESS, ESS-rate, half-energy, Half-Z, drift, tau_int, sigmas."""
        try:
            data = np.loadtxt(pullx_path, comments=['#','@'])
        except Exception:
            return {"ess": None, "ess_pass": False, "half_energy_diff_kj": None, "half_pass": False,
                    "drift_nm_per_ns": None, "drift_pass": False}
        if data.size == 0:
            return {"ess": None, "ess_pass": False, "half_energy_diff_kj": None, "half_pass": False,
                    "drift_nm_per_ns": None, "drift_pass": False}

        t = data[:, 0]  # ps
        z = data[:, 1]

        # Discard initial transient
        if self.qc_discard_frac > 0 and len(z) > 10:
            start_idx = int(len(z) * self.qc_discard_frac)
            start_idx = min(max(start_idx, 0), len(z) - 2)
            z = z[start_idx:]
            t = t[start_idx:]

        # dt and segment time
        dt_ps = float(np.median(np.diff(t))) if len(t) >= 2 else 1.0
        time_ns = float((t[-1] - t[0]) * 0.001) if len(t) >= 2 else (len(t) * dt_ps * 0.001)

        # Optional: ESS stride (supports "auto")
        stride_cfg = self.qc_config.get('ess_stride', 1)
        if isinstance(stride_cfg, str) and str(stride_cfg).lower() == "auto":
            tau_est = _integrated_autocorr_time(z)
            target_dt_ps = max(4.0, 0.2 * tau_est * dt_ps)
            s = max(1, int(round(target_dt_ps / dt_ps)))
            z = z[::s]; t = t[::s]
            dt_ps *= s
        else:
            try:
                s = int(stride_cfg)
            except Exception:
                s = 1
            if s > 1:
                z = z[::s]; t = t[::s]
                dt_ps *= s

        # ESS and tau_int
        tau_int = _integrated_autocorr_time(z)
        ess = int(_ess_pymbar(z))
        ess_rate = (ess / time_ns) if time_ns > 0 else 0.0

        # Half-energy halves
        mid = len(z) // 2
        def mean_bias(zseg):
            dz = zseg - center_nm
            return 0.5 * k_kj_mol_nm2 * float(np.mean(dz * dz))
        e1 = mean_bias(z[:mid]); e2 = mean_bias(z[mid:])
        de = abs(e1 - e2)

        # Half-Z test on means (use effective sample size to account for correlation)
        import math as _m
        m1 = float(np.mean(z[:mid])) if mid > 1 else float(np.mean(z))
        m2 = float(np.mean(z[mid:])) if len(z) - mid > 1 else float(np.mean(z))
        s1 = float(np.var(z[:mid], ddof=1)) if mid > 1 else 0.0
        s2 = float(np.var(z[mid:], ddof=1)) if len(z) - mid > 1 else 0.0
        # Effective sample sizes per half (more conservative and robust)
        try:
            n1_eff = int(max(1, _ess_pymbar(z[:mid])))
            n2_eff = int(max(1, _ess_pymbar(z[mid:])))
        except Exception:
            n1_eff = max(1, mid)
            n2_eff = max(1, len(z) - mid)
        se = _m.sqrt((s1 / n1_eff) + (s2 / n2_eff)) if (n1_eff > 1 and n2_eff > 1 and (s1 > 0 or s2 > 0)) else float('inf')
        half_z_stat = abs((m1 - m2) / se) if _m.isfinite(se) and se > 0 else 0.0

        # Drift
        try:
            tt_ns = t * 0.001
            A = np.vstack([tt_ns, np.ones_like(tt_ns)]).T
            slope, intercept = np.linalg.lstsq(A, z, rcond=None)[0]
            drift = float(abs(slope))
        except Exception:
            drift = None
        drift_tol = float(self.qc_config.get('drift_tolerance_nm_per_ns', 0.1))

        # Pass/Fail using global defaults (region-specific handled in run_qc_checks)
        ess_pass = ess >= int(self.qc_config.get('min_ess_frames', 200))
        half_pass_energy = de <= half_energy_tol_kj
        half_pass_z = half_z_stat <= float(self.qc_config.get('half_z_tol_sigma', 2.0))
        half_pass = (half_pass_energy and half_pass_z)

        sigma_meas = float(np.std(z))
        sigma_theory = float(_m.sqrt(max(1e-12, KB_KJ_MOL_K * float(self.temperature_K) / float(k_kj_mol_nm2))))

        return {
            "ess": ess, "ess_pass": ess_pass,
            "half_energy_diff_kj": de, "half_pass": half_pass,
            "drift_nm_per_ns": drift, "drift_pass": (drift is not None) and (drift <= drift_tol),
            "time_ns": time_ns, "ess_rate": ess_rate,
            "tau_int": float(tau_int),
            "half_z_stat": float(half_z_stat),
            "sigma_meas": sigma_meas, "sigma_theory": sigma_theory
        }
    
    def run_umbrella_window(self, window_center, window_dir, start_structure, 
                          replicate=1, extend_time=0,
                          prod_time_override: float | None = None,
                          force_k_override: float | None = None,
                          ref_mode_override: str | None = None,
                          seed_tag: str | None = None,
                          force_gen_vel: bool = False):
        """Run single umbrella window with all improvements"""
        
        window_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare index for this window: reuse common index if configured
        index_file = window_dir / "index.ndx"
        if (self.consistent_reference and self.common_index_file) and not ref_mode_override:
            try:
                shutil.copy(self.common_index_file, index_file)
                ref_group = self.common_ref_group or 'LocalPatch'
            except Exception:
                # Fallback to generating a window-specific index
                index_file, ref_group = self.create_dynamic_index_file(start_structure, index_file, window_center, ref_mode_override=ref_mode_override)
        else:
            index_file, ref_group = self.create_dynamic_index_file(
                start_structure,
                index_file,
                window_center,
                ref_mode_override=ref_mode_override
            )

        # Determine pbcatom for reference group
        using_common_index = bool(self.consistent_reference and self.common_index_file and not ref_mode_override)
        if ref_group in ('LocalMidplane','GlobalMidplane'):
            pbcatom1 = 0
        elif using_common_index:
            pbcatom1 = int(self.common_pbcatom1_index or 0)
        else:
            gro_atoms = self.read_gro_atoms(start_structure)
            pbcatom1 = int(self._compute_pbcatom(index_file, ref_group, gro_atoms))
        # Validate reference group and pbcatom selection early to avoid unstable pulls
        ref_ids = self._parse_ndx_group(index_file, ref_group)
        if not ref_ids:
            raise RuntimeError(f"Reference group '{ref_group}' is empty or missing in index {index_file}.")
        if pbcatom1 <= 0 and ref_group not in ('LocalMidplane','GlobalMidplane'):
            raise RuntimeError(
                "No valid pbcatom for reference group. Ensure the reference group contains atoms "
                "and that the common index contains the selected reference group."
            )

        
        # Generate deterministic seed
        peptide_id = self.run_dir.name.split('_')[0]
        seed_suffix = f"_{seed_tag}" if seed_tag else ""
        seed_string = f"{peptide_id}_rep{replicate}_z{window_center:.2f}{seed_suffix}"
        seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16) % 1000000
        
        # Production time
        prod_time = (prod_time_override if prod_time_override is not None else self.umbrella_config.get('production_ns', 60.0)) + extend_time
        force_k = float(force_k_override if force_k_override is not None else self.umbrella_config.get('force_constant', 900))
        # Optional adaptive force constant targeting desired sigma_z
        adapt_k_cfg = (self.umbrella_config.get('adaptive_k') or {})
        # Apply adaptation only if no explicit override is provided
        if bool(adapt_k_cfg.get('enabled', False)) and (force_k_override is None):
            sigma_target = float(adapt_k_cfg.get('sigma_target_nm', 0.06))
            k_min = float(adapt_k_cfg.get('k_min', 200.0))
            k_max = float(adapt_k_cfg.get('k_max', 5000.0))
            k_thermal = KB_KJ_MOL_K * float(self.temperature_K) / max(1e-6, sigma_target**2)
            k_opt = float(np.clip(k_thermal, k_min, k_max))
            if abs(k_opt - force_k) / max(1e-12, force_k) > 0.1:
                self.logger.info(f"  Adaptive k: {force_k:.1f} → {k_opt:.1f} kJ/mol/nm^2 (target σ≈{sigma_target} nm)")
            force_k = k_opt
        elif bool(adapt_k_cfg.get('enabled', False)) and (force_k_override is not None):
            self.logger.info(f"  Adaptive k: override {force_k_override:.1f} provided → skip adaptation")
        
        # Resolve structure for grompp
        structure_for_grompp = Path(start_structure)
        prev_gro = window_dir / "umbrella.gro"
        if extend_time > 0:
            # Prefer checkpoint continuation; if not available, reuse last produced structure
            if (window_dir / "umbrella.cpt").exists():
                pass  # handled via -cpi below
            elif prev_gro.exists():
                structure_for_grompp = prev_gro

        # Optionally prepare start via short SMD toward target center (reduces initial forces)
        # Skip SMD if we are extending from an existing checkpoint
        smd_prepare = bool(self.umbrella_config.get('smd_prepare', True))
        cpt_path = window_dir / "umbrella.cpt"
        if smd_prepare and float(self.umbrella_config.get('pre_smd_ns', 1.0)) > 0 and not (extend_time > 0 and cpt_path.exists()):
            smd_gro = self.run_smd_ladder(Path(start_structure), float(window_center), window_dir, index_file, ref_group, pbcatom1, seed)
            if smd_gro and Path(smd_gro).exists():
                structure_for_grompp = Path(smd_gro)

        # (moved) Reference validation is performed after optional pre-equilibration for robustness

        # --- Optional pre-equilibration at target z to damp initial forces ---
        try:
            pre_eq_ns = float(self.umbrella_config.get('pre_equil_ns', 0.0) or 0.0)
        except Exception:
            pre_eq_ns = 0.0
        do_pre_eq = (pre_eq_ns > 0.0) and not (extend_time > 0)
        if do_pre_eq:
            preeq_mdp = window_dir / "preeq.mdp"
            preeq_tpr = window_dir / "preeq.tpr"
            preeq_deffnm = window_dir / "preeq"
            # Output stride similar to production
            try:
                pull_out_ps = float(self.umbrella_config.get('pull_output_ps', 1.0))
            except Exception:
                pull_out_ps = 1.0
            out_every = max(1, int(round(pull_out_ps / 0.02)))
            preeq_mdp_content = f"""
; Short pre-equilibration at target z
integrator              = md
dt                      = 0.02
nsteps                  = {int(pre_eq_ns * 1000 / 0.02)}
nstcomm                 = 100

; Output (sparse)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstenergy               = 500
nstxout-compressed      = 500

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW (Martini 3)
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1

; Thermostat
tcoupl                  = v-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = {int(self.temperature_K)}

; Pressure coupling (softer)
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau-p                   = 12.0
ref-p                   = 1.0 1.0
compressibility         = 4.5e-5 0

; Continuation (no new velocities)
gen-vel                 = no
continuation            = yes

; Constraints
constraints             = none
constraint-algorithm    = lincs

; Pull code – same geometry as production
pull                    = yes
pull-ngroups            = 2
pull-ncoords            = 1
pull-group1-name        = {ref_group}
pull-group2-name        = Peptide
pull-group1-pbcatom     = {pbcatom1}
pull-group2-pbcatom     = 0
pull-pbc-ref-prev-step-com = yes
pull-coord1-type        = umbrella
pull-coord1-geometry    = distance
pull-coord1-dim         = N N Y
pull-coord1-groups      = 1 2
pull-coord1-init        = {window_center}
pull-coord1-k           = {int(force_k)}
pull-coord1-start       = no

; Pull output
pull-nstxout            = {out_every}
pull-nstfout            = {out_every}
"""
            with open(preeq_mdp, 'w') as f:
                f.write(preeq_mdp_content)

            # Build and run pre-equilibration TPR/MD
            project_root = Path(__file__).parent.parent.parent
            top_file = self.find_topology_file()
            rel_preeq_mdp = os.path.relpath(preeq_mdp, project_root)
            rel_structure = os.path.relpath(structure_for_grompp, project_root)
            rel_top = os.path.relpath(top_file, project_root)
            rel_index = os.path.relpath(index_file, project_root)
            rel_preeq_tpr = os.path.relpath(preeq_tpr, project_root)
            rel_preeq_deffnm = os.path.relpath(preeq_deffnm, project_root)
            rel_window_dir = os.path.relpath(window_dir, project_root)
            container_wd = f"/work/{rel_window_dir}"
            grompp_preeq = [
                "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
                "grompp",
                "-f", f"/work/{rel_preeq_mdp}",
                "-c", f"/work/{rel_structure}",
                "-p", f"/work/{rel_top}",
                "-n", f"/work/{rel_index}",
                "-o", f"/work/{rel_preeq_tpr}",
                "-maxwarn", "2"
            ]
            res_preq = subprocess.run(grompp_preeq, capture_output=True, text=True, cwd=str(project_root))
            if res_preq.returncode != 0:
                # Fallback for older GROMACS lacking pbcatom=0 support
                if ref_group in ('LocalMidplane','GlobalMidplane') and pbcatom1 == 0:
                    try:
                        gro_atoms_fb = self.read_gro_atoms(structure_for_grompp)
                        pb_fallback = int(self._compute_pbcatom(index_file, ref_group, gro_atoms_fb))
                        mdp_fb = preeq_mdp_content.replace(f"pull-group1-pbcatom     = {pbcatom1}", f"pull-group1-pbcatom     = {pb_fallback}")
                        with open(preeq_mdp, 'w') as f:
                            f.write(mdp_fb)
                        res_preq = subprocess.run(grompp_preeq, capture_output=True, text=True, cwd=str(project_root))
                    except Exception:
                        pass
                # Additional fallback: replace group2 pbcatom COM (0) with central peptide atom
                if res_preq.returncode != 0:
                    try:
                        gro_atoms_fb = self.read_gro_atoms(structure_for_grompp)
                        pb2 = int(self._compute_pbcatom(index_file, 'Peptide', gro_atoms_fb))
                        with open(preeq_mdp, 'r') as f:
                            txt = f.read()
                        txt2 = txt.replace("pull-group2-pbcatom     = 0", f"pull-group2-pbcatom     = {pb2}")
                        if txt2 != txt:
                            with open(preeq_mdp, 'w') as f:
                                f.write(txt2)
                            res_preq = subprocess.run(grompp_preeq, capture_output=True, text=True, cwd=str(project_root))
                    except Exception:
                        pass
                if res_preq.returncode != 0:
                    print(f"  ✗ Pre-equil grompp failed: {res_preq.stderr}")
            else:
                mdrun_preeq = [
                    "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
                    "mdrun",
                    "-deffnm", f"/work/{rel_preeq_deffnm}",
                    "-s", f"/work/{rel_preeq_tpr}",
                ]
                # Ensure preeq.gro is written explicitly
                rel_preeq_conf = os.path.relpath(window_dir / "preeq.gro", project_root)
                mdrun_preeq += ["-c", f"/work/{rel_preeq_conf}"]
                print(f"  Running pre-equilibration for {pre_eq_ns:.2f} ns …")
                res_mdrun_preeq = subprocess.run(mdrun_preeq, capture_output=True, text=True, cwd=str(project_root))
                if res_mdrun_preeq.returncode != 0:
                    print(f"  ✗ Pre-equil mdrun failed: {res_mdrun_preeq.stderr}")
                else:
                    preeq_gro = window_dir / "preeq.gro"
                    if preeq_gro.exists():
                        structure_for_grompp = preeq_gro
                        print("  ✓ Pre-equilibration done; using preeq.gro for production.")

        # Validate reference group after optional pre-equilibration (PBC-aware XY)
        try:
            gro_atoms_check = self.read_gro_atoms(structure_for_grompp)
            LxLy = self._box_lengths_from_gro(Path(structure_for_grompp))
            if not self._validate_reference_group(gro_atoms_check, ref_ids, threshold=0.95, group_name=ref_group, box_xy=LxLy):
                self.logger.warning("Reference group validation failed (outliers). Consider a more global PO4 reference.")
                # Optional: in hybrid mode, fallback to global PO4 only if not enforcing consistent reference
                if (self.umbrella_config.get('ref_mode', 'patch') == 'hybrid') and (not self.consistent_reference):
                    index_file, ref_group = self.create_dynamic_index_file(structure_for_grompp, index_file, ref_mode_override='global')
                    # Re-parse ids and recompute pbcatom for the switched group
                    ref_ids = self._parse_ndx_group(index_file, ref_group)
                    pbcatom1 = 0 if ref_group in ('LocalMidplane','GlobalMidplane') else int(self._compute_pbcatom(index_file, ref_group, gro_atoms_check))
        except Exception:
            pass

        # Create umbrella MDP
        mdp_file = window_dir / "umbrella.mdp"
        # Determine safe output stride and velocity/continuation settings
        try:
            pull_out_ps = float(self.umbrella_config.get('pull_output_ps', 1.0))
        except Exception:
            pull_out_ps = 1.0
        out_every = max(1, int(round(pull_out_ps / 0.02)))

        def _gro_has_velocities(path: Path) -> bool:
            try:
                with open(path, 'r') as fh:
                    lines = fh.readlines()
                count = 0; with_vel = 0
                for line in lines[2:-1][:50]:
                    if len(line) < 68:
                        continue
                    count += 1
                    try:
                        float(line[44:52]); float(line[52:60]); float(line[60:68])
                        with_vel += 1
                    except Exception:
                        pass
                return with_vel > 0 and with_vel >= max(1, int(0.2 * count))
            except Exception:
                return False

        has_vel = _gro_has_velocities(structure_for_grompp)
        gen_vel_flag = True if force_gen_vel else ((not has_vel) and (extend_time <= 0))
        continuation_flag = (not gen_vel_flag)

        mdp_content = f"""
; Umbrella sampling window at z = {window_center:.3f} nm
integrator              = md
dt                      = 0.02
nsteps                  = {int(prod_time * 1000 / 0.02)}
nstcomm                 = 100

; Output control  
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 1000
nstenergy              = 1000
nstxout-compressed     = 1000

; Neighbor searching
cutoff-scheme          = Verlet
nstlist                = 20
pbc                    = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype            = reaction-field
rcoulomb               = 1.1
epsilon_r              = 15
vdw-type               = Cut-off
vdw-modifier           = Potential-shift-verlet
rvdw                   = 1.1

; Temperature coupling
tcoupl                 = v-rescale
tc-grps                = System
tau-t                  = 1.0
ref-t                  = 310

; Pressure coupling
pcoupl                 = Parrinello-Rahman
pcoupltype            = semiisotropic
tau-p                  = 12.0
ref-p                  = 1.0 1.0
compressibility       = 4.5e-5 0

; Velocity generation
gen-vel                = {"yes" if gen_vel_flag else "no"}
continuation          = {"yes" if continuation_flag else "no"}

; Constraints
constraints            = none
constraint-algorithm   = lincs

        ; Pull code - Umbrella
        pull                   = yes
        pull-ngroups          = 2
        pull-ncoords          = 1
        pull-group1-name      = {ref_group}
        pull-group2-name      = Peptide
        pull-group1-pbcatom   = {pbcatom1}
        pull-group2-pbcatom   = 0  ; Use COM as reference
        pull-pbc-ref-prev-step-com = yes
        pull-coord1-type      = umbrella
        pull-coord1-geometry  = distance
        pull-coord1-dim       = N N Y
        pull-coord1-groups    = 1 2
        pull-coord1-init      = {window_center}
pull-coord1-k         = {int(force_k)}
pull-coord1-start     = no  ; Use init value directly

; Output pull force and coordinate
pull-nstxout          = {out_every}
pull-nstfout          = {out_every}
"""
        # Optional lipid restraints for production (disabled by default)
        if bool(self.umbrella_config.get('lipid_z_posres', False)):
            mdp_content += """

; Optional lipid restraints
define                 = -DPOSRES_LIPID_Z
"""
        
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        
        # Find topology
        top_file = self.find_topology_file()
        
        # Prepare TPR for run (convert-tpr for extensions with checkpoint; else grompp)
        tpr_file = window_dir / "umbrella.tpr"
        
        # Get relative paths from project root
        project_root = Path(__file__).parent.parent.parent
        rel_mdp = os.path.relpath(mdp_file, project_root)
        rel_structure = os.path.relpath(structure_for_grompp, project_root)
        rel_top = os.path.relpath(top_file, project_root)
        rel_index = os.path.relpath(index_file, project_root)
        rel_tpr = os.path.relpath(tpr_file, project_root)
        
        # Container workdir so emergency PDBs land in the window dir
        rel_window_dir = os.path.relpath(window_dir, project_root)
        container_wd = f"/work/{rel_window_dir}"
        grompp_cmd = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "grompp",  # Docker already provides "gmx"
            "-f", f"/work/{rel_mdp}",
            "-c", f"/work/{rel_structure}",
            "-p", f"/work/{rel_top}",
            "-n", f"/work/{rel_index}",
            "-o", f"/work/{rel_tpr}",
            "-maxwarn", "2"
        ]
        
        # If extending and checkpoint + original TPR exist, extend TPR via convert-tpr (avoid step mismatch)
        tpr_ready = False
        if extend_time > 0 and cpt_path.exists() and tpr_file.exists():
            try:
                rel_tpr = os.path.relpath(tpr_file, project_root)
                tpr_ext = window_dir / "umbrella_ext.tpr"
                rel_tpr_ext = os.path.relpath(tpr_ext, project_root)
                extend_ps = int(round(float(extend_time) * 1000.0))
                # Aim for absolute end time >= current last written time + extend_ps
                # Prefer pullx.xvg as source of current time
                pullx_path = window_dir / "pullx.xvg"
                until_ps = None
                if pullx_path.exists():
                    try:
                        last_t = None
                        with open(pullx_path, 'r') as fh:
                            for ln in fh:
                                if ln.startswith('#') or ln.startswith('@'):
                                    continue
                                parts = ln.split()
                                if len(parts) >= 2:
                                    try:
                                        last_t = float(parts[0])
                                    except Exception:
                                        pass
                        if last_t is not None:
                            until_ps = int(round(last_t + extend_ps))
                    except Exception:
                        until_ps = None
                convert_cmd = [
                    "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
                    "convert-tpr",
                    "-s", f"/work/{rel_tpr}",
                ]
                if until_ps is not None and until_ps > 0:
                    convert_cmd += ["-until", str(until_ps)]
                else:
                    convert_cmd += ["-extend", str(extend_ps)]
                convert_cmd += ["-o", f"/work/{rel_tpr_ext}"]
                print(f"  Extending TPR by +{extend_time:.1f} ns via convert-tpr…")
                rconv = subprocess.run(convert_cmd, capture_output=True, text=True, cwd=str(project_root))
                if rconv.returncode != 0:
                    print(f"  ✗ convert-tpr failed; falling back to grompp: {rconv.stderr}")
                else:
                    try:
                        # Replace original TPR with extended one
                        os.replace(tpr_ext, tpr_file)
                        tpr_ready = True
                        print("  ✓ TPR extended")
                    except Exception as e:
                        print(f"  ✗ Failed to replace TPR after convert-tpr: {e}")
            except Exception as e:
                print(f"  ✗ convert-tpr error: {e}")
        if not tpr_ready:
            print(f"  Running grompp…")
            sys.stdout.flush()
            result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
            if result.returncode != 0:
                # Fallback: if Local/Global Midplane with pbcatom1=0 is not supported by old GROMACS, retry with central atom
                if ref_group in ('LocalMidplane','GlobalMidplane') and pbcatom1 == 0:
                    try:
                        gro_atoms_fb = self.read_gro_atoms(structure_for_grompp)
                        pb_fallback = int(self._compute_pbcatom(index_file, ref_group, gro_atoms_fb))
                        # Regenerate MDP with fallback pbcatom
                        mdp_content_fb = mdp_content.replace(f"pull-group1-pbcatom   = {pbcatom1}", f"pull-group1-pbcatom   = {pb_fallback}")
                        with open(mdp_file, 'w') as f:
                            f.write(mdp_content_fb)
                        result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
                    except Exception:
                        pass
                # Additional fallback: replace group2 pbcatom COM (0) with central peptide atom if still failing
                if result.returncode != 0:
                    try:
                        gro_atoms_fb = self.read_gro_atoms(structure_for_grompp)
                        pb2 = int(self._compute_pbcatom(index_file, 'Peptide', gro_atoms_fb))
                        with open(mdp_file, 'r') as f:
                            txt = f.read()
                        txt2 = txt.replace("pull-group2-pbcatom   = 0", f"pull-group2-pbcatom   = {pb2}")
                        if txt2 != txt:
                            with open(mdp_file, 'w') as f:
                                f.write(txt2)
                            result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
                    except Exception:
                        pass
                if result.returncode != 0:
                    print(f"  ✗ Error in grompp: {result.stderr}")
                    return None
            print(f"  ✓ TPR file created")
        
        # Run mdrun via Docker
        log_file = window_dir / "umbrella.log"
        pullf_file = window_dir / "pullf.xvg"
        pullx_file = window_dir / "pullx.xvg"
        
        # Get relative paths from project root
        rel_tpr = os.path.relpath(tpr_file, project_root)
        rel_pullx = os.path.relpath(pullx_file, project_root)
        rel_pullf = os.path.relpath(pullf_file, project_root)
        rel_log = os.path.relpath(log_file, project_root)
        rel_deffnm = os.path.relpath(window_dir / "umbrella", project_root)

        # Build mdrun command with optional continuation
        mdrun_cmd = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "mdrun",
            "-deffnm", f"/work/{rel_deffnm}",
            "-s", f"/work/{rel_tpr}",
            "-px", f"/work/{rel_pullx}",
            "-pf", f"/work/{rel_pullf}",
            "-g", f"/work/{rel_log}"
        ]
        # Force confout path to guarantee umbrella.gro availability even with append/resume
        rel_confout = os.path.relpath(window_dir / "umbrella.gro", project_root)
        mdrun_cmd += ["-c", f"/work/{rel_confout}"]
        # if extending, use checkpoint input and append
        rel_cpt = os.path.relpath(cpt_path, project_root)
        if extend_time > 0 and cpt_path.exists():
            mdrun_cmd += ["-cpi", f"/work/{rel_cpt}", "-append"]

        print(f"  Starting mdrun for z={window_center:.3f} nm...")
        print(f"  Production time: {prod_time} ns")
        sys.stdout.flush()
        
        # Run with real-time output streaming
        process = subprocess.Popen(
            mdrun_cmd, 
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        
        # Stream output line by line
        last_time = time.time()
        lincs_flag = False
        for line in process.stdout:
            # Show progress lines (step updates)
            if "step" in line.lower():
                # Clear line and print progress
                sys.stdout.write(f"\r  Progress: {line.strip()[:80]}")
                sys.stdout.flush()
            # Show important messages
            elif "Back Off!" in line:
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                continue
            elif any(keyword in line for keyword in ["WARNING", "ERROR", "Note", "Writing"]):
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                if ("LINCS WARNING" in line) or ("Fatal error" in line) or ("Constraint error" in line):
                    lincs_flag = True
                    try:
                        process.terminate()
                    except Exception:
                        pass
                    try:
                        process.terminate()
                    except Exception:
                        pass
            # Show periodic heartbeat
            current_time = time.time()
            if current_time - last_time > 5:  # Every 5 seconds
                sys.stdout.write(".")
                sys.stdout.flush()
                last_time = current_time
        
        process.wait()
        sys.stdout.write("\n")  # New line after progress
        
        if process.returncode == 0:
            self._invalidate_qc_cache(pullx_file)
            print(f"✓ Window z={window_center:.3f} nm completed")
            return {
                "center": window_center,
                "pullf": str(pullf_file),
                "pullx": str(pullx_file),
                "tpr": str(tpr_file),
                "seed": seed,
                "k": float(force_k)
            }
        # Attempt auto-recovery on LINCS or failure: short stabilization run with softer settings
        print(f"✗ Window z={window_center:.3f} nm encountered instability (LINCS or failure). Attempting stabilization…")
        stab_mdp = window_dir / "stabilize.mdp"
        stab_time_ns = 0.05  # 50 ps
        k_soft = max(300, int(float(force_k) * 0.5))
        stab_mdp_content = f"""
; Short stabilization with softer pull and smaller timestep
integrator              = sd
dt                      = 0.01
nsteps                  = {int(stab_time_ns * 1000 / 0.01)}
nstcomm                 = 100

; Output (minimal)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 200
nstenergy               = 200
nstxout-compressed      = 0

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW (MARTINI defaults)
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1

; Thermostat (sd integrates stochastic thermostat)
tcoupl                  = no
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 310
ld-seed                 = {seed}

; Pressure coupling (more forgiving)
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau-p                   = 12.0
ref-p                   = 1.0 1.0
compressibility         = 4.5e-5 0

; Constraints
constraints             = none
constraint-algorithm    = lincs

        ; Pull code - Umbrella (soft)
        pull                    = yes
        pull-ngroups            = 2
        pull-ncoords            = 1
        pull-group1-name        = {ref_group}
        pull-group2-name        = Peptide
        pull-group1-pbcatom     = {pbcatom1}
        pull-group2-pbcatom     = 0
        pull-pbc-ref-prev-step-com = yes
        pull-coord1-type        = umbrella
        pull-coord1-geometry    = distance
        pull-coord1-dim         = N N Y
        pull-coord1-groups      = 1 2
        pull-coord1-init        = {window_center}
pull-coord1-k           = {k_soft}
pull-coord1-start       = yes
"""
        with open(stab_mdp, 'w') as f:
            f.write(stab_mdp_content)
        stab_tpr = window_dir / "stabilize.tpr"
        stab_deffnm = window_dir / "stabilize"
        rel_stab_mdp = os.path.relpath(stab_mdp, project_root)
        rel_stab_tpr = os.path.relpath(stab_tpr, project_root)
        rel_stab_deffnm = os.path.relpath(stab_deffnm, project_root)

        # Use the best available structure for stabilization
        stab_structure = prev_gro if prev_gro.exists() else structure_for_grompp
        rel_stab_structure = os.path.relpath(stab_structure, project_root)

        stab_grompp = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "grompp",
            "-f", f"/work/{rel_stab_mdp}",
            "-c", f"/work/{rel_stab_structure}",
            "-p", f"/work/{rel_top}",
            "-n", f"/work/{rel_index}",
            "-o", f"/work/{rel_stab_tpr}",
            "-maxwarn", "2"
        ]
        print("  Running stabilization grompp…")
        result = subprocess.run(stab_grompp, capture_output=True, text=True, cwd=str(project_root))
        if result.returncode != 0:
            print(f"  ✗ Stabilization grompp failed: {result.stderr}")
            print(f"✗ Window z={window_center:.3f} nm failed")
            return None
        print("  ✓ Stabilization TPR created")
        stab_mdrun = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "mdrun",
            "-deffnm", f"/work/{rel_stab_deffnm}",
            "-s", f"/work/{rel_stab_tpr}"
        ]
        # Ensure stabilized confout is explicitly written
        rel_stab_conf = os.path.relpath(window_dir / "stabilize.gro", project_root)
        stab_mdrun += ["-c", f"/work/{rel_stab_conf}"]
        print("  Starting stabilization mdrun…")
        res2 = subprocess.run(stab_mdrun, capture_output=True, text=True, cwd=str(project_root))
        if res2.returncode != 0:
            print(f"  ✗ Stabilization mdrun failed: {res2.stderr}")
            print(f"✗ Window z={window_center:.3f} nm failed")
            return None
        print("  ✓ Stabilization completed. Restarting production…")
        # Use stabilized structure as new input
        stabilized_gro = window_dir / "stabilize.gro"
        if stabilized_gro.exists():
            rel_structure = os.path.relpath(stabilized_gro, project_root)
        else:
            # fallback to previous gro
            rel_structure = os.path.relpath(prev_gro if prev_gro.exists() else structure_for_grompp, project_root)
        # Rebuild production TPR
        grompp_cmd = [
            "docker", "compose", "run", "--rm", "--workdir", container_wd, "gromacs",
            "grompp",
            "-f", f"/work/{rel_mdp}",
            "-c", f"/work/{rel_structure}",
            "-p", f"/work/{rel_top}",
            "-n", f"/work/{rel_index}",
            "-o", f"/work/{rel_tpr}",
            "-maxwarn", "2"
        ]
        res3 = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
        if res3.returncode != 0:
            print(f"  ✗ Re-grompp after stabilization failed: {res3.stderr}")
            print(f"✗ Window z={window_center:.3f} nm failed")
            return None
        # Rerun production
        process = subprocess.Popen(
            mdrun_cmd,
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        last_time = time.time()
        lincs_flag = False
        for line in process.stdout:
            if "step" in line.lower():
                sys.stdout.write(f"\r  Progress: {line.strip()[:80]}")
                sys.stdout.flush()
            elif "Back Off!" in line:
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                continue
            elif any(keyword in line for keyword in ["WARNING", "ERROR", "Note", "Writing"]):
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                if ("LINCS WARNING" in line) or ("Fatal error" in line) or ("Constraint error" in line):
                    lincs_flag = True
            current_time = time.time()
            if current_time - last_time > 5:
                sys.stdout.write(".")
                sys.stdout.flush()
                last_time = current_time
        process.wait()
        sys.stdout.write("\n")
        if process.returncode == 0:
            self._invalidate_qc_cache(pullx_file)
            print(f"✓ Window z={window_center:.3f} nm completed (after stabilization)")
            return {
                "center": window_center,
                "pullf": str(pullf_file),
                "pullx": str(pullx_file),
                "tpr": str(tpr_file),
                "seed": seed,
                "k": float(force_k)
            }
        print(f"✗ Window z={window_center:.3f} nm failed even after stabilization")
        # Mark this window to avoid endless retries on resume
        try:
            (window_dir / ".no_retry").write_text("unstable_after_stabilization\n")
        except Exception:
            pass
        return None
    def _extend_failed_windows(self, pmf_dir, windows_data, qc_results, start_structure, replicate):
        """Adaptive extension for ESS/half-time deficits with region-aware thresholds.

        Returns (windows_data, n_extended, extended_entries)
        where extended_entries is a list of {center, add_ns}.
        """
        if not self.auto_extend:
            return windows_data, 0, []

        extended = 0
        extended_entries = []
        new_windows = []

        # Quick maps for current QC
        ess_map = {float(it['center']): it for it in qc_results.get('ess_check', [])}
        half_map = {float(it['center']): it for it in qc_results.get('convergence_check', [])}

        for w in windows_data:
            if not w:
                new_windows.append(w)
                continue
            c = float(w['center'])
            e = ess_map.get(c, {})
            h = half_map.get(c, {})
            reg = self._region_cfg(c)

            need_ess = int(reg["min_ess_frames"])
            cur_ess = int(e.get("ess", 0))
            ess_rate = float(e.get("ess_rate") or 0.0)
            time_ns = float(e.get("time_ns") or 0.0)
            ess_fail = (cur_ess < need_ess)
            half_fail = (not bool(h.get("passed", h.get("half_pass", False))))

            if not (ess_fail or half_fail):
                new_windows.append(w)
                continue

            # Compute required extra time
            add_by_ess = 0.0
            if ess_fail:
                deficit = max(0, need_ess - cur_ess)
                if ess_rate <= 1e-6:
                    add_by_ess = self.extend_ns  # stagnating → minimal block
                else:
                    add_by_ess = max(5.0, float(np.ceil(deficit / ess_rate)))

            add_by_half = 0.0
            if half_fail and time_ns < float(reg.get("min_time_ns", 0.0)):
                add_by_half = float(reg["min_time_ns"]) - time_ns

            add = max(add_by_ess, add_by_half, 0.0)
            if add <= 0.0:
                new_windows.append(w)
                continue

            # Cap per window using marker file
            window_dir = pmf_dir / "windows" / f"z_{c:+.3f}"
            marker = window_dir / ".extended_ns"
            already = 0.0
            if marker.exists():
                try:
                    already = float(marker.read_text().strip())
                except Exception:
                    already = 0.0
            if already >= self.max_extend_ns:
                print(f"  ↷ Window z={c:+.3f} reached max extension ({already} ns), not extending further.")
                new_windows.append(w)
                continue

            add = float(min(add, self.max_extend_ns - already))
            # Enforce global time budget per window if configured
            try:
                max_time_budget = float(self.qc_config.get('max_time_per_window', 0.0) or 0.0)
            except Exception:
                max_time_budget = 0.0
            if max_time_budget > 0.0:
                base_prod = float(self.umbrella_config.get('production_ns', 60.0))
                if (base_prod + already + add) > max_time_budget:
                    add = max(0.0, max_time_budget - base_prod - already)
            print(f"  ⤴ Extending window z={c:+.3f} by +{add:.1f} ns (ESS/half-time).")
            res = self.run_umbrella_window(c, window_dir, start_structure, replicate, extend_time=add)
            if res:
                new_windows.append(res)
                extended += 1
                extended_entries.append({"center": float(c), "add_ns": float(add)})
                try:
                    marker.write_text(f"{already + add}\n")
                except Exception:
                    pass
            else:
                new_windows.append(w)

        if extended:
            print(f"Extended {extended} windows due to ESS/half-time QC failures.")
        return new_windows, extended, extended_entries

    def _replicate_low_ess_windows(self, pmf_dir: Path, windows_data, qc_results, start_structure, replicate, targets_override=None):
        """Fallback: replicate low-ESS or half-failing windows with new seeds.

        Strategy (configurable in pmf.qc.fallback):
          - Identify windows with ESS below trigger or half_pass/drift fail
          - For each such center, run N short replicates (fresh gen-vel, new seed)
          - Optionally increase k slightly to reduce relaxation time
          - Optionally force a global PO4 reference to stabilize z
          - For windows with half_fail, optionally restart from nearest passing neighbor .gro
        Appends new window dicts and returns (windows_data, n_new, replicated_entries).
        replicated_entries is a list of dicts with details for round history.
        """
        fb = (self.qc_config.get('fallback') or {})
        if not bool(fb.get('enabled', True)):
            return windows_data, 0, []
        min_ess_trig = int(fb.get('min_ess_trigger', 50))
        n_reps = int(fb.get('replicate_count', 3))
        rep_ns = float(fb.get('per_replicate_ns', 50.0))
        k_scale = float(fb.get('k_tune_scale', 1.2))
        ref_override = str(fb.get('ref_mode', 'global')).lower()  # 'global' -> all PO4/P
        restart_from_neighbor = bool(fb.get('restart_from_neighbor', True))

        # Build maps for quick lookup
        ess_map = {int(_round3(item['center'])*1000): int(item.get('ess', 0)) for item in qc_results.get('ess_check', [])}
        half_map = {int(_round3(item['center'])*1000): bool(item.get('half_pass', False)) for item in qc_results.get('convergence_check', [])}

        # Select target centers to replicate
        targets = []
        allowed = None
        if targets_override is not None:
            try:
                allowed = {_round3(float(t)) for t in targets_override}
            except Exception:
                allowed = {_round3(float(t)) for t in list(targets_override)}
        for w in sorted([w for w in windows_data if w], key=lambda x: x['center'], reverse=True):
            key = int(_round3(w['center'])*1000)
            ess = ess_map.get(key, 0)
            half_ok = half_map.get(key, True)
            if allowed is not None and _round3(float(w['center'])) not in allowed:
                continue
            if (ess < min_ess_trig) or (not half_ok):
                targets.append(float(w['center']))
        targets = _dedup_sorted_centers(targets)
        if not targets:
            return windows_data, 0, []

        self.logger.warning(f"QC fallback: replicating {len(targets)} low-ESS/half-failing windows with {n_reps}×{rep_ns:.0f} ns")

        # Helper: find neighbor gro for restart
        def _neighbor_gro(center: float):
            wins_sorted = sorted([w for w in windows_data if w], key=lambda x: x['center'], reverse=True)
            # find nearest neighbor with half_pass True
            cands = []
            for w in wins_sorted:
                k = int(_round3(w['center'])*1000)
                if half_map.get(k, True):
                    # prefer neighbors within 0.1–0.2 nm
                    cands.append((abs(float(w['center']) - float(center)), w))
            if not cands:
                return None
            cands.sort(key=lambda x: x[0])
            best = cands[0][1]
            gro = Path(best['pullx']).parent / 'umbrella.gro'
            return gro if gro.exists() else None

        base_k = float(self.umbrella_config.get('force_constant', 900))
        k_rep = base_k * k_scale if k_scale and k_scale > 0 else base_k
        # Map ref override to our create_dynamic_index modes
        ref_mode_override = None
        if ref_override in ('global', 'po4_all', 'upperpo4', 'po4'):
            ref_mode_override = 'global'
        elif ref_override in ('patch', 'local'):
            ref_mode_override = 'patch'
        elif ref_override in ('hybrid', 'midplane', 'midplane_local'):
            ref_mode_override = 'hybrid'

        new_entries = 0
        replicated_entries = []
        for c in targets:
            win_dir_base = pmf_dir / 'windows' / f"z_{float(c):+0.3f}"
            # Choose start structure
            start_for_rep = None
            if restart_from_neighbor:
                try:
                    start_for_rep = _neighbor_gro(float(c))
                except Exception:
                    start_for_rep = None
            if start_for_rep is None or not Path(start_for_rep).exists():
                # fallback to the window's own last structure if present, else global start_structure
                maybe = win_dir_base / 'umbrella.gro'
                start_for_rep = maybe if maybe.exists() else start_structure

            for r in range(1, n_reps+1):
                rep_dir = win_dir_base / f"rep_{r}"
                seed_tag = f"fb{r}"
                try:
                    res = self.run_umbrella_window(
                        c, rep_dir, start_for_rep, replicate=replicate,
                        prod_time_override=rep_ns,
                        force_k_override=k_rep,
                        ref_mode_override=ref_mode_override,
                        seed_tag=seed_tag,
                        force_gen_vel=True,
                    )
                except Exception as e:
                    self.logger.error(f"Fallback replicate failed for z={c:+.3f} rep {r}: {e}")
                    res = None
                if res:
                    windows_data.append(res)
                    new_entries += 1
                    # Reason tagging for provenance
                    reasons = []
                    key = int(_round3(float(c))*1000)
                    if ess_map.get(key, 10**9) < min_ess_trig:
                        reasons.append("low_ess")
                    if half_map.get(key, True) is False:
                        reasons.append("half_fail")
                    replicated_entries.append({
                        "center": float(c),
                        "rep": int(r),
                        "ns": float(rep_ns),
                        "k": float(k_rep),
                        "ref_mode": (ref_override or "default"),
                        "seed_tag": seed_tag,
                        "reason": reasons or ["fallback"],
                    })
        if new_entries:
            self.logger.info(f"QC fallback added {new_entries} replicate window runs.")
        return windows_data, new_entries, replicated_entries

    def _qc_gates_passed(self, qc_results):
        """Return booleans for whether QC gates pass.

        Currently enforces:
        - Overlap gate: all adjacent pairs pass min_neighbor_overlap
        - ESS gate: all windows meet min_ess_frames

        Returns (overlap_ok, ess_ok).
        """
        def _active(entries):
            return [it for it in (entries or []) if not it.get('superseded', False)]

        ol_entries = qc_results.get('overlap_check', [])
        ess_entries = _active(qc_results.get('ess_check', []))
        ol = all(item.get('passed', False) for item in ol_entries) if ol_entries else True
        ess = all(item.get('passed', False) for item in ess_entries) if ess_entries else True
        require_half = bool(self.qc_config.get('require_half_convergence', False))
        if require_half:
            half_entries = _active(qc_results.get('convergence_check', []))
            half_ok = all(entry.get('passed', entry.get('half_pass', False)) for entry in half_entries) if half_entries else True
            ess = ess and half_ok
        return ol, ess

    def _collect_problem_centers(self, qc_results, include_half=True):
        centers = []
        for entry in qc_results.get('ess_check', []) or []:
            try:
                if not entry.get('passed', False):
                    centers.append(float(entry.get('center')))
            except Exception:
                continue
        if include_half:
            for entry in qc_results.get('convergence_check', []) or []:
                try:
                    if not entry.get('passed', entry.get('half_pass', False)):
                        centers.append(float(entry.get('center')))
                except Exception:
                    continue
        return _dedup_sorted_centers(centers)
    
    def find_topology_file(self):
        """Find the appropriate topology file"""
        # Prefer replicate-specific PMF topology if available
        try:
            rep = getattr(self, 'replicate', None)
        except Exception:
            rep = None
        if rep is not None:
            specific = self.run_dir / "pmf_system" / f"replicate_{rep}" / "system.top"
            if specific.exists():
                return specific

        # Prefer explicit system/*.top next (stable, non-ambiguous)
        system_tops = list(self.run_dir.glob("system/*.top"))
        if system_tops:
            for top in system_tops:
                if "n1" in top.name:
                    return top
            for top in system_tops:
                if top.name == "system.top":
                    return top
            return sorted(system_tops)[0]

        # Else choose newest pmf_system replicate if multiple exist
        pmf_dirs = list(self.run_dir.glob("pmf_system/replicate_*/system.top"))
        if pmf_dirs:
            try:
                pmf_dirs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
            except Exception:
                pmf_dirs.sort()
            return pmf_dirs[0]

        # Fallback to other locations
        possible_paths = [
            self.run_dir / "system" / "system.top",
            self.run_dir / "equilibration" / "system.top"
        ]
        
        for path in possible_paths:
            if path.exists():
                return path
        
        # Print what we searched for debugging
        print("Error: No topology file found. Searched:")
        print(f"  - pmf_system/replicate_*/system.top")
        print(f"  - system/*.top")
        for path in possible_paths:
            print(f"  - {path.relative_to(self.run_dir)}")
        
        raise FileNotFoundError("No topology file found")
    
    def run_qc_checks(self, windows_data):
        qc_results = {"overlap_check": [], "convergence_check": [], "ess_check": []}
        val_windows = sorted([w for w in windows_data if w], key=lambda x: x['center'], reverse=True)

        # Overlap checks with region-strict thresholds and expected-vs-measured
        k_global = float(self.umbrella_config.get('force_constant', 900))
        for i in range(len(val_windows) - 1):
            w1, w2 = val_windows[i], val_windows[i+1]
            if float(w1['center']) == float(w2['center']):
                continue
            overlap = self.check_window_overlap(w1['pullx'], w2['pullx'], self.target_overlap)
            cfg1 = self._region_cfg(w1['center']); cfg2 = self._region_cfg(w2['center'])
            min_ol = max(float(cfg1["min_neighbor_overlap"]), float(cfg2["min_neighbor_overlap"]))
            dz = abs(float(w1['center']) - float(w2['center']))
            k1 = float(w1.get('k', k_global)); k2 = float(w2.get('k', k_global))
            expected = self._expected_overlap_gaussian_mixed(dz, k1, k2)
            # JS divergence as additional diagnostic (histogram-based)
            jsd = None
            try:
                d1 = np.loadtxt(w1['pullx'], comments=['#','@'])
                d2 = np.loadtxt(w2['pullx'], comments=['#','@'])
                if d1.size and d2.size:
                    z1 = d1[:,1]; z2 = d2[:,1]
                    jsd = float(_js_divergence_from_series(z1, z2, bins=100))
            except Exception:
                jsd = None
            qc_results["overlap_check"].append({
                "windows": [w1['center'], w2['center']],
                "overlap": overlap,
                "expected": expected,
                "delta": float(overlap - expected),
                "min_required": min_ol,
                "js_divergence": jsd,
                "passed": overlap >= min_ol
            })

        # Per-window checks with region thresholds
        for w in val_windows:
            reg = self._region_cfg(w['center'])
            k_w = float(w.get('k', k_global))
            qc = self._get_window_qc(w['pullx'], k_w, w['center'], reg["half_energy_tol_kj"])
            qc_results["ess_check"].append({
                "center": w['center'],
                "pullx": str(w['pullx']),
                "ess": qc['ess'],
                "ess_rate": qc.get('ess_rate'),
                "time_ns": qc.get('time_ns'),
                "min_required": int(reg["min_ess_frames"]),
                "passed": int(qc['ess'] or 0) >= int(reg["min_ess_frames"])
            })
            # Region-aware half-time decision
            half_pass_reg = (float(qc['half_energy_diff_kj']) <= float(reg["half_energy_tol_kj"])) and \
                            (float(qc.get('half_z_stat', 0.0)) <= float(reg["half_z_tol_sigma"]))
            conv_item = {
                "center": w['center'],
                "half_energy_diff_kj": qc['half_energy_diff_kj'],
                "half_z_stat": qc.get('half_z_stat'),
                "z_tol": float(reg["half_z_tol_sigma"]),
                "half_pass": bool(half_pass_reg),
                "drift_nm_per_ns": qc.get('drift_nm_per_ns'),
                "drift_pass": qc.get('drift_pass', False),
                "pullx": str(w['pullx']),
            }
            conv_item["passed"] = bool(conv_item["half_pass"]) and bool(conv_item.get('drift_pass', True))
            qc_results["convergence_check"].append(conv_item)
        return qc_results

    def write_run_manifest(self, pmf_dir: Path, window_centers, qc_results):
        """Write a RUN_INFO.yaml manifest with provenance for scientific reproducibility."""
        manifest = {}
        # Timestamp
        manifest['timestamp'] = datetime.now().isoformat()
        # Git info
        project_root = Path(__file__).parent.parent.parent
        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(project_root), text=True).strip()
            branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=str(project_root), text=True).strip()
            status = subprocess.check_output(["git", "status", "--porcelain"], cwd=str(project_root), text=True)
            dirty = bool(status.strip())
            manifest['git'] = {"commit": commit, "branch": branch, "dirty": dirty}
        except Exception:
            manifest['git'] = {"commit": None, "branch": None, "dirty": None}
        # Containers (if available)
        containers = None
        cfg_path = project_root / "config" / "config.yaml"
        if cfg_path.exists():
            try:
                with open(cfg_path, 'r') as f:
                    cfg = yaml.safe_load(f)
                    containers = cfg.get('containers', None)
            except Exception:
                containers = None
        manifest['containers'] = containers
        # Environment
        manifest['environment'] = {
            'python': sys.version.split()[0],
            'platform': platform.platform(),
            'numpy': np.__version__,
        }
        # PMF config snapshot
        manifest['pmf_config_snapshot'] = self.config.get('pmf', {})
        # Reference info
        manifest['reference'] = {
            'consistent_reference': bool(self.consistent_reference),
            'reference_group': self.common_ref_group,
            'index_file': str(Path(self.common_index_file).relative_to(self.run_dir)) if self.common_index_file else None,
        }
        # Peptide footprint proxy (if available)
        try:
            if getattr(self, '_peptide_rg_xy', None) is not None:
                manifest['peptide'] = {'rg_xy_nm': float(self._peptide_rg_xy)}
        except Exception:
            pass
        # Windows and QC summary
        manifest['windows'] = {
            'centers': window_centers,
            'n_windows': len(window_centers),
        }
        # Record effective patch radius used
        try:
            manifest['reference']['patch_radius_nm'] = float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.0)))
        except Exception:
            manifest['reference']['patch_radius_nm'] = float(self.umbrella_config.get('patch_radius', 2.0))
        # Persist
        run_info = pmf_dir / "RUN_INFO.yaml"
        try:
            with open(run_info, 'w') as f:
                yaml.dump(manifest, f, default_flow_style=False)
        except Exception as e:
            print(f"Warning: failed to write RUN_INFO.yaml: {e}")
        return run_info

    def run_preflight(self, replicate=1, tag=None):
        """Validate starting structure, common index/reference, and planned windows without running MD."""
        pmf_dir = self.run_dir / "pmf"
        pmf_dir = pmf_dir / (tag if tag else f"replicate_{replicate}")
        pmf_dir.mkdir(parents=True, exist_ok=True)
        self._setup_logging(pmf_dir)
        try:
            self.replicate = int(replicate)
        except Exception:
            self.replicate = replicate

        # Locate starting structure (same logic as run)
        equil_dir = self.run_dir / "equilibration"
        candidates = [
            equil_dir / "npt" / "npt.gro",
            equil_dir / "npt_pr.gro",
            equil_dir / "nvt" / "nvt.gro",
            equil_dir / "em" / "em.gro",
            self.run_dir / "pmf_system" / f"replicate_{replicate}" / "system.gro",
        ]
        start_structure = next((p for p in candidates if p.exists()), None)
        if not start_structure:
            raise FileNotFoundError("Preflight: no starting structure found; run equilibration first.")

        # Optional leaflet index
        leaf = self.run_dir / "index_leaflets.ndx"
        if not leaf.exists():
            try:
                self.generate_leaflet_index(start_structure, leaf)
            except Exception as e:
                self.logger.warning(f"Preflight: could not generate index_leaflets.ndx: {e}")

        # Ensure one common index/reference and validate
        idx_file, ref_group = self.ensure_common_index(start_structure, pmf_dir)
        self._preflight_reference_index(start_structure, pmf_dir)

        # Window plan and overlap geometry
        centers = self.calculate_window_positions()
        adapt_k_cfg = (self.umbrella_config.get('adaptive_k') or {})
        k_base = float(self.umbrella_config.get('force_constant', 900))
        if bool(adapt_k_cfg.get('enabled', False)):
            sigma_target = float(adapt_k_cfg.get('sigma_target_nm', 0.06))
            k_min = float(adapt_k_cfg.get('k_min', 200.0))
            k_max = float(adapt_k_cfg.get('k_max', 5000.0))
            k_eff = float(np.clip(KB_KJ_MOL_K * float(self.temperature_K) / max(1e-9, sigma_target**2), k_min, k_max))
        else:
            k_eff = k_base
        # z-critical for desired overlap (use built-in _ppf)
        try:
            zcrit = float(_ppf(1.0 - float(self.target_overlap)/2.0))
        except Exception:
            zcrit = 1.2815515655446004  # ~ppf(0.9)
        sigma = float(math.sqrt(KB_KJ_MOL_K * float(self.temperature_K) / float(k_eff)))
        dz_max = float(2.0 * sigma * zcrit)

        cs = sorted({float(f"{c:.3f}") for c in centers}, reverse=True)
        densify_pred = []
        for a, b in zip(cs[:-1], cs[1:]):
            if abs(a - b) > 1.05 * dz_max:
                densify_pred.append(float(f"{0.5*(a+b):.3f}"))

        report = {
            "start_structure": str(start_structure.relative_to(self.run_dir)),
            "common_index": str(Path(idx_file).relative_to(self.run_dir)),
            "reference_group": ref_group,
            "temperature_K": float(self.temperature_K),
            "k_base_kj_mol_nm2": float(k_base),
            "k_effective_kj_mol_nm2": float(k_eff),
            "sigma_nm": float(sigma),
            "target_overlap": float(self.target_overlap),
            "dz_max_nm_for_target_overlap": float(dz_max),
            "planned_windows": cs,
            "n_planned": len(cs),
            "auto_densify": bool(self.umbrella_config.get('auto_densify', True)),
            "predicted_midpoints_if_densify": densify_pred,
            "n_predicted_additions": len(densify_pred),
            "max_windows_budget": int((self.qc_config.get('max_windows') or 0) or 10**9),
        }
        with open(pmf_dir / "PREFLIGHT.yaml", "w") as f:
            yaml.dump(report, f, default_flow_style=False)
        print("\n=== Preflight summary ===")
        for k, v in report.items():
            if k in ("planned_windows", "predicted_midpoints_if_densify"):
                continue
            print(f"{k}: {v}")
        print("planned_windows:", report["planned_windows"])
        print("predicted_midpoints_if_densify:", report["predicted_midpoints_if_densify"])
        return report

    def _expected_overlap_gaussian(self, dz: float, k_kj_mol_nm2: float) -> float:
        """Approximate expected overlap of two Gaussians of equal sigma separated by dz.

        sigma = sqrt(kB T / k); overlap ≈ 2 * Phi(-dz/(2 sigma)).
        """
        import math
        sigma = math.sqrt(max(1e-12, KB_KJ_MOL_K * float(self.temperature_K) / float(k_kj_mol_nm2)))
        x = -abs(dz) / (2.0 * sigma)
        # Phi(x) = 0.5[1 + erf(x/sqrt(2))]
        phi = 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
        return max(0.0, min(1.0, 2.0 * phi))

    def _expected_overlap_gaussian_mixed(self, dz: float,
                                         k1_kj_mol_nm2: float,
                                         k2_kj_mol_nm2: float) -> float:
        """Approximate expected overlap for two Gaussians with possibly different sigmas.

        Uses an effective sigma: sigma_eff^2 = (sigma1^2 + sigma2^2) / 2.
        """
        import math
        s1 = math.sqrt(max(1e-12, KB_KJ_MOL_K * float(self.temperature_K) / float(k1_kj_mol_nm2)))
        s2 = math.sqrt(max(1e-12, KB_KJ_MOL_K * float(self.temperature_K) / float(k2_kj_mol_nm2)))
        sigma_eff = math.sqrt(0.5 * (s1 * s1 + s2 * s2))
        x = -abs(dz) / (2.0 * sigma_eff)
        phi = 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
        return max(0.0, min(1.0, 2.0 * phi))

    def validate_window_consistency(self, pmf_dir: Path, windows_data, tolerance=0.1):
        """Validate pairwise consistency between adjacent windows.

        Compares measured vs expected overlaps; writes YAML to pmf_dir/window_consistency.yaml.
        """
        k_global = float(self.umbrella_config.get('force_constant', 900))
        wins = sorted([w for w in windows_data if w], key=lambda x: x['center'], reverse=True)
        report = {"problematic_pairs": [], "passed": True}
        for a, b in zip(wins[:-1], wins[1:]):
            try:
                ol = self.check_window_overlap(a['pullx'], b['pullx'], self.target_overlap)
                dz = abs(float(a['center']) - float(b['center']))
                k1 = float(a.get('k', k_global)); k2 = float(b.get('k', k_global))
                ol_exp = self._expected_overlap_gaussian_mixed(dz, k1, k2)
                if abs(ol - ol_exp) > float(tolerance):
                    report["problematic_pairs"].append({
                        'z1': float(a['center']), 'z2': float(b['center']),
                        'measured': float(ol), 'expected': float(ol_exp)
                    })
            except Exception:
                continue
        report['passed'] = (len(report['problematic_pairs']) == 0)
        out = pmf_dir / "window_consistency.yaml"
        try:
            with open(out, 'w') as f:
                yaml.dump(report, f, default_flow_style=False)
        except Exception:
            pass
        return report
    
    def generate_qc_report(self, qc_results, output_dir):
        report_file = output_dir / "qc_report.yaml"
        # Attach round history if present
        try:
            if getattr(self, 'round_history', None):
                qc_results = dict(qc_results)
                qc_results["rounds"] = list(self.round_history)
            # Attach peptide Rg_xy if available for cross-run comparability
            if getattr(self, '_peptide_rg_xy', None) is not None:
                qc_results["peptide"] = {"rg_xy_nm": float(self._peptide_rg_xy)}
        except Exception:
            pass
        with open(report_file, 'w') as f:
            yaml.dump(qc_results, f, default_flow_style=False)
        def _active(entries):
            return [c for c in (entries or []) if not c.get('superseded', False)]

        overlap_entries = qc_results.get('overlap_check', [])
        ess_entries = _active(qc_results.get('ess_check', []))
        half_entries = _active(qc_results.get('convergence_check', []))

        n_ol_pass = sum(1 for c in overlap_entries if c.get('passed'))
        n_ol_total = len(overlap_entries)
        n_ess_pass = sum(1 for c in ess_entries if c.get('passed'))
        n_ess_total = len(ess_entries)
        n_half_pass = sum(1 for c in half_entries if c.get('passed'))
        n_half_total = len(half_entries)
        print("\n=== QC Report ===")
        print(f"Overlap checks: {n_ol_pass}/{n_ol_total} passed (region-aware thresholds)")
        print(f"ESS checks:     {n_ess_pass}/{n_ess_total} passed (region-aware thresholds)")
        print(f"Half-time:      {n_half_pass}/{n_half_total} passed (ΔE + Z, region-aware; drift_tol {self.qc_config.get('drift_tolerance_nm_per_ns',0.1)} nm/ns)")
        for check in qc_results.get('overlap_check', []):
            if not check['passed']:
                print(f"  ⚠ Low overlap ({check['overlap']:.3f}) between z={check['windows'][0]:.2f} and z={check['windows'][1]:.2f}")
        return report_file

    def _filter_windows_by_qc(self, windows_data, pmf_dir: Path):
        """Filter out bad replicates per center based on QC criteria.

        Criteria come from pmf.qc.filter:
          - drop_if_half_fail (bool)
          - drop_if_drift_fail (bool)
          - min_ess_keep (int)
          - min_keep_per_center (int)
          - keep_best_if_all_filtered (bool)
          - prefer_by: 'ess' | 'half_then_ess'

        Returns (selected_windows, selection_report)
        selection_report is a dict with center-wise details and reasons.
        """
        fcfg = (self.qc_config.get('filter') or {})
        if not bool(fcfg.get('enabled', True)):
            return [w for w in windows_data if w], {}
        drop_half = bool(fcfg.get('drop_if_half_fail', True))
        drop_drift = bool(fcfg.get('drop_if_drift_fail', False))
        min_ess = int(fcfg.get('min_ess_keep', 20))
        min_keep = int(fcfg.get('min_keep_per_center', 1))
        keep_best_anyway = bool(fcfg.get('keep_best_if_all_filtered', True))
        max_keep = fcfg.get('max_keep_per_center', 1)
        try:
            max_keep = int(max_keep)
        except Exception:
            max_keep = 1
        if max_keep < 1:
            max_keep = 0
        prefer_by = str(fcfg.get('prefer_by', 'ess')).lower()
        respect_gates = bool(fcfg.get('respect_gates', True))

        # Group windows by center
        centers = {}
        for w in [w for w in windows_data if w]:
            centers.setdefault(float(w['center']), []).append(w)

        k_global = float(self.umbrella_config.get('force_constant', 900))
        selection = {}
        selected_windows = []
        for c, wins in sorted(centers.items(), key=lambda kv: kv[0], reverse=True):
            items = []
            gating_candidates = []
            for w in wins:
                reg = self._region_cfg(w['center'])
                kw = float(w.get('k', k_global))
                qc = self._get_window_qc(w['pullx'], kw, w['center'], reg.get('half_energy_tol_kj', 2.0))
                # Region-aware half decision (ΔE + Z)
                half_pass_reg = (float(qc.get('half_energy_diff_kj', float('inf'))) <= float(reg["half_energy_tol_kj"])) and \
                                (float(qc.get('half_z_stat', 0.0)) <= float(reg["half_z_tol_sigma"]))
                reason = []
                drop = False
                ess_val = int(qc.get('ess', 0))
                ess_gate = int(reg.get('min_ess_frames', min_ess))
                if drop_half and not half_pass_reg:
                    drop = True; reason.append('half_fail')
                if drop_drift and not qc.get('drift_pass', True):
                    drop = True; reason.append('drift_fail')
                if ess_val < min_ess:
                    drop = True; reason.append(f"ess<{min_ess}")
                passes_gate = (ess_val >= ess_gate) and (not drop_half or half_pass_reg)
                if respect_gates and passes_gate and drop:
                    if any(r.startswith('ess<') for r in reason):
                        reason = [r for r in reason if not r.startswith('ess<')]
                    if not (drop_half and not half_pass_reg) and not (drop_drift and not qc.get('drift_pass', True)):
                        drop = False
                if passes_gate:
                    gating_candidates.append({'window': w, 'qc': qc, 'half_pass': half_pass_reg})
                items.append({
                    'center': float(c),
                    'pullx': w['pullx'],
                    'pullf': w['pullf'],
                    'tpr': w.get('tpr'),
                    'seed': w.get('seed'),
                    'qc': {
                        'ess': int(qc.get('ess', 0)),
                        'half_pass': bool(half_pass_reg),
                        'drift_pass': bool(qc.get('drift_pass', True)),
                        'half_energy_diff_kj': float(qc.get('half_energy_diff_kj', 0.0)),
                        'drift_nm_per_ns': float(qc.get('drift_nm_per_ns', 0.0)) if qc.get('drift_nm_per_ns') is not None else None,
                    },
                    'drop': bool(drop),
                    'reason': reason,
                })
            # Select survivors
            survivors = [it for it in items if not it['drop']]
            if len(survivors) < min_keep:
                if keep_best_anyway and items:
                    # sort by preference
                    if prefer_by.startswith('half'):
                        items_sorted = sorted(items, key=lambda it: ((not it['qc']['half_pass']), -it['qc']['ess']))
                    else:
                        items_sorted = sorted(items, key=lambda it: (-it['qc']['ess']))
                    survivors = items_sorted[:max(1, min_keep)]
                elif respect_gates and gating_candidates:
                    if prefer_by.startswith('half'):
                        gating_sorted = sorted(gating_candidates, key=lambda it: ((not it['half_pass']), -it['qc'].get('ess', 0)))
                    else:
                        gating_sorted = sorted(gating_candidates, key=lambda it: (-it['qc'].get('ess', 0)))
                    survivors = []
                    for cand in gating_sorted[:max(1, min_keep)]:
                        target_pullx = cand['window']['pullx']
                        match = next((it for it in items if it['pullx'] == target_pullx), None)
                        if match:
                            match['drop'] = False
                            if f"ess<{min_ess}" in match['reason']:
                                match['reason'] = [r for r in match['reason'] if r != f"ess<{min_ess}"]
                            survivors.append(match)
                else:
                    survivors = []
            # Hard cap on survivors per center if requested
            if max_keep and len(survivors) > max_keep:
                if prefer_by.startswith('half'):
                    survivors = sorted(survivors, key=lambda it: ((not it['qc']['half_pass']), -it['qc']['ess']))[:max_keep]
                else:
                    survivors = sorted(survivors, key=lambda it: (-it['qc']['ess']))[:max_keep]
            selection[float(c)] = {
                'replicates': items,
                'selected': [{'pullx': it['pullx'], 'ess': it['qc']['ess'], 'half_pass': it['qc']['half_pass']} for it in survivors],
                'n_selected': len(survivors),
                'n_total': len(items),
                'ess_sum_selected': int(sum(it['qc']['ess'] for it in survivors))
            }
            # Append survivors in original window dict form
            px_to_win = {w['pullx']: w for w in wins}
            for it in survivors:
                sel = px_to_win.get(it['pullx'])
                if sel:
                    selected_windows.append(sel)

        # Write selection report
        try:
            with open(pmf_dir / 'qc_filter.yaml', 'w') as f:
                yaml.dump(selection, f, default_flow_style=False)
        except Exception:
            pass
        return selected_windows, selection

    def _annotate_superseded_qc(self, qc_results: dict, selection: dict) -> dict:
        """Mark QC entries belonging to filtered-out replicates as superseded.

        Adds a boolean flag `superseded` to convergence/ESS entries when the
        corresponding replicate was dropped by the QC filter so that reports
        can hide or down-weight those failures.
        """
        if not selection:
            return qc_results
        try:
            selected_pullx = set()
            all_pullx = set()
            for block in selection.values():
                for sel in block.get('selected', []) or []:
                    if sel.get('pullx'):
                        selected_pullx.add(str(sel['pullx']))
                for rep in block.get('replicates', []) or []:
                    if rep.get('pullx'):
                        all_pullx.add(str(rep['pullx']))
            superseded = all_pullx - selected_pullx
            if not superseded:
                # Still annotate kept entries for completeness
                relevant = selected_pullx
            else:
                relevant = selected_pullx
            for section in ("convergence_check", "ess_check"):
                for item in qc_results.get(section, []) or []:
                    pullx = item.get('pullx')
                    if not pullx:
                        continue
                    pullx_str = str(pullx)
                    if pullx_str in superseded:
                        item['superseded'] = True
                    elif pullx_str in relevant:
                        item['superseded'] = False
        except Exception:
            pass
        return qc_results
    
    def _scan_existing_windows(self, pmf_dir: Path):
        """Return list of window dicts for existing window outputs, including replicates.

        Scans windows/z_*/ and windows/z_*/rep_*/ for pullx/pullf/tpr entries.
        """
        windows_root = pmf_dir / "windows"
        if not windows_root.exists():
            return []
        wins = []

        def _append_entry(dir_path: Path, center: float):
            pullx = dir_path / "pullx.xvg"
            pullf = dir_path / "pullf.xvg"
            tpr = dir_path / "umbrella.tpr"
            if pullx.exists():
                # Try to parse force constant k from umbrella.mdp
                k_val = None
                try:
                    mdp_path = dir_path / "umbrella.mdp"
                    if mdp_path.exists():
                        with open(mdp_path, 'r') as mf:
                            for ln in mf:
                                if 'pull-coord1-k' in ln:
                                    parts = ln.strip().split()
                                    for tok in reversed(parts):
                                        try:
                                            k_val = float(tok)
                                            break
                                        except Exception:
                                            continue
                                    break
                except Exception:
                    k_val = None
                entry = {
                    "center": float(center),
                    "pullf": str(pullf),
                    "pullx": str(pullx),
                    "tpr": str(tpr)
                }
                if k_val is not None:
                    entry["k"] = float(k_val)
                wins.append(entry)

        for zdir in sorted(windows_root.glob("z_*")):
            try:
                m = re.match(r"z_([+-]?(?:\d+(?:\.\d+)?))$", zdir.name)
                if not m:
                    continue
                center = float(m.group(1))
            except Exception:
                continue

            # Base window
            _append_entry(zdir, center)
            # Replicates
            for repdir in sorted(zdir.glob("rep_*")):
                _append_entry(repdir, center)

        # No dedup by center here; we keep all replicates
        # Optionally dedup exact duplicates by (center, pullx)
        seen = set(); kept = []
        for w in wins:
            key = (float(w['center']), w.get('pullx'))
            if key in seen:
                continue
            seen.add(key); kept.append(w)
        return kept

    def run_pmf_calculation(self, replicate=1, tag=None, resume=False, qc_only=False):
        """Main PMF calculation workflow with all improvements"""

        # Setup directories
        pmf_dir = self.run_dir / "pmf"
        if tag:
            pmf_dir = pmf_dir / tag
        else:
            pmf_dir = pmf_dir / f"replicate_{replicate}"

        pmf_dir.mkdir(parents=True, exist_ok=True)
        # Expose replicate for topology selection and provenance
        try:
            self.replicate = int(replicate)
        except Exception:
            self.replicate = replicate
        self._setup_logging(pmf_dir)
        self.logger.info(f"PMF run started for {self.run_dir} -> {pmf_dir}")
        # Reset QC cache for this run to avoid stale metrics across runs
        self._window_qc_cache = {}

        # Get starting structure - try different locations
        equil_dir = self.run_dir / "equilibration"

        # Try different equilibration outputs
        possible_structures = [
            equil_dir / "npt" / "npt.gro",     # Standard NPT output
            equil_dir / "npt_pr.gro",          # Old format
            equil_dir / "nvt" / "nvt.gro",     # Fallback to NVT
            equil_dir / "em" / "em.gro",       # Last resort: EM
            self.run_dir / "pmf_system" / f"replicate_{replicate}" / "system.gro",  # PMF system
        ]

        start_structure = None
        for structure in possible_structures:
            if structure.exists():
                start_structure = structure
                print(f"Using starting structure: {structure.relative_to(self.run_dir)}")
                break

        if not start_structure:
            self.logger.error("No starting structure found. Checked:")
            for structure in possible_structures:
                self.logger.error(f"  - {structure.relative_to(self.run_dir)}: {'✓' if structure.exists() else '✗'}")
            raise FileNotFoundError("No starting structure found")

        # Estimate peptide footprint proxy (Rg_xy) from starting structure (nm)
        peptide_rg_xy = None
        try:
            atoms0 = self.read_gro_atoms(start_structure)
            pep_idx0 = self.find_peptide_indices(atoms0)
            if pep_idx0:
                xy = np.array([[x, y] for i, (_r,_a,x,y,_z) in enumerate(atoms0) if (i+1) in pep_idx0], dtype=float)
                com = xy.mean(axis=0)
                d2 = ((xy - com)**2).sum(axis=1)
                peptide_rg_xy = float(np.sqrt(d2.mean())) if d2.size else None
        except Exception:
            peptide_rg_xy = None

        # Expose Rg_xy to other writers
        try:
            self._peptide_rg_xy = float(peptide_rg_xy) if peptide_rg_xy is not None else None
        except Exception:
            self._peptide_rg_xy = None

        # Size-aware reference & QC tuning
        try:
            rgxy = float(self._peptide_rg_xy) if getattr(self, '_peptide_rg_xy', None) is not None else 0.0
        except Exception:
            rgxy = 0.0
        try:
            a = float((self.umbrella_config.get('patch_radius_scale') or {}).get('a_nm', 0.5))
            b = float((self.umbrella_config.get('patch_radius_scale') or {}).get('b_x', 1.2))
            rmin = float(self.umbrella_config.get('patch_radius_min_nm', 2.0))
            if bool(self.umbrella_config.get('patch_radius_auto', True)):
                self.autopatch_radius_nm = max(rmin, a + b * rgxy)
            else:
                self.autopatch_radius_nm = float(self.umbrella_config.get('patch_radius', 2.5))
            self.logger.info(f"Auto patch radius set to {self.autopatch_radius_nm:.2f} nm (Rg_xy={rgxy:.2f} nm)")
        except Exception:
            self.autopatch_radius_nm = float(self.umbrella_config.get('patch_radius', 2.5))
        # Gentle QC bump in interface region for large peptides
        try:
            tun = (self.qc_config.get('size_tuning') or {})
            rg_thr = float(tun.get('rgxy_threshold_nm', 1.6))
            if rgxy >= rg_thr:
                bump = float(tun.get('interface_overlap_bump', 0.02))
                self.qc_config.setdefault('region_thresholds', {})
                self.qc_config['region_thresholds'].setdefault('interface', {})
                cur = float(self.qc_config['region_thresholds']['interface'].get('min_neighbor_overlap', 0.20))
                self.qc_config['region_thresholds']['interface']['min_neighbor_overlap'] = max(cur, cur + bump)
                self.logger.info(
                    f"Bumped interface min_neighbor_overlap to "
                    f"{self.qc_config['region_thresholds']['interface']['min_neighbor_overlap']:.2f} due to size"
                )
        except Exception:
            pass

        # Calculate baseline window positions
        window_centers = self.calculate_window_positions()
        self.logger.info(f"Planning {len(window_centers)} umbrella windows")
        self.logger.info(f"Starting structure: {start_structure.relative_to(self.run_dir)}")

        # Create a leaflet index once if missing (improves patch/midplane reference)
        leaflet_ndx = self.run_dir / "index_leaflets.ndx"
        try:
            if not leaflet_ndx.exists():
                self.generate_leaflet_index(start_structure, leaflet_ndx)
        except Exception as e:
            self.logger.warning(f"Could not generate index_leaflets.ndx automatically: {e}")

        # Create a single, consistent index for all windows (avoids reference bias)
        if self.consistent_reference:
            self.ensure_common_index(start_structure, pmf_dir)
            # Preflight the common index and reference group to fail fast
            try:
                self._preflight_reference_index(start_structure, pmf_dir)
            except Exception as e:
                self.logger.error(f"Reference preflight failed: {e}")
                raise

        # Fast paths
        if qc_only:
            self.logger.info("qc-only mode: scanning existing windows and generating reports")
            windows_data = self._scan_existing_windows(pmf_dir)
            if not windows_data:
                self.logger.error("No existing windows found to analyze.")
                raise SystemExit(2)
            qc_results = self.run_qc_checks(windows_data)
            windows_selected, selection = self._filter_windows_by_qc(windows_data, pmf_dir)
            qc_results = self._annotate_superseded_qc(qc_results, selection)
            if windows_selected:
                windows_data = windows_selected
            self.generate_qc_report(qc_results, pmf_dir)
            run_info = self.write_run_manifest(pmf_dir, [w['center'] for w in windows_data if w], qc_results)
            metadata = {
                "replicate": replicate,
                "n_windows": len(windows_data),
                "window_centers": [w['center'] for w in windows_data if w],
                "force_constant": self.umbrella_config.get('force_constant', 900),
                "production_time": self.umbrella_config.get('production_ns', 20.0),
                "reference_mode": self.umbrella_config.get('ref_mode', 'patch'),
                "consistent_reference": bool(self.consistent_reference),
                "reference_group": self.common_ref_group if self.common_ref_group else self.umbrella_config.get('ref_mode', 'patch'),
                "index_file": str(Path(self.common_index_file).relative_to(self.run_dir)) if self.common_index_file else None,
                "patch_radius": float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.0))),
                "windows": windows_data,
                "qc_results": qc_results,
                "peptide": ({"rg_xy_nm": float(peptide_rg_xy)} if peptide_rg_xy is not None else None)
            }
            if selection:
                metadata['qc_selection'] = selection
            with open(pmf_dir / "pmf_metadata.yaml", 'w') as f:
                yaml.dump(metadata, f, default_flow_style=False)
            self.logger.info(f"QC + metadata written. Manifest: {run_info.name}")
            return metadata

        # Run windows
        windows_data = []

        # Resume: consider existing windows first
        if resume:
            existing = self._scan_existing_windows(pmf_dir)
            if existing:
                self.logger.info(f"Resume mode: found {len(existing)} existing windows; will skip them")
                windows_data.extend(existing)
            else:
                self.logger.info("Resume mode: no existing windows found")
        for i, center in enumerate(window_centers):
            self.logger.info(f"[Window {i+1}/{len(window_centers)}] z = {center:+.3f} nm")

            window_dir = pmf_dir / "windows" / f"z_{center:+.3f}"

            # Check if window already exists
            if (window_dir / "pullf.xvg").exists():
                self.logger.info("  ✓ Already completed, skipping")
                windows_data.append({
                    "center": center,
                    "pullf": str(window_dir / "pullf.xvg"),
                    "pullx": str(window_dir / "pullx.xvg"),
                    "tpr": str(window_dir / "umbrella.tpr")
                })
                continue

            # Show progress
            self.logger.info("  Setting up window...")
            sys.stdout.flush()

            # Run window
            result = self.run_umbrella_window(
                center, window_dir, start_structure, replicate
            )
            windows_data.append(result)

            # Show overall progress
            completed = sum(1 for w in windows_data if w is not None)
            self.logger.info(f"  Overall progress: {completed}/{len(window_centers)} windows completed")
            sys.stdout.flush()

            # Update start structure for next window (optional)
            if result and self.umbrella_config.get('use_previous_frame', False):
                start_structure = window_dir / "umbrella.gro"

        # Optional: deduplicate exact duplicates by (center, pullx) while keeping replicates
        seen = set(); kept = []
        for w in windows_data:
            if not w:
                continue
            key = (float(w['center']), w.get('pullx'))
            if key in seen:
                continue
            seen.add(key)
            kept.append(w)
        windows_data = kept

        # QC loop with auto-densify/extend until gates pass (ESS + overlap)
        # If qc.max_qc_rounds <= 0, run indefinitely until gates pass or no further changes are possible.
        max_rounds = int(self.qc_config.get('max_qc_rounds', 5))
        # If qc.until_gates_pass is true (default), ignore max_rounds and run until gates are satisfied
        until_pass = bool(self.qc_config.get('until_gates_pass', True))
        unlimited = until_pass or (max_rounds <= 0)
        all_centers = [w['center'] for w in windows_data if w]
        round_idx = 0
        while True:
            # Initialize per-round log
            round_log = {
                "round": round_idx + 1,
                "added": [],
                "extended": [],
                "replicated": [],
                "reasons": {"overlap_pairs": [], "gaps": []},
            }
            qc_results = self.run_qc_checks(windows_data)
            overlap_ok, ess_ok = self._qc_gates_passed(qc_results)
            if overlap_ok and ess_ok:
                self.logger.info(f"QC gates satisfied after {round_idx} round(s).")
                # Still persist a round entry showing no actions
                self.round_history.append(round_log)
                break

            to_add = []
            auto_densify = bool(self.umbrella_config.get('auto_densify', True))
            if auto_densify:
                # Add midpoint windows for low-overlap neighbors
                for check in qc_results.get('overlap_check', []):
                    if not check.get('passed', False):
                        mid = _round3(0.5*(check['windows'][0] + check['windows'][1]))
                        if mid not in all_centers:
                            to_add.append(mid)
                            try:
                                round_log["reasons"]["overlap_pairs"].append({
                                    "z1": float(check['windows'][0]),
                                    "z2": float(check['windows'][1]),
                                    "overlap": float(check.get('overlap', 0.0)),
                                    "min_required": float(check.get('min_required', 0.0)),
                                    "mid": float(mid),
                                })
                            except Exception:
                                pass
                # Add midpoint windows for large geometric gaps (safety)
                k = float(self.umbrella_config.get('force_constant', 900))
                dz_max = _dz_max_for_overlap(float(self.temperature_K), k, float(self.target_overlap))
                all_centers_sorted = _dedup_sorted_centers([w['center'] for w in windows_data if w])
                for a,b in zip(all_centers_sorted[:-1], all_centers_sorted[1:]):
                    if abs(a-b) > dz_max*1.05:
                        mid = _round3(0.5*(a+b))
                        if mid not in all_centers:
                            to_add.append(mid)
                            try:
                                round_log["reasons"]["gaps"].append({
                                    "z1": float(a), "z2": float(b), "mid": float(mid), "dz": float(abs(a-b))
                                })
                            except Exception:
                                pass
            to_add = _dedup_sorted_centers(to_add)
            # Respect max windows budget if configured
            try:
                max_w = int(self.qc_config.get('max_windows', 10**9))
                current_n = len(_dedup_sorted_centers([w['center'] for w in windows_data if w]))
                remaining = max(0, max_w - current_n)
                if len(to_add) > remaining:
                    to_add = to_add[:remaining]
            except Exception:
                pass
            added = 0
            if to_add:
                self.logger.info(f"QC round {round_idx+1}: adding {len(to_add)} windows to improve coverage/overlap")
                for new_center in to_add:
                    window_dir = pmf_dir / "windows" / f"z_{new_center:+.3f}"
                    res = self.run_umbrella_window(new_center, window_dir, start_structure, replicate)
                    if res:
                        windows_data.append(res)
                        all_centers.append(new_center)
                        added += 1
                        round_log["added"].append({"center": float(new_center), "reason": "overlap/gap"})

            # After densification in this round, try extending windows that failed ESS/half-time
            qc_results = self.run_qc_checks(windows_data)
            windows_data, extended, extended_entries = self._extend_failed_windows(
                pmf_dir, windows_data, qc_results, start_structure, replicate
            )
            if extended_entries:
                round_log["extended"].extend(extended_entries)
            # Round summary logging
            try:
                self.logger.info(f"QC round {round_idx+1} summary: added={added}, extended={len(extended_entries)}")
            except Exception:
                pass

            replicated_now = 0
            fallback_cfg = (self.qc_config.get('fallback') or {})
            if bool(fallback_cfg.get('auto_after_extend', True)):
                include_half = bool(fallback_cfg.get('include_half_convergence', True))
                qc_results = self.run_qc_checks(windows_data)
                early_targets = self._collect_problem_centers(qc_results, include_half=include_half)
                if early_targets:
                    self.logger.info(f"QC round {round_idx+1}: early fallback for {len(early_targets)} window(s)")
                    windows_data, replicated_now, repl_entries_now = self._replicate_low_ess_windows(
                        pmf_dir, windows_data, qc_results, start_structure, replicate, targets_override=early_targets
                    )
                    if repl_entries_now:
                        round_log["replicated"].extend(repl_entries_now)
                        all_centers.extend([float(entry['center']) for entry in repl_entries_now])
                    qc_results = self.run_qc_checks(windows_data)
                    overlap_ok, ess_ok = self._qc_gates_passed(qc_results)
                else:
                    overlap_ok, ess_ok = self._qc_gates_passed(qc_results)
            else:
                overlap_ok, ess_ok = self._qc_gates_passed(qc_results)

            # Decide whether to continue
            round_idx += 1
            if not unlimited and round_idx >= max_rounds and not (overlap_ok and ess_ok):
                self.logger.warning("Reached max_qc_rounds before all QC gates passed.")
                self.round_history.append(round_log)
                break
            if added == 0 and extended == 0 and replicated_now == 0:
                # No more actions possible; avoid infinite loop
                self.logger.warning("QC gates not fully satisfied, but no windows were added or extended in this round. Stopping.")
                self.round_history.append(round_log)
                break
            # Persist per-round log
            self.round_history.append(round_log)
        # Final QC report (pre-fallback)
        qc_results = self.run_qc_checks(windows_data)
        overlap_ok, ess_ok = self._qc_gates_passed(qc_results)
        if not (overlap_ok and ess_ok):
            self.logger.warning("QC gates not fully satisfied after QC rounds. Engaging fallback replicates…")
            windows_data, added, replicated_entries = self._replicate_low_ess_windows(
                pmf_dir, windows_data, qc_results, start_structure, replicate
            )
            # Fallback logged as a dedicated round
            if replicated_entries:
                self.round_history.append({
                    "round": round_idx + 1,  # next sequential
                    "added": [],
                    "extended": [],
                    "replicated": replicated_entries,
                    "reasons": {"overlap_pairs": [], "gaps": []},
                    "note": "fallback",
                })
                try:
                    self.logger.info(f"QC fallback summary: replicated={len(replicated_entries)}")
                except Exception:
                    pass
            # Re-check QC after fallback
            qc_results = self.run_qc_checks(windows_data)
        # Select good replicates per center and update metadata payload
        windows_selected, selection = self._filter_windows_by_qc(windows_data, pmf_dir)
        qc_results = self._annotate_superseded_qc(qc_results, selection)
        overlap_ok_final, ess_ok_final = self._qc_gates_passed(qc_results)
        strict_info = qc_results.setdefault('strict_mode', {})
        strict_info['enabled'] = bool(self.strict_qc)
        strict_info['ess_pass'] = bool(ess_ok_final)
        strict_info['overlap_pass'] = bool(overlap_ok_final)
        self.generate_qc_report(qc_results, pmf_dir)
        if windows_selected:
            windows_data = windows_selected
        # Consistency validation
        self.validate_window_consistency(pmf_dir, windows_data)

        # Write run manifest for provenance
        run_info = self.write_run_manifest(pmf_dir, [w['center'] for w in windows_data if w], qc_results)

        # Generate metadata for analysis
        metadata = {
            "replicate": replicate,
            "n_windows": len(windows_data),
            "window_centers": [w['center'] for w in windows_data if w],
            "force_constant": self.umbrella_config.get('force_constant', 900),
            "production_time": self.umbrella_config.get('production_ns', 60.0),
            "reference_mode": self.umbrella_config.get('ref_mode', 'patch'),
            "consistent_reference": bool(self.consistent_reference),
            "reference_group": self.common_ref_group if self.common_ref_group else self.umbrella_config.get('ref_mode', 'patch'),
            "index_file": str(Path(self.common_index_file).relative_to(self.run_dir)) if self.common_index_file else None,
                "patch_radius": float(getattr(self, 'autopatch_radius_nm', self.umbrella_config.get('patch_radius', 2.0))),
            "windows": windows_data,
            "qc_results": qc_results,
            "peptide": ({"rg_xy_nm": float(peptide_rg_xy)} if peptide_rg_xy is not None else None)
        }
        # Persist selection summary for traceability
        if selection:
            metadata['qc_selection'] = selection

        metadata_file = pmf_dir / "pmf_metadata.yaml"
        with open(metadata_file, 'w') as f:
            yaml.dump(metadata, f, default_flow_style=False)

        # Final summary with more details + completion gating
        successful_windows = len([w for w in windows_data if w])
        # Re-evaluate gates for final status
        overlap_ok, ess_ok = overlap_ok_final, ess_ok_final
        gate_cfg = (self.qc_config.get('completion_gate') or {})
        try:
            min_frac = float(gate_cfg.get('min_fraction', 0.6))
        except Exception:
            min_frac = 0.6
        completion_ok = (successful_windows / max(1, len(window_centers)) >= min_frac) and overlap_ok and ess_ok
        strict_info['failed'] = bool(self.strict_qc and not (overlap_ok and ess_ok))
        if strict_info['failed']:
            self.strict_failure = True
            self.logger.error("Strict QC mode enabled: ESS/overlap gates failed; aborting downstream analysis.")
        if completion_ok:
            self.logger.info("✅ PMF CALCULATION COMPLETED")
        else:
            self.logger.warning("❌ PMF CALCULATION INCOMPLETE (coverage/QC gates not satisfied)")
        self.logger.info(f"  Peptide: {self.run_dir.name.split('_')[0]}")
        self.logger.info(f"  Replicate: {replicate}")
        self.logger.info(f"  Windows completed: {successful_windows}/{len(window_centers)}")
        try:
            rel_out = pmf_dir.relative_to(self.run_dir)
        except Exception:
            rel_out = pmf_dir
        self.logger.info(f"  Output directory: {rel_out}")
        self.logger.info(f"  Metadata: pmf_metadata.yaml")
        self.logger.info(f"  QC Report: qc_report.yaml")
        self.logger.info(f"  Run manifest: {run_info.name}")
        if successful_windows < len(window_centers):
            self.logger.warning(f"⚠ {len(window_centers) - successful_windows} windows failed")
            self.logger.warning("  Check the output for error messages")
        if completion_ok:
            self.logger.info("Next step: Run MBAR analysis")
            self.logger.info(f"Command: python scripts/analysis/pmf_mbar_analysis.py {pmf_dir}")
        else:
            self.logger.warning("Skip MBAR: insufficient coverage/QC. Inspect qc_report.yaml and logs.")

        return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced PMF runner with local patch reference"
    )
    parser.add_argument("run_dir", help="Simulation run directory")
    parser.add_argument("--replicate", type=int, default=1, help="Replicate number")
    parser.add_argument("--tag", help="Optional tag for output directory")
    parser.add_argument("--resume", action="store_true", help="Resume an incomplete run (skip existing windows)")
    parser.add_argument("--qc-only", action="store_true", help="Only generate qc_report.yaml + pmf_metadata.yaml from existing windows")
    parser.add_argument("--preflight-only", action="store_true", help="Validate reference and plan windows; no MD")
    parser.add_argument("--strict-qc", action="store_true", help="Abort downstream analysis when ESS/overlap gates fail")
    
    args = parser.parse_args()
    
    # Initialize runner
    runner = PMFRunner(args.run_dir)
    if args.preflight_only:
        runner.run_preflight(replicate=args.replicate, tag=args.tag)
        return
    if args.strict_qc:
        runner.strict_qc = True
    # Run PMF calculation
    runner.run_pmf_calculation(replicate=args.replicate, tag=args.tag, resume=args.resume, qc_only=args.qc_only)
    if runner.strict_failure:
        sys.exit(2)

if __name__ == "__main__":
    main()
