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

# --- math/helpers for QC and spacing ---
def _dz_max_for_overlap(T_K: float, k_kj_mol_nm2: float, overlap: float) -> float:
    """Given target histogram overlap O and harmonic k, return max allowed Δz."""
    sigma = math.sqrt(KB_KJ_MOL_K * T_K / k_kj_mol_nm2)
    # 1D Gaussian: choose z so that two Gaussians separated by Δz have min-overlap ~O
    # This is approximated by Δz_max ≈ 2 * sigma * z_alpha with alpha = 1 - O/2
    # using inverse normal CDF; numeric approx here via scipy-free Acklam like in old script
    def _norm_ppf(p: float) -> float:
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
        if p < plow:
            q = math.sqrt(-2.0 * math.log(p))
            return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                   ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        if p > phigh:
            q = math.sqrt(-2.0 * math.log(1.0 - p))
            return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                     ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        q = p - 0.5; r = q*q
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0)
    z = _norm_ppf(1.0 - overlap / 2.0)
    return 2.0 * sigma * z

def _round3(x: float) -> float:
    return float(f"{x:.3f}")

def _dedup_sorted_centers(centers):
    """Round to 3 decimals, drop duplicates, return sorted(desc)."""
    s = sorted({_round3(c) for c in centers}, reverse=True)
    return s

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
        # QC: optional stride for ESS to mitigate oversampling bias
        self.qc_ess_stride = int(self.qc_config.get('ess_stride', 1))
        # Ensure same reference group across all windows to avoid analysis bias
        self.consistent_reference = bool(self.umbrella_config.get('consistent_reference', True))
        self.common_index_file = None
        self.common_ref_group = None
        self.logger = logging.getLogger("pmf")
        self.logger.setLevel(logging.INFO)
        # Handlers are attached lazily in run_pmf_calculation when pmf_dir is known

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

    def create_midplane_patch_reference_group(self, gro_atoms, peptide_indices, patch_radius=2.5, min_per_leaflet=10):
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
            return [(i, x, y) for (i, x, y) in lst if (x - pep_xy[0])**2 + (y - pep_xy[1])**2 <= R*R]
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
            sel = sorted(sel, key=lambda t: (t[1]-pep_xy[0])**2 + (t[2]-pep_xy[1])**2)
            return [i for (i, _x, _y) in sel[:m]]
        o_ids = take_nearest(o_sel, m)
        i_ids = take_nearest(i_sel, m)
        return o_ids + i_ids

    def create_patch_reference_group(self, gro_atoms, peptide_indices, patch_radius=2.0):
        """
        Create local patch reference from phosphates near peptide using OuterPO4 indices
        if available; fallback to all PO4/P atoms.
        """
        # peptide COM (XYZ), project to XY
        pep_xyz = np.array([ [x,y,z] for idx, (_,_,x,y,z) in enumerate(gro_atoms) if (idx+1) in peptide_indices ])
        if pep_xyz.size == 0:
            raise ValueError("No peptide atoms found")
        pep_xy = pep_xyz[:, :2].mean(axis=0)

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
            patch = [idx for idx,x,y in candidates if ( (x - pep_xy[0])**2 + (y - pep_xy[1])**2 ) <= R*R ]
            if len(patch) >= 10 or R >= patch_radius + 2.0:
                break
            R += 0.5
        if len(patch) < 10:
            print(f"Warning: only {len(patch)} phosphates in patch (R={R:.1f} nm)")
        return patch
    
    def create_dynamic_index_file(self, gro_path, output_ndx, window_center=None):
        """Create index file with (by default) a local patch reference.

        For PMF, we prefer to create this once from a single reference structure
        and reuse it for all windows (see ensure_common_index) so that the
        biased coordinate is defined against the same atoms across windows.
        """
        gro_atoms = self.read_gro_atoms(gro_path)

        # Find peptide
        peptide_indices = self.find_peptide_indices(gro_atoms)
        if not peptide_indices:
            raise ValueError("No peptide found in structure")

        # Get reference mode
        ref_mode = self.umbrella_config.get('ref_mode', 'patch')

        if ref_mode == 'midplane_local':
            patch_radius = self.umbrella_config.get('patch_radius', 2.5)
            reference_indices = self.create_midplane_patch_reference_group(
                gro_atoms, peptide_indices, patch_radius
            )
            ref_group_name = "LocalMidplane"
        elif ref_mode == 'hybrid':
            # Try midplane_local first, validate; fall back to local patch, then global PO4
            ref_group_name = "LocalMidplane"
            patch_radius = self.umbrella_config.get('patch_radius', 2.5)
            try:
                reference_indices = self.create_midplane_patch_reference_group(
                    gro_atoms, peptide_indices, patch_radius
                )
            except Exception:
                reference_indices = []
            if not reference_indices or not self._validate_reference_group(gro_atoms, reference_indices, 0.95):
                ref_group_name = "LocalPatch"
                patch_radius = self.umbrella_config.get('patch_radius', 2.0)
                try:
                    reference_indices = self.create_patch_reference_group(
                        gro_atoms, peptide_indices, patch_radius
                    )
                except Exception:
                    reference_indices = []
            if not reference_indices or not self._validate_reference_group(gro_atoms, reference_indices, 0.95):
                ref_group_name = "UpperPO4"
                reference_indices = [i+1 for i, (_r, a, *_xyz) in enumerate(gro_atoms) if a in ("PO4","P")]
        elif ref_mode == 'patch':
            # LOCAL PATCH REFERENCE (key improvement)
            patch_radius = self.umbrella_config.get('patch_radius', 2.0)
            reference_indices = self.create_patch_reference_group(
                gro_atoms, peptide_indices, patch_radius
            )
            ref_group_name = "LocalPatch"
        else:
            reference_indices = [i+1 for i, (_r, a, *_xyz) in enumerate(gro_atoms) if a in ("PO4","P")]
            ref_group_name = "UpperPO4"

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
                for candidate in ("LocalMidplane", "LocalPatch", "UpperPO4"):
                    if candidate in groups:
                        self.common_ref_group = candidate
                        break
                else:
                    # Fallback if none found
                    self.common_ref_group = "LocalPatch"
                # Compute a central pbcatom for the reference group
                gro_atoms = self.read_gro_atoms(gro_path)
                self.common_pbcatom1_index = self._compute_pbcatom(idx_path, self.common_ref_group, gro_atoms)
            except Exception:
                self.common_ref_group = 'LocalPatch'
            return self.common_index_file, self.common_ref_group
        # Otherwise, create new according to config
        idx_file, ref_group = self.create_dynamic_index_file(gro_path, idx_path)
        self.common_index_file = str(idx_file)
        self.common_ref_group = ref_group
        # Compute central pbcatom
        gro_atoms = self.read_gro_atoms(gro_path)
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

    def _validate_reference_group(self, gro_atoms, ref_indices, threshold: float = 0.95) -> bool:
        """Validate reference group spatial consistency (simple outlier test)."""
        try:
            if not ref_indices or len(ref_indices) < 5:
                return False
            coords = np.array([gro_atoms[i-1][2:5] for i in ref_indices], dtype=float)
            d = np.linalg.norm(coords - coords.mean(axis=0), axis=1)
            q1, q3 = np.percentile(d, [25, 75])
            iqr = q3 - q1 if q3 > q1 else (np.std(d) * 1.349)
            outliers = np.sum((d < q1 - 1.5*iqr) | (d > q3 + 1.5*iqr))
            return (outliers / len(d)) < (1.0 - float(threshold))
        except Exception:
            return True

    def robust_pbc_distance(self, pos1: np.ndarray, pos2: np.ndarray, box: np.ndarray | None, pbc_type: str = 'orthorhombic') -> float:
        """Robust PBC distance between two positions.

        For orthorhombic boxes: minimal image per component. For other types,
        try MDAnalysis if available; else fallback to Euclidean.
        """
        try:
            p = np.asarray(pos1, dtype=float); q = np.asarray(pos2, dtype=float)
            d = q - p
            if pbc_type == 'orthorhombic' and box is not None and len(box) >= 3:
                Lx, Ly, Lz = float(box[0]), float(box[1]), float(box[2])
                if Lx > 0: d[0] -= round(d[0]/Lx) * Lx
                if Ly > 0: d[1] -= round(d[1]/Ly) * Ly
                if Lz > 0: d[2] -= round(d[2]/Lz) * Lz
                return float(np.linalg.norm(d))
            # Try MDAnalysis
            try:
                import MDAnalysis as mda
                from MDAnalysis.lib.distances import distance_array
                aa = pos1.reshape(1, -1); bb = pos2.reshape(1, -1)
                dist = float(distance_array(aa, bb, box=box)[0,0])
                return dist
            except Exception:
                return float(np.linalg.norm(d))
        except Exception:
            return float(np.linalg.norm(np.asarray(pos2) - np.asarray(pos1)))

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
            dz = dz - round(dz / zbox) * zbox
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
        k = float(self.umbrella_config.get('force_constant', 900))
        dz_max = _dz_max_for_overlap(float(self.temperature_K), k, float(self.target_overlap))
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
                overlap = float(np.trapezoid(np.minimum(p1, p2), x=xs))
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

    def compute_window_qc(self, pullx_path: str, k_kj_mol_nm2: float, center_nm: float,
                           half_energy_tol_kj: float = 2.0):
        """Umfassendere QC: ESS (robust), Halbenergie, Drift.

        Verwirft initiale Transienten (qc.discard_fraction) vor der Auswertung.
        """
        try:
            data = np.loadtxt(pullx_path, comments=['#','@'])
        except Exception:
            return {"ess": None, "ess_pass": False, "half_energy_diff_kj": None, "half_pass": False,
                    "drift_nm_per_ns": None, "drift_pass": False}
        if data.size == 0:
            return {"ess": None, "ess_pass": False, "half_energy_diff_kj": None, "half_pass": False,
                    "drift_nm_per_ns": None, "drift_pass": False}
        t = data[:,0]  # ps
        z = data[:,1]
        # Discard initial transient for QC
        if self.qc_discard_frac > 0 and len(z) > 10:
            start_idx = int(len(z) * self.qc_discard_frac)
            start_idx = min(max(start_idx, 0), len(z)-2)
            z = z[start_idx:]
            t = t[start_idx:]
        # sampling interval (ps)
        if len(t) >= 2:
            dt_ps = float(np.median(np.diff(t)))
        else:
            dt_ps = 1.0
        # Optional stride for ESS to reduce autocorrelation bias from oversampling
        if getattr(self, 'qc_ess_stride', 1) > 1:
            stride = max(1, int(self.qc_ess_stride))
            z = z[::stride]
            t = t[::stride]
        # ESS robust
        ess = int(_robust_ess(z))
        ess_pass = ess >= int(self.qc_config.get('min_ess_frames', 200))
        # Half-time convergence proxy: harmonic bias energy difference between halves
        mid = len(z)//2
        def mean_bias(zseg):
            dz = zseg - center_nm
            # 1/2 k (dz)^2 in kJ/mol
            return 0.5 * k_kj_mol_nm2 * float(np.mean(dz*dz))
        e1 = mean_bias(z[:mid])
        e2 = mean_bias(z[mid:])
        de = abs(e1 - e2)
        half_pass = de <= half_energy_tol_kj
        # Drift-Analyse (nm/ns)
        try:
            # Linear fit z(t)
            tt_ns = t * 0.001
            A = np.vstack([tt_ns, np.ones_like(tt_ns)]).T
            slope, intercept = np.linalg.lstsq(A, z, rcond=None)[0]
            drift = float(abs(slope))  # nm/ns
        except Exception:
            drift = None
        drift_tol = float(self.qc_config.get('drift_tolerance_nm_per_ns', 0.1))
        drift_pass = (drift is not None) and (drift <= drift_tol)
        return {"ess": ess, "ess_pass": ess_pass, "half_energy_diff_kj": de, "half_pass": half_pass,
                "drift_nm_per_ns": drift, "drift_pass": drift_pass}
    
    def run_umbrella_window(self, window_center, window_dir, start_structure, 
                          replicate=1, extend_time=0):
        """Run single umbrella window with all improvements"""
        
        window_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare index for this window: reuse common index if configured
        index_file = window_dir / "index.ndx"
        if self.consistent_reference and self.common_index_file:
            try:
                shutil.copy(self.common_index_file, index_file)
                ref_group = self.common_ref_group or 'LocalPatch'
            except Exception:
                # Fallback to generating a window-specific index
                index_file, ref_group = self.create_dynamic_index_file(start_structure, index_file, window_center)
        else:
            index_file, ref_group = self.create_dynamic_index_file(
                start_structure,
                index_file,
                window_center
            )

        # Determine pbcatom for reference group
        pbcatom1 = 0
        if self.consistent_reference and getattr(self, 'common_pbcatom1_index', None):
            pbcatom1 = int(self.common_pbcatom1_index)
        else:
            gro_atoms = self.read_gro_atoms(start_structure)
            pbcatom1 = int(self._compute_pbcatom(index_file, ref_group, gro_atoms))
        # Validate reference group and pbcatom selection early to avoid unstable pulls
        ref_ids = self._parse_ndx_group(index_file, ref_group)
        if not ref_ids:
            raise RuntimeError(f"Reference group '{ref_group}' is empty or missing in index {index_file}.")
        if pbcatom1 <= 0:
            raise RuntimeError(
                "No valid pbcatom for reference group. Ensure the reference group contains atoms "
                "and that index_leaflets.ndx/patch selection are correct."
            )

        
        # Generate deterministic seed
        peptide_id = self.run_dir.name.split('_')[0]
        seed_string = f"{peptide_id}_rep{replicate}_z{window_center:.2f}"
        seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16) % 1000000
        
        # Production time
        prod_time = self.umbrella_config.get('production_ns', 60.0) + extend_time
        force_k = float(self.umbrella_config.get('force_constant', 900))
        # Optional adaptive force constant targeting desired sigma_z
        adapt_k_cfg = (self.umbrella_config.get('adaptive_k') or {})
        if bool(adapt_k_cfg.get('enabled', False)):
            sigma_target = float(adapt_k_cfg.get('sigma_target_nm', 0.06))
            k_min = float(adapt_k_cfg.get('k_min', 200.0))
            k_max = float(adapt_k_cfg.get('k_max', 5000.0))
            k_thermal = KB_KJ_MOL_K * float(self.temperature_K) / max(1e-6, sigma_target**2)
            k_opt = float(np.clip(k_thermal, k_min, k_max))
            if abs(k_opt - force_k) / force_k > 0.1:
                self.logger.info(f"  Adaptive k: {force_k:.1f} → {k_opt:.1f} kJ/mol/nm^2 (target σ≈{sigma_target} nm)")
            force_k = k_opt
        
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
        smd_prepare = bool(self.umbrella_config.get('smd_prepare', True))
        if smd_prepare and float(self.umbrella_config.get('pre_smd_ns', 1.0)) > 0:
            smd_gro = self.run_smd_ladder(Path(start_structure), float(window_center), window_dir, index_file, ref_group, pbcatom1, seed)
            if smd_gro and Path(smd_gro).exists():
                structure_for_grompp = Path(smd_gro)

        # Validate reference group consistency (warn if poor) after input structure is known
        try:
            gro_atoms_check = self.read_gro_atoms(structure_for_grompp)
            if not self._validate_reference_group(gro_atoms_check, ref_ids, threshold=0.95):
                self.logger.warning("Reference group validation failed (outliers). Consider a more global PO4 reference.")
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
        gen_vel_flag = (not has_vel) and (extend_time <= 0)
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
ns_type                = grid
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
ld-seed                = {seed}

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
pull-pbc-ref-prev-step-com = yes  ; Use COM from previous step for PBC
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
        
        # Run grompp via Docker
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
        
        print(f"  Running grompp...")
        sys.stdout.flush()
        result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
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
        # if extending, use checkpoint input and append
        cpt_path = window_dir / "umbrella.cpt"
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
            if "step" in line.lower() or "Step" in line:
                # Clear line and print progress
                sys.stdout.write(f"\r  Progress: {line.strip()[:80]}")
                sys.stdout.flush()
            # Show important messages
            elif any(keyword in line for keyword in ["WARNING", "ERROR", "Note", "Writing", "Back Off!"]):
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                if ("LINCS WARNING" in line) or ("Back Off!" in line and re.search(r"step\\d+.*\\.pdb", line)):
                    lincs_flag = True
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
        
        if process.returncode == 0 and not lincs_flag:
            print(f"✓ Window z={window_center:.3f} nm completed")
            return {
                "center": window_center,
                "pullf": str(pullf_file),
                "pullx": str(pullx_file),
                "tpr": str(tpr_file),
                "seed": seed
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
            if "step" in line.lower() or "Step" in line:
                sys.stdout.write(f"\r  Progress: {line.strip()[:80]}")
                sys.stdout.flush()
            elif any(keyword in line for keyword in ["WARNING", "ERROR", "Note", "Writing", "Back Off!"]):
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
                if ("LINCS WARNING" in line) or ("Back Off!" in line and re.search(r"step\\d+.*\\.pdb", line)):
                    lincs_flag = True
            current_time = time.time()
            if current_time - last_time > 5:
                sys.stdout.write(".")
                sys.stdout.flush()
                last_time = current_time
        process.wait()
        sys.stdout.write("\n")
        if process.returncode == 0 and not lincs_flag:
            print(f"✓ Window z={window_center:.3f} nm completed (after stabilization)")
            return {
                "center": window_center,
                "pullf": str(pullf_file),
                "pullx": str(pullx_file),
                "tpr": str(tpr_file),
                "seed": seed
            }
        print(f"✗ Window z={window_center:.3f} nm failed even after stabilization")
        # Mark this window to avoid endless retries on resume
        try:
            (window_dir / ".no_retry").write_text("unstable_after_stabilization\n")
        except Exception:
            pass
        return None
    def _extend_failed_windows(self, pmf_dir, windows_data, qc_results, start_structure, replicate):
        """Extend windows failing ESS or half-time convergence, up to max_extend_ns per window."""
        if not self.auto_extend:
            return windows_data
        # Map centers to current total extra time by inspecting existing metadata in memory
        extended = 0
        new_windows = []
        for w in windows_data:
            if not w:
                new_windows.append(w)
                continue
            center = w['center']
            # Find QC for this center
            ess_item = next((x for x in qc_results.get('ess_check', []) if x.get('center') == center), None)
            half_item = next((x for x in qc_results.get('convergence_check', []) if x.get('center') == center), None)
            ess_fail = (ess_item is None) or (not ess_item.get('passed', False))
            half_fail = (half_item is None) or (not half_item.get('passed', False))
            if not (ess_fail or half_fail):
                new_windows.append(w)
                continue
            # Determine current extension by reading log nsteps? We simply add configured extend_ns but cap at max_extend_ns
            # Store current cumulative extension in a small marker file
            window_dir = pmf_dir / "windows" / f"z_{center:+.3f}"
            marker = window_dir / ".extended_ns"
            already = 0.0
            if marker.exists():
                try:
                    already = float(marker.read_text().strip())
                except Exception:
                    already = 0.0
            if already >= self.max_extend_ns:
                print(f"  ↷ Window z={center:+.3f} reached max extension ({already} ns), not extending further.")
                new_windows.append(w)
                continue
            add = min(self.extend_ns, self.max_extend_ns - already)
            print(f"  ⤴ Extending window z={center:+.3f} by +{add} ns (ESS/half-time fail; total so far {already+add} ns)")
            res = self.run_umbrella_window(center, window_dir, start_structure, replicate, extend_time=add)
            if res:
                # overwrite existing entry with updated paths (same files) and keep the center
                new_windows.append(res)
                extended += 1
                try:
                    marker.write_text(f"{already+add}\n")
                except Exception:
                    pass
            else:
                new_windows.append(w)
        if extended:
            print(f"Extended {extended} windows due to ESS/half-time QC failures.")
        return new_windows
    
    def find_topology_file(self):
        """Find the appropriate topology file"""
        # Look in PMF system directory first
        pmf_dirs = list(self.run_dir.glob("pmf_system/replicate_*/system.top"))
        if pmf_dirs:
            return pmf_dirs[0]
        
        # Check system directory for any .top file
        system_tops = list(self.run_dir.glob("system/*.top"))
        if system_tops:
            # Prefer system_n1.top if it exists (single peptide)
            for top in system_tops:
                if "n1" in top.name:
                    return top
            # Otherwise use the first one found
            return system_tops[0]
        
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
        # Overlap between neighbors (sorted by center so adjacency is correct)
        val_windows = [w for w in windows_data if w]
        val_windows = sorted(val_windows, key=lambda x: x['center'], reverse=True)
        for i in range(len(val_windows) - 1):
            w1 = val_windows[i]
            w2 = val_windows[i+1]
            overlap = self.check_window_overlap(w1['pullx'], w2['pullx'], self.target_overlap)
            qc_results["overlap_check"].append({
                "windows": [w1['center'], w2['center']],
                "overlap": overlap,
                "passed": overlap >= self.min_overlap
            })
        # Per-window ESS & convergence
        k = float(self.umbrella_config.get('force_constant', 900))
        for w in val_windows:
            qc = self.compute_window_qc(w['pullx'], k, w['center'], self.qc_config.get('half_energy_tol_kj', 2.0))
            qc_results["ess_check"].append({"center": w['center'], "ess": qc['ess'], "passed": qc['ess_pass']})
            qc_results["convergence_check"].append({
                "center": w['center'],
                "half_energy_diff_kj": qc['half_energy_diff_kj'],
                "half_pass": qc['half_pass'],
                "drift_nm_per_ns": qc.get('drift_nm_per_ns'),
                "drift_pass": qc.get('drift_pass', False),
                "passed": qc['half_pass'] and qc.get('drift_pass', True)
            })
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
        # Windows and QC summary
        manifest['windows'] = {
            'centers': window_centers,
            'n_windows': len(window_centers),
        }
        # Persist
        run_info = pmf_dir / "RUN_INFO.yaml"
        try:
            with open(run_info, 'w') as f:
                yaml.dump(manifest, f, default_flow_style=False)
        except Exception as e:
            print(f"Warning: failed to write RUN_INFO.yaml: {e}")
        return run_info

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

    def validate_window_consistency(self, pmf_dir: Path, windows_data, tolerance=0.1):
        """Validate pairwise consistency between adjacent windows.

        Compares measured vs expected overlaps; writes YAML to pmf_dir/window_consistency.yaml.
        """
        k = float(self.umbrella_config.get('force_constant', 900))
        wins = sorted([w for w in windows_data if w], key=lambda x: x['center'], reverse=True)
        report = {"problematic_pairs": [], "passed": True}
        for a, b in zip(wins[:-1], wins[1:]):
            try:
                ol = self.check_window_overlap(a['pullx'], b['pullx'], self.target_overlap)
                dz = abs(float(a['center']) - float(b['center']))
                ol_exp = self._expected_overlap_gaussian(dz, k)
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
        with open(report_file, 'w') as f:
            yaml.dump(qc_results, f, default_flow_style=False)
        n_ol_pass = sum(1 for c in qc_results.get('overlap_check', []) if c.get('passed'))
        n_ol_total = len(qc_results.get('overlap_check', []))
        n_ess_pass = sum(1 for c in qc_results.get('ess_check', []) if c.get('passed'))
        n_ess_total = len(qc_results.get('ess_check', []))
        n_half_pass = sum(1 for c in qc_results.get('convergence_check', []) if c.get('passed'))
        n_half_total = len(qc_results.get('convergence_check', []))
        print("\n=== QC Report ===")
        print(f"Overlap checks: {n_ol_pass}/{n_ol_total} passed")
        print(f"ESS checks:     {n_ess_pass}/{n_ess_total} passed (threshold {self.qc_config.get('min_ess_frames',200)})")
        print(f"Half-time:      {n_half_pass}/{n_half_total} passed (ΔE_tol {self.qc_config.get('half_energy_tol_kj',2.0)} kJ/mol, drift_tol {self.qc_config.get('drift_tolerance_nm_per_ns',0.1)} nm/ns)")
        for check in qc_results.get('overlap_check', []):
            if not check['passed']:
                print(f"  ⚠ Low overlap ({check['overlap']:.3f}) between z={check['windows'][0]:.2f} and z={check['windows'][1]:.2f}")
        return report_file
    
    def _scan_existing_windows(self, pmf_dir: Path):
        """Return list of window dicts for existing window outputs."""
        windows_root = pmf_dir / "windows"
        if not windows_root.exists():
            return []
        wins = []
        for d in sorted(windows_root.glob("z_*")):
            pullx = d / "pullx.xvg"
            pullf = d / "pullf.xvg"
            tpr = d / "umbrella.tpr"
            if not pullx.exists():
                continue
            # parse center from folder name 'z_+0.700'
            try:
                m = re.match(r"z_([+-]?[0-9]+\.[0-9]+)", d.name)
                center = float(m.group(1)) if m else None
            except Exception:
                center = None
            if center is None:
                # try fallback: read center from mdp header
                center = None
            wins.append({
                "center": center,
                "pullf": str(pullf),
                "pullx": str(pullx),
                "tpr": str(tpr)
            })
        # Filter entries without center
        wins = [w for w in wins if w.get("center") is not None]
        # Deduplicate by center
        out = {}
        for w in wins:
            out[_round3(w['center'])] = w
        return [out[k] for k in sorted(out.keys(), reverse=True)]

    def run_pmf_calculation(self, replicate=1, tag=None, resume=False, qc_only=False):
        """Main PMF calculation workflow with all improvements"""

        # Setup directories
        pmf_dir = self.run_dir / "pmf"
        if tag:
            pmf_dir = pmf_dir / tag
        else:
            pmf_dir = pmf_dir / f"replicate_{replicate}"

        pmf_dir.mkdir(parents=True, exist_ok=True)
        self._setup_logging(pmf_dir)
        self.logger.info(f"PMF run started for {self.run_dir} -> {pmf_dir}")

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

        # Fast paths
        if qc_only:
            self.logger.info("qc-only mode: scanning existing windows and generating reports")
            windows_data = self._scan_existing_windows(pmf_dir)
            if not windows_data:
                self.logger.error("No existing windows found to analyze.")
                raise SystemExit(2)
            qc_results = self.run_qc_checks(windows_data)
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
                "patch_radius": self.umbrella_config.get('patch_radius', 2.0),
                "windows": windows_data,
                "qc_results": qc_results
            }
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

        # Deduplicate windows by center (keep last occurrence, typically most recent/extended)
        uniq = {}
        for w in windows_data:
            if not w:
                continue
            uniq[_round3(w['center'])] = w
        windows_data = [uniq[k] for k in sorted(uniq.keys(), reverse=True)]

        # QC loop with auto-densify/extend
        max_rounds = int(self.qc_config.get('max_qc_rounds', 5))
        all_centers = [w['center'] for w in windows_data if w]
        for round_idx in range(max_rounds):
            qc_results = self.run_qc_checks(windows_data)
            to_add = []
            # Add midpoint windows for low-overlap neighbors
            for check in qc_results['overlap_check']:
                if not check['passed']:
                    mid = _round3(0.5*(check['windows'][0] + check['windows'][1]))
                    if mid not in all_centers:
                        to_add.append(mid)
            # Add midpoint windows for large geometric gaps (safety)
            k = float(self.umbrella_config.get('force_constant', 900))
            dz_max = _dz_max_for_overlap(float(self.temperature_K), k, float(self.target_overlap))
            all_centers_sorted = _dedup_sorted_centers([w['center'] for w in windows_data if w])
            for a,b in zip(all_centers_sorted[:-1], all_centers_sorted[1:]):
                if abs(a-b) > dz_max*1.05:
                    mid = _round3(0.5*(a+b))
                    if mid not in all_centers:
                        to_add.append(mid)
            to_add = _dedup_sorted_centers(to_add)
            if to_add:
                self.logger.info(f"QC round {round_idx+1}: adding {len(to_add)} windows to improve coverage/overlap")
                for new_center in to_add:
                    window_dir = pmf_dir / "windows" / f"z_{new_center:+.3f}"
                    res = self.run_umbrella_window(new_center, window_dir, start_structure, replicate)
                    if res:
                        windows_data.append(res)
                        all_centers.append(new_center)
            # After densification in this round, try extending windows that failed ESS/half-time
            qc_results = self.run_qc_checks(windows_data)
            windows_data = self._extend_failed_windows(pmf_dir, windows_data, qc_results, start_structure, replicate)
            # Continue to next round; QC will be recomputed at top
            if not to_add:
                break
        # Final QC report
        qc_results = self.run_qc_checks(windows_data)
        self.generate_qc_report(qc_results, pmf_dir)
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
            "production_time": self.umbrella_config.get('production_ns', 20.0),
            "reference_mode": self.umbrella_config.get('ref_mode', 'patch'),
            "consistent_reference": bool(self.consistent_reference),
            "reference_group": self.common_ref_group if self.common_ref_group else self.umbrella_config.get('ref_mode', 'patch'),
            "index_file": str(Path(self.common_index_file).relative_to(self.run_dir)) if self.common_index_file else None,
            "patch_radius": self.umbrella_config.get('patch_radius', 2.0),
            "windows": windows_data,
            "qc_results": qc_results
        }

        metadata_file = pmf_dir / "pmf_metadata.yaml"
        with open(metadata_file, 'w') as f:
            yaml.dump(metadata, f, default_flow_style=False)

        # Final summary with more details
        successful_windows = len([w for w in windows_data if w])
        self.logger.info("✅ PMF CALCULATION COMPLETED")
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
        self.logger.info("Next step: Run MBAR analysis")
        self.logger.info(f"Command: python scripts/analysis/pmf_mbar_analysis.py {pmf_dir}")

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
    
    args = parser.parse_args()
    
    # Initialize runner
    runner = PMFRunner(args.run_dir)
    
    # Run PMF calculation
    runner.run_pmf_calculation(replicate=args.replicate, tag=args.tag, resume=args.resume, qc_only=args.qc_only)

if __name__ == "__main__":
    main()
