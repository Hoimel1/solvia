#!/usr/bin/env python3
"""
SOLVIA Peptide Insertion
Inserts peptides in controlled 3-2-3 arrangement above membrane
"""

import os
import sys
import yaml
import shutil
import argparse
import numpy as np
import math
import zlib
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def load_config(config_path=None):
    """Load SOLVIA configuration"""
    if config_path is None:
        config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
    config_path = Path(config_path)
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def read_gro_file(filename):
    """Read a GRO file and return its components"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    header = lines[0]
    n_atoms = int(lines[1])
    atom_lines = lines[2:2+n_atoms]
    box_line = lines[-1]
    
    return header, n_atoms, atom_lines, box_line

def write_gro_file(filename, header, atom_lines, box_line):
    """Write a GRO file"""
    with open(filename, 'w') as f:
        f.write(header)
        f.write(f"{len(atom_lines):5d}\n")
        for line in atom_lines:
            f.write(line)
        f.write(box_line)

def get_peptide_dimensions(peptide_gro):
    """Calculate peptide dimensions from GRO file"""
    _, _, atom_lines, _ = read_gro_file(peptide_gro)
    
    coords = []
    for line in atom_lines:
        if len(line) > 44:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            coords.append([x, y, z])
    
    coords = np.array(coords)
    dimensions = coords.max(axis=0) - coords.min(axis=0)
    center = coords.mean(axis=0)
    
    return dimensions, center

def find_membrane_top(membrane_gro):
    """Find the top surface (upper leaflet) z using PO4 positions dynamically.

    Uses median PO4 z to split leaflets; returns mean z of upper-leaflet PO4s.
    Falls back to box mid-plane if PO4 not found.
    """
    _, _, atom_lines, box_line = read_gro_file(membrane_gro)
    po4_z = []
    for line in atom_lines:
        if len(line) > 44 and 'PO4' in line:
            try:
                z = float(line[36:44])
                po4_z.append(z)
            except Exception:
                continue
    if po4_z:
        zs = sorted(po4_z)
        zmed = zs[len(zs)//2]
        upper = [z for z in po4_z if z >= zmed]
        if upper:
            return float(np.mean(upper))
    # Fallback: estimate based on box height
    try:
        box_z = float(box_line.split()[2])
        return box_z * 0.5
    except Exception as e:
        logging.warning(f"Could not parse box line for fallback membrane top: {e}")
        return 8.5

def _pbc_d2(x: float, y: float, px: float, py: float, width: float, height: float) -> float:
    dx = abs(x - px)
    dy = abs(y - py)
    dx = min(dx, width - dx)
    dy = min(dy, height - dy)
    return dx*dx + dy*dy

def _poisson_disk_sample(n_points: int, width: float, height: float, min_dist: float, 
                         rng: np.random.Generator, max_attempts: int = 10000) -> list:
    """Simple Poisson-disk sampling with PBC-aware distance in XY.

    Samples in [margin, width-margin] × [margin, height-margin].
    """
    margin = min_dist
    points: list[list[float]] = []
    attempts = 0
    r2 = (min_dist ** 2)
    while len(points) < n_points and attempts < max_attempts:
        attempts += 1
        x = float(rng.uniform(margin, max(margin, width - margin)))
        y = float(rng.uniform(margin, max(margin, height - margin)))
        ok = True
        for px, py in points:
            if _pbc_d2(x, y, px, py, width, height) < r2:
                ok = False
                break
        if ok:
            points.append([x, y])
    return points


def calculate_peptide_positions(n_peptides, box_size, spacing, placement: str, rng: np.random.Generator,
                                *, run_dir: str | None = None, config: dict | None = None):
    """Calculate XY positions for peptide placement.

    - For 1 peptide: place at box center.
    - For known counts (8, 12, 16): use predefined grids.
    - For other counts: distribute on a simple grid within [1.5, box-1.5].
    """
    positions = []
    
    if placement in ("poisson", "cooperative", "thermo") and n_peptides > 0:
        pts = _poisson_disk_sample(n_peptides, box_size[0], box_size[1], spacing, rng)
        # If not enough points found, fall back to grid below
        if len(pts) >= n_peptides:
            pts = pts[:n_peptides]
            if placement in ("cooperative", "thermo") and n_peptides > 1:
                # Cooperative/thermo refinement
                d_opt = 1.8
                cutoff = max(2.0, spacing)
                steps = 2000
                step_sigma = 0.1
                scale = 1.0
                sigma_frac = 0.15
                rep_exp = 12.0
                rep_scale = 1.0
                att_scale = 1.0
                if placement == "thermo":
                    # Sequence-based scaling (more hydrophobic => stronger cooperative attraction)
                    seq = _load_peptide_sequence(run_dir) if run_dir else None
                    energies = _get_aa_transfer_energies(config or {})
                    delta_g = _sequence_hydrophobic_transfer(seq, energies)
                    # Map delta_g to a positive scale: more negative -> larger scale
                    # Heuristic: scale = clip(1 + (-delta_g)/50, 0.5, 3.0)
                    try:
                        scale = float(np.clip(1.0 + (-delta_g) / 50.0, 0.5, 3.0))
                    except Exception:
                        scale = 1.0
                    # Allow config overrides
                    thermo_cfg = (config or {}).get('peptide_insertion', {}).get('thermo', {})
                    d_opt = float(thermo_cfg.get('d_opt', d_opt))
                    cutoff = float(thermo_cfg.get('cutoff', cutoff))
                    steps = int(thermo_cfg.get('mc_steps', steps))
                    step_sigma = float(thermo_cfg.get('step_sigma', step_sigma))
                    sigma_frac = float(thermo_cfg.get('sigma_frac', sigma_frac))
                    rep_exp = float(thermo_cfg.get('repulsion_exponent', rep_exp))
                    rep_scale = float(thermo_cfg.get('repulsion_scale', rep_scale))
                    att_scale = float(thermo_cfg.get('attraction_scale', att_scale))
                else:
                    coop_cfg = (config or {}).get('peptide_insertion', {}).get('cooperative', {})
                    d_opt = float(coop_cfg.get('d_opt', d_opt))
                    cutoff = float(coop_cfg.get('cutoff', cutoff))
                    steps = int(coop_cfg.get('mc_steps', steps))
                    step_sigma = float(coop_cfg.get('step_sigma', step_sigma))
                    sigma_frac = float(coop_cfg.get('sigma_frac', sigma_frac))
                    rep_exp = float(coop_cfg.get('repulsion_exponent', rep_exp))
                    rep_scale = float(coop_cfg.get('repulsion_scale', rep_scale))
                    att_scale = float(coop_cfg.get('attraction_scale', att_scale))
                return _adjust_positions_cooperative(
                    pts, (box_size[0], box_size[1]), rng,
                    d_opt=d_opt, cutoff=cutoff, d_min=spacing,
                    steps=steps, step_sigma=step_sigma, scale=scale,
                    sigma_frac=sigma_frac, rep_exp=rep_exp,
                    rep_scale=rep_scale, att_scale=att_scale
                )
            return pts

    if n_peptides == 1:
        positions = [[box_size[0] / 2.0, box_size[1] / 2.0]]
    elif n_peptides == 8:
        # 3-2-3 arrangement
        rows = [3, 2, 3]
        y_positions = [2.5, 5.0, 7.5]  # nm
        
        for row_idx, (n_in_row, y_pos) in enumerate(zip(rows, y_positions)):
            if n_in_row == 3:
                x_positions = [2.5, 5.0, 7.5]
            else:  # n_in_row == 2
                x_positions = [3.5, 6.5]  # Offset for staggered arrangement
            
            for x_pos in x_positions:
                positions.append([x_pos, y_pos])
    
    elif n_peptides == 12:
        # 4-4-4 arrangement
        for y in [2.0, 4.0, 6.0, 8.0]:
            for x in [2.0, 4.0, 6.0, 8.0]:
                positions.append([x, y])
                if len(positions) >= n_peptides:
                    break
    
    elif n_peptides == 16:
        # 4x4 grid
        for y in np.linspace(1.5, 8.5, 4):
            for x in np.linspace(1.5, 8.5, 4):
                positions.append([x, y])
    
    else:
        # Generic fallback grid
        grid_n = int(np.ceil(np.sqrt(n_peptides)))
        xs = np.linspace(1.5, max(1.6, box_size[0] - 1.5), grid_n)
        ys = np.linspace(1.5, max(1.6, box_size[1] - 1.5), grid_n)
        for y in ys:
            for x in xs:
                positions.append([float(x), float(y)])
                if len(positions) >= n_peptides:
                    break
            if len(positions) >= n_peptides:
                break

    return positions[:n_peptides]


def _rotation_matrix(angle_deg: float, axis: str) -> np.ndarray:
    angle = math.radians(angle_deg)
    if axis == 'x':
        return np.array([[1, 0, 0],
                         [0, math.cos(angle), -math.sin(angle)],
                         [0, math.sin(angle),  math.cos(angle)]], dtype=float)
    if axis == 'y':
        return np.array([[ math.cos(angle), 0, math.sin(angle)],
                         [0, 1, 0],
                         [-math.sin(angle), 0, math.cos(angle)]], dtype=float)
    # default z
    return np.array([[math.cos(angle), -math.sin(angle), 0],
                     [math.sin(angle),  math.cos(angle), 0],
                     [0, 0, 1]], dtype=float)


def _orientation_to_angle(orientation: str, rng: np.random.Generator) -> float:
    orientation = (orientation or "random").lower()
    if orientation == "parallel":
        return 0.0
    if orientation in ("45", "tilt45", "tilt_45"):
        return 45.0
    if orientation in ("perpendicular", "orthogonal", "90"):
        return 90.0
    # random choice from pools
    return float(rng.choice([0.0, 45.0, 90.0]))

def _sample_tilt_roll(orientation: str, rng: np.random.Generator) -> tuple[float, float]:
    """Sample tilt and roll angles.

    - Discrete legacy modes return fixed tilt and zero roll.
    - "continuous": tilt ~ truncated N(30°, 20°) within [0, 90°], roll ~ U(0, 360°).
    - "random" (legacy): tilt from {0,45,90} with roll uniform.
    """
    o = (orientation or "random").lower()
    if o in ("parallel", "45", "tilt45", "tilt_45", "perpendicular", "orthogonal", "90"):
        return _orientation_to_angle(o, rng), 0.0
    if o == "continuous":
        mu, sigma = 30.0, 20.0
        tilt = None
        for _ in range(50):
            t = float(rng.normal(mu, sigma))
            if 0.0 <= t <= 90.0:
                tilt = t
                break
        if tilt is None:
            tilt = float(np.clip(mu, 0.0, 90.0))
        roll = float(rng.uniform(0.0, 360.0))
        return tilt, roll
    # legacy random pool with random roll
    return float(rng.choice([0.0, 45.0, 90.0])), float(rng.uniform(0.0, 360.0))

def _rotation_matrix_general(tilt_deg: float, roll_deg: float) -> np.ndarray:
    """Rotation matrix R = Rz(roll) @ Ry(tilt) with Z the membrane normal."""
    tilt = math.radians(tilt_deg)
    roll = math.radians(roll_deg)
    Ry = np.array([[ math.cos(tilt), 0, math.sin(tilt)],
                   [ 0,             1, 0            ],
                   [-math.sin(tilt), 0, math.cos(tilt)]], dtype=float)
    Rz = np.array([[math.cos(roll), -math.sin(roll), 0],
                   [math.sin(roll),  math.cos(roll), 0],
                   [0,               0,              1]], dtype=float)
    return Rz @ Ry

def _pair_energy(d: float, d_opt: float, d_min: float, cutoff: float,
                 sigma_frac: float = 0.15, rep_exp: float = 12.0,
                 rep_scale: float = 1.0, att_scale: float = 1.0) -> float:
    """Physically motivated pair potential (repulsion + short-range attraction).

    - Strong repulsion for d < d_min (soft 12-power).
    - Shallow Gaussian well centered at d_opt within cutoff.
    """
    if d <= 1e-6:
        return 1e6
    V_rep = 0.0
    if d < d_min:
        ratio = max(1e-6, d_min / d)
        V_rep = rep_scale * (ratio ** rep_exp)
    V_att = 0.0
    if d <= cutoff:
        sigma = max(0.01, sigma_frac * d_opt)
        V_att = -att_scale * math.exp(-0.5 * ((d - d_opt) / sigma) ** 2)
    return V_rep + V_att

def _adjust_positions_cooperative(positions: list[list[float]], box_size_xy: tuple[float, float],
                                  rng: np.random.Generator,
                                  d_opt: float = 1.8,
                                  cutoff: float = 2.5,
                                  d_min: float = 2.5,
                                  steps: int = 2000,
                                  step_sigma: float = 0.1,
                                  scale: float = 1.0,
                                  sigma_frac: float = 0.15,
                                  rep_exp: float = 12.0,
                                  rep_scale: float = 1.0,
                                  att_scale: float = 1.0) -> list[list[float]]:
    if len(positions) <= 1:
        return positions
    W, H = float(box_size_xy[0]), float(box_size_xy[1])
    margin = float(d_min)
    kbT = 0.008314 * 310.0
    pos = np.array(positions, dtype=float)
    def e_for_idx(idx: int, arr: np.ndarray) -> float:
        e = 0.0
        xi, yi = arr[idx]
        for j in range(len(arr)):
            if j == idx:
                continue
            xj, yj = arr[j]
            dx = abs(xi - xj); dy = abs(yi - yj)
            dx = min(dx, W - dx); dy = min(dy, H - dy)
            d = math.hypot(dx, dy)
            e += scale * _pair_energy(d, d_opt, d_min, cutoff,
                                      sigma_frac=sigma_frac, rep_exp=rep_exp,
                                      rep_scale=rep_scale, att_scale=att_scale)
        return e
    for _ in range(max(0, steps)):
        i = int(rng.integers(0, len(pos)))
        e_old = e_for_idx(i, pos)
        dx, dy = rng.normal(0.0, step_sigma, size=2)
        new_x = float((pos[i, 0] + dx) % W)
        new_y = float((pos[i, 1] + dy) % H)
        old = pos[i].copy()
        pos[i, 0], pos[i, 1] = new_x, new_y
        e_new = e_for_idx(i, pos)
        dE = e_new - e_old
        if not (dE <= 0.0 or float(rng.random()) < math.exp(-dE / kbT)):
            pos[i] = old
    return [[float(x), float(y)] for x, y in pos]

def _load_peptide_sequence(run_dir: str) -> str | None:
    """Load peptide sequence from FASTA in run_dir/input/peptide.fasta"""
    fasta = os.path.join(run_dir, "input", "peptide.fasta")
    if not os.path.exists(fasta):
        return None
    seq_parts = []
    with open(fasta, 'r') as f:
        for line in f:
            if not line:
                continue
            if line.startswith('>'):
                continue
            seq_parts.append(line.strip().upper())
    seq = ''.join(seq_parts)
    return seq if seq else None

def _get_aa_transfer_energies(config: dict) -> dict:
    """Return per-residue transfer energies (kJ/mol). Config may override."""
    # Approximate Wimley-White interface-like scale in kJ/mol (sign: negative favors transfer)
    default_map = {
        'A': 0.7,  'R': 12.3, 'N': 8.6,  'D': 10.5, 'C': -0.9,
        'Q': 9.2,  'E': 8.2,  'G': 1.0,  'H': 4.6,  'I': -4.5,
        'L': -4.2, 'K': 10.0, 'M': -2.5, 'F': -7.8, 'P': 3.0,
        'S': 3.7,  'T': 2.7,  'W': -12.3,'Y': -5.3,'V': -3.5
    }
    try:
        return config.get('peptide_insertion', {}).get('thermo', {}).get('energies_kj', default_map) or default_map
    except Exception:
        return default_map

def _sequence_hydrophobic_transfer(seq: str, energies: dict) -> float:
    """Sum of per-residue transfer energies (kJ/mol); more negative => more favorable."""
    return float(sum(energies.get(aa, 0.0) for aa in seq)) if seq else 0.0

def insert_peptides(run_dir, occupancy="low", n_peptides_override=None, placement: str = "poisson", orientation: str = "random", replicates: int = 1, config_path=None):
    """Insert peptides above membrane"""
    config = load_config(config_path)
    metadata = load_run_metadata(run_dir)
    
    # Get number of peptides based on occupancy or override
    if n_peptides_override is not None and n_peptides_override > 0:
        n_peptides = int(n_peptides_override)
        occupancy_label = f"n{n_peptides}"
    else:
        n_peptides = config['peptide_insertion']['occupancy_levels'][occupancy]
        occupancy_label = occupancy
    
    # Input files
    peptide_gro = os.path.join(run_dir, "coarse_grain", f"{metadata['peptide_id']}_cg.pdb")
    membrane_gro = os.path.join(run_dir, "membrane_template", "membrane_template.gro")
    
    # Check if membrane exists, if not create it
    if not os.path.exists(membrane_gro):
        logging.info("Membrane template not found. Building it first...")
        import subprocess
        subprocess.run([
            sys.executable,
            os.path.join(Path(__file__).parent, "04_build_membrane.py"),
            run_dir
        ], check=True)
    
    # Output directory
    output_dir = os.path.join(run_dir, "system")
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Inserting {n_peptides} peptide(s) ({occupancy_label}), replicates={replicates}...")
    
    # Convert PDB to GRO if needed
    if peptide_gro.endswith('.pdb'):
        peptide_gro_temp = peptide_gro.replace('.pdb', '.gro')
        convert_pdb_to_gro(peptide_gro, peptide_gro_temp)
        peptide_gro = peptide_gro_temp
    
    # Read files
    membrane_header, membrane_natoms, membrane_atoms, box_line = read_gro_file(membrane_gro)
    peptide_header, peptide_natoms, peptide_atoms, _ = read_gro_file(peptide_gro)
    
    # Get dimensions
    peptide_dims, peptide_center = get_peptide_dimensions(peptide_gro)
    membrane_top = find_membrane_top(membrane_gro)
    box_size = [float(x) for x in box_line.split()][:3]
    
    logging.info(f"Membrane top at z={membrane_top:.2f} nm")
    logging.info(f"Peptide dimensions: {peptide_dims[0]:.2f} x {peptide_dims[1]:.2f} x {peptide_dims[2]:.2f} nm")
    
    last_gro, last_top = None, None
    for rep in range(1, max(1, replicates) + 1):
        rep_tag = f"{occupancy_label}_rep{rep}" if replicates > 1 else occupancy_label
        output_gro = os.path.join(output_dir, f"system_{rep_tag}.gro")
        output_top = os.path.join(output_dir, f"system_{rep_tag}.top")

        if os.path.exists(output_gro) and os.path.exists(output_top):
            logging.info(f"System already exists: {rep_tag}")
            last_gro, last_top = output_gro, output_top
            continue

        # RNG seed (deterministic per peptide/run/rep)
        seed_str = f"{metadata['peptide_id']}|{n_peptides}|{occupancy_label}|rep{rep}"
        seed = zlib.crc32(seed_str.encode('utf-8')) & 0xFFFFFFFF
        rng = np.random.default_rng(seed)

        # Calculate positions
        positions = calculate_peptide_positions(
            n_peptides,
            box_size,
            config['peptide_insertion']['minimum_spacing'],
            placement,
            rng,
            run_dir=run_dir,
            config=config
        )
        # Base z position above membrane
        z_position_base = membrane_top + config['peptide_insertion']['distance_from_membrane']
        # Optional z-jitter (per-replicate or per-peptide)
        zj_cfg = (config.get('peptide_insertion') or {}).get('z_jitter', {})
        zj_mode = str(zj_cfg.get('mode', 'none')).lower()
        per_peptide = bool(zj_cfg.get('per_peptide', False))
        z_position = z_position_base
        if not per_peptide:
            if zj_mode == 'normal':
                std = float(zj_cfg.get('std_nm', 0.0))
                if std > 0:
                    z_position = z_position_base + float(rng.normal(0.0, std))
            elif zj_mode == 'uniform':
                half = float(zj_cfg.get('range_nm', 0.0))
                if half > 0:
                    z_position = z_position_base + float(rng.uniform(-half, half))

        # Orientation handling (per-peptide tilt/roll relative to membrane plane)
        legacy_tilt = None
        legacy_rot = None
        if orientation.lower() in ("parallel", "45", "tilt45", "tilt_45", "perpendicular", "orthogonal", "90"):
            legacy_tilt = _orientation_to_angle(orientation, rng)
            legacy_rot = _rotation_matrix_general(legacy_tilt, 0.0)

        # Place peptides
        all_atoms = []
        atom_counter = 0

        # Add peptides first (to match topology order)
        for i, (x, y) in enumerate(positions):
            z_here = z_position
            if per_peptide:
                if zj_mode == 'normal':
                    std = float(zj_cfg.get('std_nm', 0.0))
                    if std > 0:
                        z_here = z_position_base + float(rng.normal(0.0, std))
                elif zj_mode == 'uniform':
                    half = float(zj_cfg.get('range_nm', 0.0))
                    if half > 0:
                        z_here = z_position_base + float(rng.uniform(-half, half))
            if legacy_tilt is not None:
                tilt_deg, roll_deg = legacy_tilt, 0.0
                rot = legacy_rot
            else:
                tilt_deg, roll_deg = _sample_tilt_roll(orientation, rng)
                rot = _rotation_matrix_general(tilt_deg, roll_deg)
            logging.info(f"  [{rep_tag}] Placing peptide {i+1} at ({x:.1f}, {y:.1f}, {z_here:.1f}); tilt={tilt_deg:.1f}°, roll={roll_deg:.1f}°")
            
            for line in peptide_atoms:
                if len(line) > 44:
                    # Parse atom line
                    try:
                        orig_resid = int(line[0:5])
                    except Exception as e:
                        logging.warning(f"Could not parse residue ID in peptide atom line: {e}")
                        orig_resid = 1
                    resid = i * 1000 + orig_resid
                    resname = line[5:10]
                    atomname = line[10:15]
                    atomid = atom_counter + 1
                    
                    # Get original coordinates
                    orig_x = float(line[20:28])
                    orig_y = float(line[28:36])
                    orig_z = float(line[36:44])

                    # Translate to peptide-centered coords
                    vx = orig_x - peptide_center[0]
                    vy = orig_y - peptide_center[1]
                    vz = orig_z - peptide_center[2]

                    # Rotate by combined tilt/roll
                    rx, ry, rz = rot @ np.array([vx, vy, vz])

                    # Translate to target position (x,y,z)
                    new_x = rx + x
                    new_y = ry + y
                    new_z = rz + z_here
                    
                    # Write new atom line
                    new_line = f"{resid:5d}{resname}{atomname}{atomid:5d}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}\n"
                    all_atoms.append(new_line)
                    atom_counter += 1
        
        # Simple steric validation (XY footprint) before writing
        try:
            footprint_r = 0.5 * float(max(peptide_dims[0], peptide_dims[1]))
            min_xy = max(config['peptide_insertion']['minimum_spacing'], 2.0 * footprint_r * 0.9)
            W, H = box_size[0], box_size[1]
            for a in range(len(positions)):
                xa, ya = positions[a]
                for b in range(a + 1, len(positions)):
                    xb, yb = positions[b]
                    dx = abs(xa - xb); dy = abs(ya - yb)
                    dx = min(dx, W - dx); dy = min(dy, H - dy)
                    dxy = math.hypot(dx, dy)
                    if dxy < min_xy:
                        logging.warning(f"Potential steric overlap between peptides {a+1} and {b+1} (d_xy={dxy:.2f} nm < {min_xy:.2f} nm)")
        except Exception:
            pass

        # Add membrane atoms
        for line in membrane_atoms:
            if len(line) > 44:
                atomid = atom_counter + 1
                atomid_str = f"{atomid:5d}"
                line = line[:15] + atomid_str + line[20:]
                all_atoms.append(line)
                atom_counter += 1
        
        # Write combined GRO file
        header = f"System with {n_peptides} {metadata['peptide_id']} peptides; {rep_tag}\n"
        write_gro_file(output_gro, header, all_atoms, box_line)
        logging.info(f"✓ Created system: {output_gro}")

        # Create topology for this replicate
        create_system_topology(run_dir, metadata['peptide_id'], n_peptides, rep_tag, config_path=config_path)
        # Move to replicate-specific filename if created under occupancy-only name
        default_top = os.path.join(output_dir, f"system_{occupancy_label}.top")
        if os.path.exists(default_top) and default_top != output_top:
            try:
                shutil.move(default_top, output_top)
            except Exception as e:
                logging.warning(f"Could not move topology file: {e}")
        last_gro, last_top = output_gro, output_top

    return last_gro, last_top

def convert_pdb_to_gro(pdb_file, gro_file, box_line=None):
    """Convert PDB to GRO format (simple conversion)"""
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()
    
    atoms = []
    for line in pdb_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Parse PDB
            atomid = int(line[6:11])
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            resid = int(line[22:26])
            x = float(line[30:38]) / 10.0  # Convert Å to nm
            y = float(line[38:46]) / 10.0
            z = float(line[46:54]) / 10.0
            
            # Format as GRO
            gro_line = f"{resid:5d}{resname:5s}{atomname:>5s}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
            atoms.append(gro_line)
    
    if box_line is not None:
        if not box_line.endswith('\n'):
            box_line += '\n'
    else:
        box_line = "  10.000  10.000  14.000\n"  # Default box
    # Write GRO
    with open(gro_file, 'w') as f:
        f.write("Converted from PDB\n")
        f.write(f"{len(atoms):5d}\n")
        for atom in atoms:
            f.write(atom)
        f.write(box_line)

def create_system_topology(run_dir, peptide_id, n_peptides, tag, config_path=None):
    """Create system topology file"""
    config = load_config(config_path)
    
    # Read membrane composition from template topology
    template_top = os.path.join(run_dir, "membrane_template", "membrane_template.top")
    membrane_molecules = []
    
    with open(template_top, 'r') as f:
        in_molecules = False
        for line in f:
            line = line.strip()
            if line == '[ molecules ]':
                in_molecules = True
                continue
            if in_molecules and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    mol_count = parts[1]
                    if mol_name not in ['INSANE!']:  # Skip system name
                        membrane_molecules.append((mol_name, int(mol_count)))
    
    # Determine peptide moleculetype name from ITP
    peptide_itp = os.path.join(run_dir, "coarse_grain", f"{peptide_id}.itp")
    peptide_moltype = peptide_id
    try:
        with open(peptide_itp, 'r') as f:
            lines = [l.strip() for l in f]
        for i, l in enumerate(lines):
            if l.startswith('[') and 'moleculetype' in l:
                # next non-empty, non-comment line
                for j in range(i+1, len(lines)):
                    s = lines[j]
                    if not s or s.startswith(';'):
                        continue
                    peptide_moltype = s.split()[0]
                    raise StopIteration
    except StopIteration:
        pass
    except Exception as e:
        logging.warning(f"Could not determine peptide moleculetype: {e}")

    # Create topology
    top_content = f"""; System topology for {peptide_id} with tag {tag}
; Generated by SOLVIA peptide insertion script

; Include force field
#include "{config['directories']['force_fields']}/martini_v3.0.0.itp"
#include "{config['directories']['force_fields']}/martini_v3.0.0_ffbonded_v2.itp"
#include "{config['directories']['force_fields']}/martini_v3.0.0_phospholipids_v1.itp"
#include "{config['directories']['force_fields']}/martini_v3.0.0_phospholipids_SM_v2.itp"
#include "{config['directories']['force_fields']}/martini_v3.0_sterols_v1.0.itp"
#include "{config['directories']['force_fields']}/martini_v3.0.0_solvents_v1.itp"
#include "{config['directories']['force_fields']}/martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "../coarse_grain/{peptide_id}.itp"

[ system ]
{n_peptides} {peptide_id} peptides in RBC membrane

[ molecules ]
; Peptides first (to match GRO order)
{peptide_moltype:12s} {n_peptides}
"""
    
    # Add membrane molecules
    for mol_name, mol_count in membrane_molecules:
        top_content += f"{mol_name:8s} {mol_count:6d}\n"
    
    # Write topology
    output_top = os.path.join(run_dir, "system", f"system_{tag}.top")
    with open(output_top, 'w') as f:
        f.write(top_content)
    
    logging.info(f"✓ Created topology: {output_top}")

def main():
    parser = argparse.ArgumentParser(
        description="Insert peptides into membrane for SOLVIA"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory"
    )
    parser.add_argument(
        "--occupancy",
        choices=["low", "medium", "high"],
        default="low",
        help="Peptide occupancy level (default: low)"
    )
    parser.add_argument(
        "--n-peptides",
        type=int,
        help="Override number of peptides to insert (e.g., 1 for single peptide)"
    )
    parser.add_argument(
        "--placement",
        choices=["poisson", "grid", "cooperative", "thermo"],
        default="poisson",
        help="XY placement: 'poisson' (min spacing), 'grid', 'cooperative' (MC), or 'thermo' (sequence-weighted MC)"
    )
    parser.add_argument(
        "--orientation",
        choices=["random", "parallel", "tilt45", "perpendicular", "continuous"],
        default="random",
        help="Orientation model: legacy pool (parallel/45/perpendicular/random) or 'continuous' (tilt~N(30°,20°), roll~U)"
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Anzahl unabhängiger Replikate (separate Systeme mit unterschiedlichen Seeds)"
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Optional path to configuration YAML (default: config/solvia_config.yaml)"
    )
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Insert peptides
    system_gro, system_top = insert_peptides(
        args.run_dir,
        args.occupancy,
        n_peptides_override=args.n_peptides,
        placement=args.placement,
        orientation=args.orientation,
        replicates=args.replicates,
        config_path=args.config
    )
    
    print(f"\nNext step: Equilibrate system")
    # Determine the correct tag for the equilibration command
    if args.n_peptides is not None and args.n_peptides > 0:
        eq_tag = f"--tag n{args.n_peptides}"
    else:
        eq_tag = f"--occupancy {args.occupancy}"
    print(f"Command: python 06_equilibrate.py {args.run_dir} {eq_tag}")

if __name__ == "__main__":
    main()
