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

def load_config():
    """Load SOLVIA configuration"""
    config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
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
    """Find the top surface of the membrane"""
    _, _, atom_lines, _ = read_gro_file(membrane_gro)
    
    z_coords_upper = []
    for line in atom_lines:
        if len(line) > 44:
            # Look for phosphate beads in upper leaflet lipids
            if any(lipid in line for lipid in ["POPC", "PSM"]) and "PO4" in line:
                z = float(line[36:44])
                if z > 7.0:  # Upper leaflet
                    z_coords_upper.append(z)
    
    if z_coords_upper:
        return np.mean(z_coords_upper)
    else:
        # Fallback: estimate based on box
        return 8.5

def _poisson_disk_sample(n_points: int, width: float, height: float, min_dist: float, 
                         rng: np.random.Generator, max_attempts: int = 10000) -> list:
    """Simple dart-throw Poisson-disk sampling in [margin, width-margin]×[margin, height-margin]."""
    margin = min_dist
    points: list[list[float]] = []
    attempts = 0
    while len(points) < n_points and attempts < max_attempts:
        attempts += 1
        x = rng.uniform(margin, max(margin, width - margin))
        y = rng.uniform(margin, max(margin, height - margin))
        ok = True
        for px, py in points:
            if (x - px) ** 2 + (y - py) ** 2 < (min_dist ** 2):
                ok = False
                break
        if ok:
            points.append([float(x), float(y)])
    return points


def calculate_peptide_positions(n_peptides, box_size, spacing, placement: str, rng: np.random.Generator):
    """Calculate XY positions for peptide placement.

    - For 1 peptide: place at box center.
    - For known counts (8, 12, 16): use predefined grids.
    - For other counts: distribute on a simple grid within [1.5, box-1.5].
    """
    positions = []
    
    if placement == "poisson" and n_peptides > 0:
        pts = _poisson_disk_sample(n_peptides, box_size[0], box_size[1], spacing, rng)
        # If not enough points found, fall back to grid below
        if len(pts) >= n_peptides:
            return pts[:n_peptides]

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

def insert_peptides(run_dir, occupancy="low", n_peptides_override=None, placement: str = "poisson", orientation: str = "random", replicates: int = 1):
    """Insert peptides above membrane"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Get number of peptides based on occupancy or override
    if n_peptides_override is not None and n_peptides_override > 0:
        n_peptides = int(n_peptides_override)
        occupancy_label = f"n{n_peptides}"
    else:
        n_peptides = config['peptide_insertion']['occupancy_levels'][occupancy]
        occupancy_label = occupancy
    
    # Input files
    peptide_gro = os.path.join(run_dir, "cg_pdb", f"{metadata['peptide_id']}_cg.pdb")
    membrane_gro = os.path.join(run_dir, "membrane_template", "membrane_template.gro")
    
    # Check if membrane exists, if not create it
    if not os.path.exists(membrane_gro):
        print("Membrane template not found. Building it first...")
        import subprocess
        subprocess.run([
            sys.executable,
            os.path.join(Path(__file__).parent, "04_build_membrane.py"),
            run_dir
        ], check=True)
    
    # Output directory
    output_dir = os.path.join(run_dir, "system")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Inserting {n_peptides} peptide(s) ({occupancy_label}), replicates={replicates}...")
    
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
    
    print(f"Membrane top at z={membrane_top:.2f} nm")
    print(f"Peptide dimensions: {peptide_dims[0]:.2f} x {peptide_dims[1]:.2f} x {peptide_dims[2]:.2f} nm")
    
    last_gro, last_top = None, None
    for rep in range(1, max(1, replicates) + 1):
        rep_tag = f"{occupancy_label}_rep{rep}" if replicates > 1 else occupancy_label
        output_gro = os.path.join(output_dir, f"system_{rep_tag}.gro")
        output_top = os.path.join(output_dir, f"system_{rep_tag}.top")

        if os.path.exists(output_gro) and os.path.exists(output_top):
            print(f"System already exists: {rep_tag}")
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
            rng
        )
        z_position = membrane_top + config['peptide_insertion']['distance_from_membrane']

        # Orientation handling (tilt around x-axis relative to membrane plane z)
        tilt_deg = _orientation_to_angle(orientation, rng)
        rot = _rotation_matrix(tilt_deg, axis='x')

        # Place peptides
        all_atoms = []
        atom_counter = 0

        # Add peptides first (to match topology order)
        for i, (x, y) in enumerate(positions):
            print(f"  [{rep_tag}] Placing peptide {i+1} at ({x:.1f}, {y:.1f}, {z_position:.1f}); tilt={tilt_deg:.1f}°")
            
            for line in peptide_atoms:
                if len(line) > 44:
                    # Parse atom line
                    resid = int(line[0:5])
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

                    # Rotate around x-axis by tilt
                    rx, ry, rz = rot @ np.array([vx, vy, vz])

                    # Translate to target position (x,y,z)
                    new_x = rx + x
                    new_y = ry + y
                    new_z = rz + z_position
                    
                    # Write new atom line
                    new_line = f"{resid:5d}{resname}{atomname}{atomid:5d}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}\n"
                    all_atoms.append(new_line)
                    atom_counter += 1
        
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
        print(f"✓ Created system: {output_gro}")

        # Create topology for this replicate
        create_system_topology(run_dir, metadata['peptide_id'], n_peptides, rep_tag)
        # Move to replicate-specific filename if created under occupancy-only name
        default_top = os.path.join(output_dir, f"system_{occupancy_label}.top")
        if os.path.exists(default_top) and default_top != output_top:
            try:
                shutil.move(default_top, output_top)
            except Exception:
                pass
        last_gro, last_top = output_gro, output_top

    return last_gro, last_top

def convert_pdb_to_gro(pdb_file, gro_file):
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
    
    # Write GRO
    with open(gro_file, 'w') as f:
        f.write("Converted from PDB\n")
        f.write(f"{len(atoms):5d}\n")
        for atom in atoms:
            f.write(atom)
        f.write("  10.000  10.000  14.000\n")  # Default box

def create_system_topology(run_dir, peptide_id, n_peptides, tag):
    """Create system topology file"""
    config = load_config()
    
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
    peptide_itp = os.path.join(run_dir, "cg_pdb", f"{peptide_id}.itp")
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
    except Exception:
        pass

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
#include "../cg_pdb/{peptide_id}.itp"

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
    
    print(f"✓ Created topology: {output_top}")

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
        choices=["poisson", "grid"],
        default="poisson",
        help="XY-Platzierung: 'poisson' (min spacing) oder 'grid'"
    )
    parser.add_argument(
        "--orientation",
        choices=["random", "parallel", "tilt45", "perpendicular"],
        default="random",
        help="Orientierungspool: parallel / 45° / senkrecht / random"
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Anzahl unabhängiger Replikate (separate Systeme mit unterschiedlichen Seeds)"
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
        replicates=args.replicates
    )
    
    print(f"\nNext step: Equilibrate system")
    print(f"Command: python 06_equilibrate.py {args.run_dir} --occupancy {args.occupancy}")

if __name__ == "__main__":
    main()
