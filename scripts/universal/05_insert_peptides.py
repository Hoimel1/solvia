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

def calculate_peptide_positions(n_peptides, box_size, spacing):
    """Calculate positions for 3-2-3 arrangement"""
    positions = []
    
    if n_peptides == 8:
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
    
    return positions[:n_peptides]

def insert_peptides(run_dir, occupancy="low"):
    """Insert peptides above membrane"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Get number of peptides based on occupancy
    n_peptides = config['peptide_insertion']['occupancy_levels'][occupancy]
    
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
    
    # Output files
    output_dir = os.path.join(run_dir, "system")
    output_gro = os.path.join(output_dir, f"system_{occupancy}.gro")
    output_top = os.path.join(output_dir, f"system_{occupancy}.top")
    
    # Check if already done
    if os.path.exists(output_gro) and os.path.exists(output_top):
        print(f"System with {occupancy} occupancy already exists.")
        return output_gro, output_top
    
    print(f"Inserting {n_peptides} peptides ({occupancy} occupancy)...")
    
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
    
    # Calculate positions
    positions = calculate_peptide_positions(n_peptides, box_size, 
                                          config['peptide_insertion']['minimum_spacing'])
    z_position = membrane_top + config['peptide_insertion']['distance_from_membrane']
    
    # Place peptides
    all_atoms = []
    atom_counter = 0
    
    # Add peptides first (to match topology order)
    for i, (x, y) in enumerate(positions):
        print(f"  Placing peptide {i+1} at ({x:.1f}, {y:.1f}, {z_position:.1f})")
        
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
                
                # Calculate new coordinates (center peptide and move to position)
                new_x = orig_x - peptide_center[0] + x
                new_y = orig_y - peptide_center[1] + y
                new_z = orig_z - peptide_center[2] + z_position
                
                # Write new atom line
                new_line = f"{resid:5d}{resname}{atomname}{atomid:5d}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}\n"
                all_atoms.append(new_line)
                atom_counter += 1
    
    # Add membrane atoms
    for line in membrane_atoms:
        if len(line) > 44:
            # Update atom number
            parts = list(line)
            atomid = atom_counter + 1
            atomid_str = f"{atomid:5d}"
            line = line[:15] + atomid_str + line[20:]
            all_atoms.append(line)
            atom_counter += 1
    
    # Write combined GRO file
    header = f"System with {n_peptides} {metadata['peptide_id']} peptides; {occupancy} occupancy\n"
    write_gro_file(output_gro, header, all_atoms, box_line)
    
    print(f"✓ Created system with {n_peptides} peptides: {output_gro}")
    
    # Create topology
    create_system_topology(run_dir, metadata['peptide_id'], n_peptides, occupancy)
    
    return output_gro, output_top

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

def create_system_topology(run_dir, peptide_id, n_peptides, occupancy):
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
    
    # Create topology
    top_content = f"""; System topology for {peptide_id} with {occupancy} occupancy
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
{peptide_id}_0    {n_peptides}
"""
    
    # Add membrane molecules
    for mol_name, mol_count in membrane_molecules:
        top_content += f"{mol_name:8s} {mol_count:6d}\n"
    
    # Write topology
    output_top = os.path.join(run_dir, "system", f"system_{occupancy}.top")
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
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Insert peptides
    system_gro, system_top = insert_peptides(args.run_dir, args.occupancy)
    
    print(f"\nNext step: Equilibrate system")
    print(f"Command: python 06_equilibrate.py {args.run_dir} --occupancy {args.occupancy}")

if __name__ == "__main__":
    main()
