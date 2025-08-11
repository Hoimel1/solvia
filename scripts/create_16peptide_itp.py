#!/usr/bin/env python3
"""
Create a single ITP file for 16 SOLVIA_1 peptides combined.
This is needed because insane treats the 16x array as one molecule.
"""

import sys
from pathlib import Path

def create_combined_itp(input_itp, output_itp, n_copies=16):
    """Create ITP for 16 peptides as one molecule."""
    
    # Read original ITP
    with open(input_itp, 'r') as f:
        lines = f.readlines()
    
    # Find sections
    atoms_start = -1
    atoms_end = -1
    bonds_start = -1
    bonds_end = -1
    
    for i, line in enumerate(lines):
        if '[ atoms ]' in line:
            atoms_start = i + 1
        elif atoms_start > 0 and atoms_end < 0 and '[' in line:
            atoms_end = i
        elif '[ bonds ]' in line:
            bonds_start = i + 1
        elif bonds_start > 0 and bonds_end < 0 and '[' in line:
            bonds_end = i
    
    # Handle end of file
    if atoms_end < 0 and atoms_start > 0:
        atoms_end = len(lines)
    if bonds_end < 0 and bonds_start > 0:
        bonds_end = len(lines)
    
    # Extract atoms and bonds
    atom_lines = lines[atoms_start:atoms_end]
    bond_lines = lines[bonds_start:bonds_end] if bonds_start > 0 else []
    
    # Count atoms per peptide
    n_atoms = len([l for l in atom_lines if l.strip() and not l.startswith(';')])
    
    # Create combined ITP
    output_lines = []
    
    # Header
    output_lines.append("; Combined ITP for 16x SOLVIA_1\n")
    output_lines.append("; Generated from single peptide ITP\n")
    output_lines.append("\n")
    output_lines.append("[ moleculetype ]\n")
    output_lines.append("; Name    nrexcl\n") 
    output_lines.append("SOLVIA_1    1\n")
    output_lines.append("\n")
    
    # Atoms section
    output_lines.append("[ atoms ]\n")
    output_lines.append(";   nr  type  resnr residue  atom   cgnr     charge       mass\n")
    
    atom_counter = 0
    for copy in range(n_copies):
        res_offset = copy * 100  # Offset residue numbers
        
        for line in atom_lines:
            if line.strip() and not line.startswith(';') and not line.startswith('['):
                # Handle the special format where atom number might be separate
                line = line.strip()
                parts = line.split()
                if len(parts) >= 6:
                    atom_counter += 1
                    atom_idx = atom_counter
                    atom_type = parts[1]
                    res_num = int(parts[2]) + res_offset
                    res_name = parts[3]
                    atom_name = parts[4]
                    cgnr = int(parts[5]) + copy * n_atoms
                    # Handle charge and mass which might be missing
                    charge = parts[6] if len(parts) > 6 else "0.0"
                    mass = parts[7] if len(parts) > 7 else ""
                    
                    # Reconstruct the line preserving format
                    if mass:
                        output_lines.append(f"{atom_idx:2d} {atom_type:<5s} {res_num:2d} {res_name} {atom_name:<4s} {cgnr:2d} {charge:>4s} {mass:>5s}\n")
                    else:
                        output_lines.append(f"{atom_idx:2d} {atom_type:<5s} {res_num:2d} {res_name} {atom_name:<4s} {cgnr:2d} {charge:>6s}\n")
    
    output_lines.append("\n")
    
    # Bonds section - only within each peptide
    if bond_lines:
        output_lines.append("[ bonds ]\n")
        output_lines.append(";  ai    aj funct            c0            c1            c2            c3\n")
        
        for copy in range(n_copies):
            atom_offset = copy * n_atoms
            
            for line in bond_lines:
                if line.strip() and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) >= 2:
                        ai = int(parts[0]) + atom_offset
                        aj = int(parts[1]) + atom_offset
                        rest = ' '.join(parts[2:]) if len(parts) > 2 else ""
                        output_lines.append(f"{ai:6d} {aj:6d} {rest}\n")
    
    # Write output
    with open(output_itp, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Created combined ITP with {atom_counter} atoms")
    print(f"Output: {output_itp}")

if __name__ == "__main__":
    input_itp = Path("data/processed/topologies/SOLVIA_1/SOLVIA_1.itp")
    output_itp = Path("simulations/rbc_simple/SOLVIA_1_16x.itp")
    
    create_combined_itp(input_itp, output_itp, n_copies=16)
