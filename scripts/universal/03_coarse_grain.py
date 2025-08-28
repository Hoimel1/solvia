#!/usr/bin/env python3
"""
SOLVIA Peptide Coarse-Graining
Converts atomistic structures to Martini 3 coarse-grained representation
Handles all ITP processing and symlink creation automatically
"""

import os
import sys
import yaml
import argparse
import subprocess
import shutil
import re
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

def get_best_structure(run_dir):
    """Get best structure from ColabFold or other source"""
    # Check for model selection file
    selection_file = os.path.join(run_dir, "colabfold", "model_selection.yaml")
    if not os.path.exists(selection_file):
        print("Error: No model selection found. Run ColabFold first.")
        sys.exit(1)
    
    with open(selection_file, 'r') as f:
        selection = yaml.safe_load(f)
    
    # Use best_model (filename) or best_pdb (full path)
    if 'best_model' in selection and selection['best_model']:
        pdb_path = os.path.join(run_dir, "colabfold", selection['best_model'])
    elif 'best_pdb' in selection:
        # best_pdb might be a full path or relative path
        if os.path.isabs(selection['best_pdb']):
            pdb_path = selection['best_pdb']
        else:
            pdb_path = selection['best_pdb']
    else:
        print("Error: No best model found in selection file")
        sys.exit(1)
    
    # Create best_model.pdb symlink for martinize
    best_model_link = os.path.join(run_dir, "colabfold", "best_model.pdb")
    if not os.path.exists(best_model_link):
        os.symlink(os.path.basename(pdb_path), best_model_link)
    
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found: {pdb_path}")
        sys.exit(1)
    
    return pdb_path, selection.get('best_plddt', 80)

def extract_peptide_id(run_dir):
    """Extract peptide ID from run directory name"""
    run_dir_name = os.path.basename(run_dir)  # e.g., "solvia_1_run_1"
    # Extract peptide_id: solvia_1_run_1 -> SOLVIA_1
    parts = run_dir_name.split('_')
    if len(parts) >= 2:
        peptide_id = '_'.join(parts[:2]).upper()  # "SOLVIA_1"
    else:
        peptide_id = "PEPTIDE"
    return peptide_id

def run_martinize_docker(run_dir, pdb_path, config):
    """Run Martinize via Docker and handle all ITP processing"""
    print(f"\n{'='*50}")
    print("Running Martinize via Docker...")
    print(f"{'='*50}")
    
    # Extract configuration
    cg_config = config.get('coarse_graining', {})
    c_terminal = cg_config.get('c_terminal', 'NH2-ter')
    force_field = cg_config.get('force_field', 'martini3001')
    
    # Get peptide ID
    peptide_id = extract_peptide_id(run_dir)
    
    # Get relative path from project root
    project_root = Path(__file__).parent.parent.parent
    rel_run_dir = os.path.relpath(run_dir, project_root)
    
    # Create output directory
    output_dir = os.path.join(run_dir, "coarse_grain")
    os.makedirs(output_dir, exist_ok=True)
    
    # Build Docker command
    cmd = [
        'docker', 'compose', 'run', '--rm', 'martinize',
        '-f', f'/work/{rel_run_dir}/colabfold/best_model.pdb',
        '-ff', force_field,
        '-x', f'/work/{rel_run_dir}/coarse_grain/{peptide_id}_cg.pdb',
        '-o', f'/work/{rel_run_dir}/coarse_grain/{peptide_id}.itp',
        '-name', peptide_id,
        '-cter', c_terminal,  # C-terminal amidation
        '-dssp',               # Use DSSP for secondary structure
        '-elastic',            # Add elastic network
        '-ef', '700',          # Elastic network force constant
        '-eu', '0.9',          # Elastic network upper cutoff
        '-p', 'backbone',      # Position restraints on backbone
        '-pf', '1000',         # Position restraint force constant
        '-maxwarn', '1',
        '-cys', 'auto'         # Auto-detect disulfide bonds
    ]
    
    print(f"Peptide ID: {peptide_id}")
    print(f"Command: {' '.join(cmd)}")
    
    # Run Martinize via Docker
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(project_root))
    
    if result.returncode != 0:
        print("\nError running Martinize:")
        print(result.stderr)
        sys.exit(1)
    
    print("\nMartinize output:")
    print(result.stdout)
    
    # Handle the ITP files created by martinize
    handle_itp_files(project_root, output_dir, peptide_id)
    
    # Create symlinks
    create_symlinks(output_dir, peptide_id)
    
    # Create position restraint file
    create_position_restraints(output_dir, peptide_id)
    
    return peptide_id

def handle_itp_files(project_root, output_dir, peptide_id):
    """Handle the two ITP files that martinize creates"""
    print("\nProcessing ITP files...")
    
    # Martinize creates both NAME.itp (master) and NAME_0.itp (molecule)
    # We only need the molecule ITP to avoid double counting
    
    molecule_itp_src = os.path.join(project_root, f"{peptide_id}_0.itp")
    molecule_itp_dst = os.path.join(output_dir, f"{peptide_id}.itp")
    master_itp = os.path.join(project_root, f"{peptide_id}.itp")
    
    # Check if _0.itp exists and move it
    if os.path.exists(molecule_itp_src):
        # Move the molecule ITP to the correct location
        shutil.move(molecule_itp_src, molecule_itp_dst)
        print(f"  ✓ Moved {peptide_id}_0.itp to coarse_grain/{peptide_id}.itp")
        
        # Remove the master ITP (causes double counting)
        if os.path.exists(master_itp):
            os.remove(master_itp)
            print(f"  ✓ Removed master ITP (prevents double counting)")
    
    # Also check if files were created directly in output_dir (newer martinize versions)
    alt_molecule_itp = os.path.join(output_dir, f"{peptide_id}_0.itp")
    alt_master_itp = os.path.join(output_dir, f"{peptide_id}.itp")
    
    if os.path.exists(alt_molecule_itp) and not os.path.exists(molecule_itp_dst):
        shutil.move(alt_molecule_itp, molecule_itp_dst)
        print(f"  ✓ Renamed {peptide_id}_0.itp to {peptide_id}.itp")
    
    # Fix molecule name if needed (martinize might use wrong name)
    if os.path.exists(molecule_itp_dst):
        with open(molecule_itp_dst, 'r') as f:
            content = f.read()
        
        # Replace any wrong molecule names
        original_content = content
        content = re.sub(r'SIMULATIONS(_0)?', peptide_id, content)
        content = re.sub(f'{peptide_id}_0', peptide_id, content)  # Remove _0 suffix from molecule name
        
        if content != original_content:
            with open(molecule_itp_dst, 'w') as f:
                f.write(content)
            print(f"  ✓ Fixed molecule name to {peptide_id}")

def create_symlinks(output_dir, peptide_id):
    """Create symlinks for backward compatibility"""
    print("\nCreating symlinks...")
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    # Create symlinks to the actual files
    symlink_map = {
        'peptide_cg.pdb': f'{peptide_id}_cg.pdb',
        'peptide_cg.gro': f'{peptide_id}_cg.pdb',  # GROMACS can read PDB as GRO
        'peptide_cg.itp': f'{peptide_id}.itp'
    }
    
    for link, target in symlink_map.items():
        if os.path.exists(link) or os.path.islink(link):
            os.remove(link)
        if os.path.exists(target):
            os.symlink(target, link)
            print(f"  ✓ {link} -> {target}")
    
    # Return to original directory
    os.chdir(original_dir)

def create_position_restraints(output_dir, peptide_id):
    """Create position restraint file"""
    print("\nGenerating position restraints...")
    
    posre_file = os.path.join(output_dir, "peptide_posre.itp")
    itp_file = os.path.join(output_dir, f"{peptide_id}.itp")
    
    if not os.path.exists(itp_file):
        print(f"  Warning: ITP file not found: {itp_file}")
        return
    
    with open(posre_file, 'w') as f:
        f.write(f"; Position restraint file for {peptide_id}\n")
        f.write("; Include this file in topology with:\n")
        f.write("; #ifdef POSRES\n")
        f.write("; #include \"peptide_posre.itp\"\n")
        f.write("; #endif\n\n")
        f.write("[ position_restraints ]\n")
        f.write("; atom  type      fx      fy      fz\n")
        
        # Extract backbone atoms from ITP
        with open(itp_file, 'r') as itp:
            in_atoms = False
            for line in itp:
                if '[ atoms ]' in line:
                    in_atoms = True
                    continue
                elif in_atoms and line.startswith('['):
                    break
                elif in_atoms and 'BB' in line and not line.startswith(';'):
                    parts = line.split()
                    if len(parts) > 0:
                        atom_id = parts[0]
                        f.write(f"{atom_id:>6}     1  1000  1000  1000\n")
    
    print(f"  ✓ Created {posre_file}")

def print_summary(run_dir, peptide_id, config):
    """Print summary of coarse-graining results"""
    rel_run_dir = os.path.relpath(run_dir, Path(__file__).parent.parent.parent)
    c_terminal = config.get('coarse_graining', {}).get('c_terminal', 'NH2-ter')
    
    print(f"\n{'='*50}")
    print("✅ Coarse-graining complete!")
    print(f"{'='*50}")
    print(f"Output files:")
    print(f"  PDB: {rel_run_dir}/coarse_grain/{peptide_id}_cg.pdb")
    print(f"  ITP: {rel_run_dir}/coarse_grain/{peptide_id}.itp (single file, no double counting)")
    print(f"  Position restraints: {rel_run_dir}/coarse_grain/peptide_posre.itp")
    print(f"\nSymlinks created:")
    print(f"  peptide_cg.pdb -> {peptide_id}_cg.pdb")
    print(f"  peptide_cg.itp -> {peptide_id}.itp")
    print(f"  peptide_cg.gro -> {peptide_id}_cg.pdb")
    print(f"\nImportant notes:")
    print(f"  Molecule name: {peptide_id}")
    print(f"  C-terminal: {c_terminal} (amidated)")
    print(f"  N-terminal: Protonated (NH3+)")
    print(f"{'='*50}")

def main():
    parser = argparse.ArgumentParser(
        description="Coarse-grain peptide structure using Martinize via Docker"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory containing ColabFold output"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run even if output exists"
    )
    
    args = parser.parse_args()
    
    # Convert to absolute path
    run_dir = os.path.abspath(args.run_dir)
    
    # Check run directory
    if not os.path.exists(run_dir):
        print(f"Error: Run directory not found: {run_dir}")
        sys.exit(1)
    
    # Check if already done (unless forced)
    peptide_id = extract_peptide_id(run_dir)
    output_dir = os.path.join(run_dir, "coarse_grain")
    output_pdb = os.path.join(output_dir, f"{peptide_id}_cg.pdb")
    output_itp = os.path.join(output_dir, f"{peptide_id}.itp")
    
    if not args.force and os.path.exists(output_pdb) and os.path.exists(output_itp):
        print("Coarse-graining already completed. Use --force to re-run.")
        print_summary(run_dir, peptide_id, load_config())
        return
    
    # Load configuration
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Get best structure
    pdb_path, plddt = get_best_structure(run_dir)
    print(f"Using best structure: {os.path.basename(pdb_path)} (pLDDT: {plddt:.1f})")
    
    # Run martinize via Docker
    peptide_id = run_martinize_docker(run_dir, pdb_path, config)
    
    # Print summary
    print_summary(run_dir, peptide_id, config)

if __name__ == "__main__":
    main()