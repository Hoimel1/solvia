#!/usr/bin/env python3
"""
SOLVIA Coarse-Graining with Martinize2
Converts atomistic structure to Martini 3 coarse-grained representation
"""

import os
import sys
import yaml
import argparse
import subprocess
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
    """Get best structure from ColabFold"""
    selection_file = os.path.join(run_dir, "colabfold", "model_selection.yaml")
    if not os.path.exists(selection_file):
        print("Error: ColabFold model selection not found. Run 02_run_colabfold.py first.")
        sys.exit(1)
    
    with open(selection_file, 'r') as f:
        selection = yaml.safe_load(f)
    
    pdb_path = os.path.join(run_dir, "colabfold", selection['best_pdb'])
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found: {pdb_path}")
        sys.exit(1)
    
    return pdb_path, selection['best_plddt']

def run_martinize2(run_dir):
    """Run Martinize2 for coarse-graining"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Get input structure
    pdb_path, plddt = get_best_structure(run_dir)
    
    # Output paths
    output_dir = os.path.join(run_dir, "cg_pdb")
    output_top = os.path.join(output_dir, f"{metadata['peptide_id']}.itp")
    output_pdb = os.path.join(output_dir, f"{metadata['peptide_id']}_cg.pdb")
    
    # Check if already done
    if os.path.exists(output_top) and os.path.exists(output_pdb):
        print("Coarse-graining already completed.")
        return output_pdb, output_top
    
    # Build Martinize2 command
    cmd = [
        "martinize2",
        "-f", pdb_path,
        "-o", output_top,
        "-x", output_pdb,
        "-ff", config['coarse_graining']['force_field'],
        "-name", metadata['peptide_id'],
        "-maxwarn", "10",
        "-p", "backbone",
        "-pf", str(config['coarse_graining']['restraint_force']),
    ]
    
    # Add C-terminal modification (always amidated for AMPs)
    if config['coarse_graining']['c_terminal']:
        cmd.extend(["-cter", config['coarse_graining']['c_terminal']])
    
    # Use DSSP for secondary structure
    if config['coarse_graining']['dssp']:
        cmd.append("-dssp")
    
    # For low pLDDT structures, use IDP flag
    if plddt < 70:
        cmd.append("--martini3-idp")
        print(f"Note: Using IDP flag due to low pLDDT ({plddt:.1f})")
    
    # Log file
    log_file = os.path.join(run_dir, "logs", "martinize2.log")
    
    print(f"Running Martinize2 for {metadata['peptide_id']}...")
    print(f"Command: {' '.join(cmd)}")
    
    # Run Martinize2
    with open(log_file, 'w') as log:
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True
            )
            log.write(result.stdout)
            print("✓ Martinize2 completed successfully")
        except subprocess.CalledProcessError as e:
            log.write(e.stdout if e.stdout else "")
            print(f"✗ Martinize2 failed. Check log: {log_file}")
            sys.exit(1)
    
    # Extract secondary structure from log
    extract_secondary_structure(log_file, output_dir, metadata['peptide_id'])
    
    # Fix the ITP file to remove martini.itp include
    fix_itp_file(output_top)
    
    # Move molecule ITP file if created in wrong directory
    molecule_itp = f"{metadata['peptide_id']}_0.itp"
    wrong_path = os.path.join(os.getcwd(), molecule_itp)
    correct_path = os.path.join(output_dir, molecule_itp)
    if os.path.exists(wrong_path) and not os.path.exists(correct_path):
        import shutil
        shutil.move(wrong_path, correct_path)
        print(f"✓ Moved {molecule_itp} to correct directory")
    
    # Create a complete topology file that includes the peptide
    create_complete_topology(output_dir, metadata['peptide_id'], config)
    
    return output_pdb, output_top

def extract_secondary_structure(log_file, output_dir, peptide_id):
    """Extract secondary structure from Martinize2 log"""
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    ss_info = {}
    for i, line in enumerate(lines):
        if "The following sequence of secondary structure" in line:
            # Secondary structure is usually 2 lines after
            if i + 2 < len(lines):
                ss_line = lines[i + 2].strip()
                if ss_line and not ss_line.startswith(';'):
                    ss_info['secondary_structure'] = ss_line
                    print(f"Secondary structure: {ss_line}")
    
    # Save secondary structure info
    if ss_info:
        with open(os.path.join(output_dir, "secondary_structure.yaml"), 'w') as f:
            yaml.dump(ss_info, f, default_flow_style=False)

def fix_itp_file(itp_path):
    """Fix ITP file to remove martini.itp include"""
    # Read the ITP file
    with open(itp_path, 'r') as f:
        lines = f.readlines()
    
    # Replace martini.itp include with a comment
    with open(itp_path, 'w') as f:
        for line in lines:
            if '#include "martini.itp"' in line:
                f.write('; Martini force field is included in the main topology file\n')
            else:
                f.write(line)
    
    print("✓ Fixed ITP file includes")

def create_complete_topology(output_dir, peptide_id, config):
    """Create a complete topology file"""
    top_content = f"""; Complete topology for {peptide_id}
; Generated by SOLVIA coarse-graining script

; Include Martini 3 force field
#include "{config['directories']['force_fields']}/martini_v3.0.0.itp"

; Include peptide topology
#include "{peptide_id}.itp"

[ system ]
{peptide_id} in solution

[ molecules ]
{peptide_id}    1
"""
    
    with open(os.path.join(output_dir, f"{peptide_id}_complete.top"), 'w') as f:
        f.write(top_content)
    
    print(f"✓ Created complete topology: {peptide_id}_complete.top")

def main():
    parser = argparse.ArgumentParser(
        description="Run Martinize2 coarse-graining for SOLVIA"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory"
    )
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Run coarse-graining
    cg_pdb, cg_top = run_martinize2(args.run_dir)
    
    print(f"\n✓ Coarse-grained structure: {os.path.basename(cg_pdb)}")
    print(f"✓ Topology file: {os.path.basename(cg_top)}")
    
    print(f"\nNext step: Build RBC membrane template")
    print(f"Command: python 04_build_membrane.py {args.run_dir}")

if __name__ == "__main__":
    main()
