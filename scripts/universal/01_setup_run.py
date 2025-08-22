#!/usr/bin/env python3
"""
SOLVIA Setup Script
Creates directory structure for a new peptide run
"""

import os
import sys
import yaml
import shutil
import argparse
from datetime import datetime
from pathlib import Path

def load_config():
    """Load SOLVIA configuration"""
    config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def parse_fasta_header(fasta_file):
    """Extract peptide ID from FASTA file"""
    with open(fasta_file, 'r') as f:
        header = f.readline().strip()
        if not header.startswith('>'):
            raise ValueError("Invalid FASTA file format")
        # Extract SOLVIA ID (e.g., SOLVIA_1)
        peptide_id = header.split('|')[0].replace('>', '').strip()
        return peptide_id

def get_next_run_number(base_dir, peptide_id):
    """Get the next available run number for a peptide"""
    runs = []
    if os.path.exists(base_dir):
        for d in os.listdir(base_dir):
            if d.startswith(f"{peptide_id.lower()}_run_"):
                try:
                    run_num = int(d.split('_')[-1])
                    runs.append(run_num)
                except:
                    pass
    return max(runs) + 1 if runs else 1

def create_directory_structure(run_dir):
    """Create the standard directory structure for a run"""
    directories = [
        "input",
        "colabfold",
        "cg_pdb",
        "membrane_template",
        "system",
        "equilibration/em",
        "equilibration/nvt", 
        "equilibration/npt",
        "production",
        "analysis",
        "logs"
    ]
    
    for d in directories:
        os.makedirs(os.path.join(run_dir, d), exist_ok=True)
        
    # Create README
    readme_content = f"""# {os.path.basename(run_dir)}

Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Directory Structure:
- `input/`: Input FASTA file
- `colabfold/`: ColabFold structure predictions
- `cg_pdb/`: Coarse-grained structures
- `membrane_template/`: RBC membrane template
- `system/`: Combined peptide-membrane system
- `equilibration/`: Energy minimization and equilibration
- `production/`: Production MD simulation
- `analysis/`: Analysis results and ML predictions
- `logs/`: Process logs

## Workflow:
1. Structure prediction (ColabFold)
2. Coarse-graining (Martinize2)
3. Membrane building (INSANE)
4. Peptide insertion
5. Equilibration (EM, NVT, NPT)
6. Production simulation (200ns)
7. Analysis and toxicity prediction
"""
    
    with open(os.path.join(run_dir, "README.md"), 'w') as f:
        f.write(readme_content)

def setup_run(fasta_file, force=False):
    """Main setup function"""
    # Load configuration
    config = load_config()
    
    # Parse peptide ID
    peptide_id = parse_fasta_header(fasta_file)
    
    # Get run number
    base_dir = config['directories']['simulations_base']
    run_number = get_next_run_number(base_dir, peptide_id)
    
    # Create run directory
    run_name = f"{peptide_id.lower()}_run_{run_number}"
    run_dir = os.path.join(base_dir, run_name)
    
    if os.path.exists(run_dir) and not force:
        print(f"Error: Directory {run_dir} already exists. Use --force to overwrite.")
        sys.exit(1)
    
    # Create directory structure
    print(f"Setting up run: {run_name}")
    create_directory_structure(run_dir)
    
    # Copy input FASTA
    input_dir = os.path.join(run_dir, "input")
    shutil.copy(fasta_file, os.path.join(input_dir, "peptide.fasta"))
    
    # Save run metadata
    metadata = {
        'peptide_id': peptide_id,
        'run_number': run_number,
        'run_name': run_name,
        'created': datetime.now().isoformat(),
        'fasta_file': os.path.basename(fasta_file),
        'config_version': config['project']['version']
    }
    
    with open(os.path.join(run_dir, "metadata.yaml"), 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)
    
    print(f"✓ Created directory structure at: {run_dir}")
    print(f"✓ Copied input FASTA to: {input_dir}/peptide.fasta")
    print(f"✓ Run metadata saved")
    print(f"\nNext step: Run ColabFold structure prediction")
    print(f"Command: python 02_run_colabfold.py {run_dir}")
    
    return run_dir

def main():
    parser = argparse.ArgumentParser(
        description="Setup SOLVIA run for peptide hemolytic toxicity prediction"
    )
    parser.add_argument(
        "fasta_file",
        help="Input FASTA file containing peptide sequence"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force overwrite if directory exists"
    )
    
    args = parser.parse_args()
    
    # Check if FASTA file exists
    if not os.path.exists(args.fasta_file):
        print(f"Error: FASTA file not found: {args.fasta_file}")
        sys.exit(1)
    
    # Setup run
    run_dir = setup_run(args.fasta_file, args.force)

if __name__ == "__main__":
    main()
