#!/usr/bin/env python3
"""
SOLVIA ColabFold Structure Prediction
Runs ColabFold for peptide structure prediction or selects best existing model
"""

import os
import sys
import yaml
import json
import argparse
import subprocess
from pathlib import Path
import numpy as np

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

def check_colabfold_available():
    """Check if colabfold_batch is available in PATH"""
    try:
        result = subprocess.run(
            ["which", "colabfold_batch"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result.returncode == 0
    except:
        return False

def run_colabfold(run_dir):
    """Run ColabFold structure prediction or select best model if already run"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Paths
    input_fasta = os.path.join(run_dir, "input", "peptide.fasta")
    output_dir = os.path.join(run_dir, "colabfold")
    
    # Check if ColabFold outputs already exist
    pdb_files = list(Path(output_dir).glob("*_unrelaxed_rank_*.pdb")) if os.path.exists(output_dir) else []
    relaxed_files = list(Path(output_dir).glob("*_relaxed_rank_*.pdb")) if os.path.exists(output_dir) else []
    
    if pdb_files or relaxed_files:
        print(f"Found {len(pdb_files + relaxed_files)} ColabFold models. Selecting best...")
        return select_best_model(output_dir, config)
    
    # Check if we can run ColabFold natively
    colabfold_available = check_colabfold_available()
    
    if not colabfold_available:
        print("")
        print("="*70)
        print("ColabFold not found in PATH.")
        print("="*70)
        print("")
        print("Please run ColabFold using Docker with this command:")
        print("")
        print(f"docker compose run --rm \\")
        print(f"  -e XDG_CACHE_HOME=/cache \\")
        print(f"  -e MPLCONFIGDIR=/cache/mpl \\")
        print(f"  colabfold \\")
        print(f"  /work/{run_dir}/input/peptide.fasta /work/{run_dir}/colabfold \\")
        print(f"  --num-seeds {config['colabfold']['num_replicates']} \\")
        print(f"  --num-models {config['colabfold']['num_models']} \\")
        print(f"  --msa-mode {config['colabfold']['msa_mode']} \\")
        print(f"  --pair-mode {config['colabfold']['pair_mode']}")
        print("")
        print("After ColabFold completes, run this script again to select the best model.")
        print("="*70)
        return None, 0
    
    # Build ColabFold command (native installation)
    cmd = [
        "colabfold_batch",
        input_fasta,
        output_dir,
        "--num-seeds", str(config['colabfold']['num_replicates']),
        "--num-models", str(config['colabfold']['num_models']),
        "--msa-mode", config['colabfold']['msa_mode'],
        "--pair-mode", config['colabfold']['pair_mode'],
    ]
    
    if config['colabfold']['relax']:
        cmd.append("--amber")
        cmd.append("--use-gpu-relax")
    
    # Log file
    log_dir = os.path.join(run_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "colabfold.log")
    
    print(f"Running ColabFold for {metadata['peptide_id']}...")
    print(f"Command: {' '.join(cmd)}")
    print(f"This may take 10-30 minutes depending on sequence length...")
    
    # Run ColabFold
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
            print("✓ ColabFold completed successfully")
        except subprocess.CalledProcessError as e:
            log.write(e.stdout if e.stdout else "")
            print(f"✗ ColabFold failed. Check log: {log_file}")
            sys.exit(1)
    
    # Select best model
    return select_best_model(output_dir, config)

def select_best_model_by_bfactor(pdb_files):
    """Select best model based on B-factors when score files are missing"""
    best_pdb = None
    best_plddt = 0
    
    for pdb_file in pdb_files:
        # Read PDB and extract B-factors (pLDDT scores)
        plddt_values = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        b_factor = float(line[60:66])
                        plddt_values.append(b_factor)
                    except:
                        pass
        
        if plddt_values:
            avg_plddt = np.mean(plddt_values)
            print(f"  {os.path.basename(pdb_file)}: pLDDT = {avg_plddt:.1f}")
            if avg_plddt > best_plddt:
                best_plddt = avg_plddt
                best_pdb = str(pdb_file)
    
    return best_pdb, best_plddt

def select_best_model(output_dir, config):
    """Select best model based on pLDDT score"""
    # Find all PDB files first
    import glob
    pdb_files = glob.glob(os.path.join(output_dir, "*_unrelaxed_rank_*.pdb"))
    relaxed_pdb_files = glob.glob(os.path.join(output_dir, "*_relaxed_rank_*.pdb"))
    all_pdb_files = (relaxed_pdb_files if relaxed_pdb_files else pdb_files)
    
    if not all_pdb_files:
        print("Error: No PDB files found in", output_dir)
        print("Please run ColabFold first using the Docker command shown above.")
        return None, 0
    
    # Find all score files
    score_files = glob.glob(os.path.join(output_dir, "*_scores_rank_*.json"))
    
    if not score_files:
        print("Warning: No score files found, extracting pLDDT from B-factors...")
        best_pdb, best_plddt = select_best_model_by_bfactor(all_pdb_files)
    else:
        # Extract pLDDT scores from score files
        plddt_scores = {}
        for score_file in score_files:
            with open(score_file, 'r') as f:
                data = json.load(f)
            
            # Extract rank from filename
            import re
            match = re.search(r'rank_(\d+)', os.path.basename(score_file))
            if match:
                rank = int(match.group(1))
                model_name = f"rank_{rank:03d}"
                
                # pLDDT is usually stored as 'plddt' or 'plddts' 
                if 'plddt' in data:
                    if isinstance(data['plddt'], list):
                        plddt = np.mean(data['plddt'])
                    else:
                        plddt = data['plddt']
                elif 'plddts' in data:
                    if isinstance(data['plddts'], list):
                        plddt = np.mean(data['plddts'])
                    else:
                        plddt = data['plddts']
                else:
                    continue
                
                plddt_scores[model_name] = plddt
                print(f"  Model {model_name}: pLDDT = {plddt:.1f}")
        
        if not plddt_scores:
            print("Warning: Could not extract pLDDT from score files, using B-factors...")
            best_pdb, best_plddt = select_best_model_by_bfactor(all_pdb_files)
        else:
            # Select best model
            best_model = max(plddt_scores, key=plddt_scores.get)
            best_plddt = plddt_scores[best_model]
            
            # Find corresponding PDB file
            rank_num = int(best_model.split('_')[1])
            best_pdb = None
            
            for pdb_file in all_pdb_files:
                if f"rank_{rank_num:03d}" in pdb_file or f"rank_{rank_num:01d}" in pdb_file:
                    best_pdb = pdb_file
                    break
            
            if not best_pdb:
                print(f"Warning: Could not find PDB file for {best_model}")
                best_pdb = all_pdb_files[0]  # Fallback to first file
    
    # Check minimum pLDDT threshold
    if best_plddt < config['colabfold']['min_plddt']:
        print(f"⚠ Warning: Best pLDDT ({best_plddt:.1f}) is below threshold ({config['colabfold']['min_plddt']})")
        print("Consider using a longer sequence or different modeling approach")
    
    # Save selection metadata
    selection = {
        'best_model': os.path.basename(best_pdb) if best_pdb else None,
        'best_plddt': float(best_plddt),
        'best_pdb': best_pdb,
        'min_plddt_threshold': config['colabfold']['min_plddt']
    }
    
    selection_file = os.path.join(output_dir, "model_selection.yaml")
    with open(selection_file, 'w') as f:
        yaml.dump(selection, f, default_flow_style=False)
    
    # Create symbolic link to best model
    if best_pdb and os.path.exists(best_pdb):
        best_link = os.path.join(output_dir, "best_model.pdb")
        if os.path.exists(best_link):
            os.remove(best_link)
        os.symlink(os.path.basename(best_pdb), best_link)
    
    print(f"\n✓ Model selection complete")
    print(f"  Best model: {os.path.basename(best_pdb) if best_pdb else 'None'}")
    print(f"  pLDDT: {best_plddt:.1f}")
    print(f"  Selection saved: {selection_file}")
    
    return best_pdb, best_plddt

def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold structure prediction or select best model"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory containing input FASTA"
    )
    
    args = parser.parse_args()
    
    # Check run directory
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Run ColabFold or select model
    best_pdb, best_plddt = run_colabfold(args.run_dir)
    
    if best_pdb:
        print(f"\n✓ Ready for next step: Coarse-graining")
        print(f"  Command: python3 scripts/universal/03_coarse_grain.py {args.run_dir}")

if __name__ == "__main__":
    main()
