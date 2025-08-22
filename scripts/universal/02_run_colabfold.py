#!/usr/bin/env python3
"""
SOLVIA ColabFold Structure Prediction
Runs ColabFold for peptide structure prediction
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

def run_colabfold(run_dir):
    """Run ColabFold structure prediction"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Paths
    input_fasta = os.path.join(run_dir, "input", "peptide.fasta")
    output_dir = os.path.join(run_dir, "colabfold")
    
    # Check if already run
    if os.path.exists(os.path.join(output_dir, "ranking_debug.json")):
        print("ColabFold already completed. Checking results...")
        return select_best_model(output_dir, config)
    
    # Build ColabFold command
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
    log_file = os.path.join(run_dir, "logs", "colabfold.log")
    
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

def select_best_model(output_dir, config):
    """Select best model based on pLDDT score"""
    # Find all score files
    import glob
    score_files = glob.glob(os.path.join(output_dir, "*_scores_rank_*.json"))
    
    if not score_files:
        print("Error: No score files found")
        sys.exit(1)
    
    # Extract pLDDT scores from individual score files
    plddt_scores = {}
    for score_file in score_files:
        with open(score_file, 'r') as f:
            data = json.load(f)
        
        # Extract model name from filename
        basename = os.path.basename(score_file)
        # Extract rank and model info
        parts = basename.split('_')
        rank_idx = parts.index('rank')
        rank = parts[rank_idx + 1]
        model_info = '_'.join(parts[rank_idx+2:]).replace('_scores_', '').replace('.json', '')
        
        # Use rank as key for sorting
        key = f"rank_{rank}_{model_info}"
        # plddt is a list of values per residue, calculate mean
        plddt_values = data.get('plddt', [])
        if plddt_values:
            plddt_scores[key] = np.mean(plddt_values)
        else:
            plddt_scores[key] = 0
    
    # Sort by pLDDT
    sorted_models = sorted(plddt_scores.items(), key=lambda x: x[1], reverse=True)
    
    # Check minimum pLDDT
    best_model, best_plddt = sorted_models[0]
    
    print(f"\nModel ranking by pLDDT:")
    for i, (model, plddt) in enumerate(sorted_models[:5]):
        print(f"  {i+1}. {model}: {plddt:.1f}")
    
    if best_plddt < config['colabfold']['min_plddt']:
        print(f"\n✗ Best pLDDT ({best_plddt:.1f}) is below minimum threshold ({config['colabfold']['min_plddt']})")
        print("Structure quality too low for reliable predictions.")
        sys.exit(1)
    
    # Find corresponding PDB file
    # Extract rank number from best model key
    rank_num = best_model.split('_')[1]
    
    # Find PDB file with matching rank
    best_pdb = None
    for pdb_file in os.listdir(output_dir):
        if pdb_file.endswith('.pdb') and f'rank_{rank_num}' in pdb_file and 'unrelaxed' in pdb_file:
            best_pdb = os.path.join(output_dir, pdb_file)
            break
    
    if not best_pdb:
        print(f"Error: Could not find PDB file for {best_model}")
        sys.exit(1)
    
    # Save selection info
    selection = {
        'best_model': best_model,
        'best_plddt': float(best_plddt),
        'best_pdb': os.path.basename(best_pdb),
        'all_scores': {k: float(v) for k, v in plddt_scores.items()}
    }
    
    with open(os.path.join(output_dir, "model_selection.yaml"), 'w') as f:
        yaml.dump(selection, f, default_flow_style=False)
    
    print(f"\n✓ Selected best model: {best_model} (pLDDT: {best_plddt:.1f})")
    print(f"✓ PDB file: {os.path.basename(best_pdb)}")
    
    return best_pdb, best_plddt

def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold structure prediction for SOLVIA"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory created by setup script"
    )
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Check if metadata exists
    if not os.path.exists(os.path.join(args.run_dir, "metadata.yaml")):
        print(f"Error: metadata.yaml not found in {args.run_dir}")
        print("Make sure to run 01_setup_run.py first")
        sys.exit(1)
    
    # Run ColabFold
    best_pdb, best_plddt = run_colabfold(args.run_dir)
    
    print(f"\nNext step: Coarse-graining with Martinize2")
    print(f"Command: python 03_coarse_grain.py {args.run_dir}")

if __name__ == "__main__":
    main()
