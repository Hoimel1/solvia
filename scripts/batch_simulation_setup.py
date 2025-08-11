#!/usr/bin/env python3
"""
Batch setup of membrane simulations for all coarse-grained peptides.
"""

import os
import sys
import subprocess
import json
import argparse
from pathlib import Path
from datetime import datetime
import logging
import multiprocessing as mp
from functools import partial

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/batch_simulation_setup.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class BatchSimulationSetup:
    def __init__(self, cg_pdb_dir, topology_dir, output_dir, force_field_dir, 
                 num_workers=1):
        self.cg_pdb_dir = Path(cg_pdb_dir)
        self.topology_dir = Path(topology_dir)
        self.output_dir = Path(output_dir)
        self.force_field_dir = Path(force_field_dir)
        self.num_workers = num_workers
        
        # Progress tracking
        self.progress_file = self.output_dir / "simulation_setup_progress.json"
        self.load_progress()
        
    def load_progress(self):
        """Load processing progress from file."""
        if self.progress_file.exists():
            with open(self.progress_file, 'r') as f:
                self.progress = json.load(f)
        else:
            self.progress = {
                'completed': [],
                'failed': [],
                'timestamp': datetime.now().isoformat()
            }
    
    def save_progress(self):
        """Save processing progress to file."""
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def get_peptides_to_process(self):
        """Get list of peptides that haven't been processed yet."""
        # Get all CG PDB files
        cg_pdbs = list(self.cg_pdb_dir.glob("*_cg.pdb"))
        
        # Extract peptide IDs
        peptides = []
        for pdb in cg_pdbs:
            peptide_id = pdb.stem.replace('_cg', '')
            
            # Check if already processed
            if peptide_id not in self.progress['completed']:
                # Check if topology exists
                topology_dir = self.topology_dir / peptide_id
                if topology_dir.exists():
                    peptides.append({
                        'id': peptide_id,
                        'cg_pdb': str(pdb),
                        'topology_dir': str(self.topology_dir)
                    })
                else:
                    logger.warning(f"No topology found for {peptide_id}, skipping")
        
        return peptides
    
    def setup_single_peptide(self, peptide_info):
        """Setup simulation for a single peptide."""
        peptide_id = peptide_info['id']
        
        logger.info(f"Setting up simulation for {peptide_id}")
        
        # Build command
        cmd = [
            sys.executable,
            'scripts/setup_membrane_simulation.py',
            '--peptide-id', peptide_id,
            '--cg-pdb', peptide_info['cg_pdb'],
            '--topology-dir', peptide_info['topology_dir'],
            '--output-dir', str(self.output_dir),
            '--force-field-dir', str(self.force_field_dir)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                self.progress['completed'].append({
                    'peptide_id': peptide_id,
                    'timestamp': datetime.now().isoformat()
                })
                logger.info(f"Successfully set up {peptide_id}")
                return True
            else:
                self.progress['failed'].append({
                    'peptide_id': peptide_id,
                    'error': result.stderr,
                    'timestamp': datetime.now().isoformat()
                })
                logger.error(f"Failed to set up {peptide_id}: {result.stderr}")
                return False
                
        except Exception as e:
            self.progress['failed'].append({
                'peptide_id': peptide_id,
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            })
            logger.error(f"Exception setting up {peptide_id}: {e}")
            return False
    
    def run(self):
        """Run batch simulation setup."""
        peptides = self.get_peptides_to_process()
        total = len(peptides)
        
        logger.info(f"Found {total} peptides to process")
        logger.info(f"Already completed: {len(self.progress['completed'])}")
        logger.info(f"Failed: {len(self.progress['failed'])}")
        
        if total == 0:
            logger.info("No new peptides to process")
            return
        
        # Process peptides
        if self.num_workers > 1:
            # Parallel processing
            logger.info(f"Using {self.num_workers} workers")
            with mp.Pool(self.num_workers) as pool:
                results = pool.map(self.setup_single_peptide, peptides)
        else:
            # Sequential processing
            results = []
            for i, peptide in enumerate(peptides):
                logger.info(f"Processing {i+1}/{total}: {peptide['id']}")
                success = self.setup_single_peptide(peptide)
                results.append(success)
                self.save_progress()
        
        # Summary
        successful = sum(results)
        failed = total - successful
        
        logger.info("\n" + "="*50)
        logger.info("BATCH SETUP COMPLETE")
        logger.info(f"Successfully set up: {successful}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Total completed: {len(self.progress['completed'])}")
        
        # Save final progress
        self.save_progress()
        
        # Create summary report
        self.create_summary_report()
    
    def create_summary_report(self):
        """Create a summary report of the setup process."""
        report_path = self.output_dir / "setup_summary.txt"
        
        with open(report_path, 'w') as f:
            f.write("SOLVIA Simulation Setup Summary\n")
            f.write("="*50 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"Total peptides processed: {len(self.progress['completed'])}\n")
            f.write(f"Failed setups: {len(self.progress['failed'])}\n\n")
            
            f.write("Successfully set up peptides:\n")
            for item in self.progress['completed']:
                f.write(f"  - {item['peptide_id']}\n")
            
            if self.progress['failed']:
                f.write("\nFailed peptides:\n")
                for item in self.progress['failed']:
                    f.write(f"  - {item['peptide_id']}: {item['error'][:80]}...\n")
            
            f.write("\nNext steps:\n")
            f.write("1. Upload Martini force field files to force_fields/martini3/\n")
            f.write("2. Complete membrane building when INSANE is available\n")
            f.write("3. Submit simulations using the generated submit.sh scripts\n")
        
        logger.info(f"Summary report saved to {report_path}")

def main():
    parser = argparse.ArgumentParser(description='Batch setup membrane simulations')
    parser.add_argument('--cg-pdb-dir', default='data/processed/cg_pdb',
                        help='Directory containing CG PDB files')
    parser.add_argument('--topology-dir', default='data/processed/topologies',
                        help='Directory containing topology files')
    parser.add_argument('--output-dir', default='simulations/systems',
                        help='Output directory for simulations')
    parser.add_argument('--force-field-dir', default='force_fields/martini3',
                        help='Force field directory')
    parser.add_argument('--num-workers', type=int, default=1,
                        help='Number of parallel workers')
    parser.add_argument('--test-run', action='store_true',
                        help='Process only first 5 peptides for testing')
    
    args = parser.parse_args()
    
    # Create necessary directories
    Path('logs').mkdir(exist_ok=True)
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Initialize batch processor
    processor = BatchSimulationSetup(
        cg_pdb_dir=args.cg_pdb_dir,
        topology_dir=args.topology_dir,
        output_dir=args.output_dir,
        force_field_dir=args.force_field_dir,
        num_workers=args.num_workers
    )
    
    if args.test_run:
        # Limit to first 5 peptides for testing
        peptides = processor.get_peptides_to_process()[:5]
        processor.get_peptides_to_process = lambda: peptides
    
    processor.run()

if __name__ == "__main__":
    main()
