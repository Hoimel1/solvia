#!/usr/bin/env python3
"""
Batch processing script for Martinize2 coarse-graining.
Processes all PDB files and generates CG-PDB structures and topologies.
"""

import os
import sys
import subprocess
import time
import json
import argparse
from pathlib import Path
from datetime import datetime
import logging
import shutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/martinize2_batch.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class Martinize2Batch:
    def __init__(self, input_dir, output_dir, force_field='martini3001', 
                 secondary_structure='auto'):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.force_field = force_field
        self.secondary_structure = secondary_structure
        
        # Create output directories
        self.cg_pdb_dir = self.output_dir / "cg_pdb"
        self.topology_dir = self.output_dir / "topologies"
        self.cg_pdb_dir.mkdir(parents=True, exist_ok=True)
        self.topology_dir.mkdir(parents=True, exist_ok=True)
        
        # Track progress
        self.progress_file = self.output_dir / "martinize2_progress.json"
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
                'stats': {}
            }
    
    def save_progress(self):
        """Save processing progress to file."""
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def get_pdb_files(self):
        """Get list of PDB files to process."""
        pdb_files = sorted(self.input_dir.glob("*.pdb"))
        
        # Filter out already processed files
        remaining = []
        for pdb in pdb_files:
            seq_id = pdb.stem
            if seq_id not in self.progress['completed']:
                remaining.append(pdb)
        
        return remaining
    
    def get_secondary_structure(self, pdb_path):
        """Determine secondary structure for the peptide."""
        if self.secondary_structure == 'auto':
            # Count number of residues
            n_residues = 0
            with open(pdb_path, 'r') as f:
                residues = set()
                for line in f:
                    if line.startswith('ATOM'):
                        resid = line[22:26].strip()
                        residues.add(resid)
                n_residues = len(residues)
            
            # For short peptides, use coil structure
            # For longer ones, let martinize2 predict
            if n_residues <= 30:
                return 'C' * n_residues  # All coil
            else:
                return None  # Let martinize2 predict
        else:
            return self.secondary_structure
    
    def run_martinize2(self, pdb_path):
        """Run martinize2 on a single PDB file."""
        seq_id = pdb_path.stem
        cg_pdb_path = self.cg_pdb_dir / f"{seq_id}_cg.pdb"
        topology_path = self.topology_dir / seq_id / "topol.top"
        topology_path.parent.mkdir(exist_ok=True)
        
        # Determine secondary structure
        ss = self.get_secondary_structure(pdb_path)
        
        # Build martinize2 command
        cmd = [
            'martinize2',
            '-f', str(pdb_path),
            '-o', str(topology_path),
            '-x', str(cg_pdb_path),
            '-ff', self.force_field,
            '-ignh'  # Ignore hydrogens
        ]
        
        if ss:
            cmd.extend(['-ss', ss])
        
        logger.info(f"Processing {seq_id}...")
        start_time = time.time()
        
        try:
            # Run martinize2
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True,
                timeout=60  # 1 minute timeout per sequence
            )
            
            if result.returncode != 0:
                logger.error(f"Martinize2 failed for {seq_id}")
                logger.error(f"STDERR: {result.stderr}")
                self.progress['failed'].append({
                    'seq_id': seq_id,
                    'error': result.stderr,
                    'timestamp': datetime.now().isoformat()
                })
                return False
            
            # Check if output files were created
            if cg_pdb_path.exists() and topology_path.exists():
                # Count beads in CG structure
                n_beads = 0
                with open(cg_pdb_path, 'r') as f:
                    for line in f:
                        if line.startswith('ATOM'):
                            n_beads += 1
                
                # Copy molecule ITP file to topology directory
                molecule_itp = Path(topology_path).parent / "molecule_0.itp"
                if molecule_itp.exists():
                    pass  # Already in correct location
                else:
                    # Look in current directory
                    if Path("molecule_0.itp").exists():
                        shutil.move("molecule_0.itp", molecule_itp)
                
                elapsed = time.time() - start_time
                logger.info(f"Completed {seq_id} in {elapsed:.1f}s - {n_beads} beads")
                
                self.progress['completed'].append({
                    'seq_id': seq_id,
                    'n_beads': n_beads,
                    'time': elapsed,
                    'timestamp': datetime.now().isoformat()
                })
                
                # Update stats
                if 'total_beads' not in self.progress['stats']:
                    self.progress['stats']['total_beads'] = 0
                self.progress['stats']['total_beads'] += n_beads
                
                return True
            else:
                logger.error(f"Output files not created for {seq_id}")
                self.progress['failed'].append({
                    'seq_id': seq_id,
                    'error': 'Output files not created',
                    'timestamp': datetime.now().isoformat()
                })
                return False
                
        except subprocess.TimeoutExpired:
            logger.error(f"Timeout for {seq_id}")
            self.progress['failed'].append({
                'seq_id': seq_id,
                'error': 'Timeout',
                'timestamp': datetime.now().isoformat()
            })
            return False
        except Exception as e:
            logger.error(f"Error processing {seq_id}: {type(e).__name__}")
            logger.error(f"Details: {str(e)}")
            self.progress['failed'].append({
                'seq_id': seq_id,
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            })
            return False
    
    def run(self):
        """Run batch processing."""
        pdb_files = self.get_pdb_files()
        total = len(pdb_files)
        completed = len(self.progress['completed'])
        
        logger.info(f"Found {total} PDB files to process")
        logger.info(f"Already completed: {completed}")
        logger.info(f"Force field: {self.force_field}")
        logger.info(f"Secondary structure: {self.secondary_structure}")
        
        if total == 0:
            logger.info("All files already processed!")
            return
        
        # Process files
        for i, pdb in enumerate(pdb_files):
            logger.info(f"\nProcessing {i+1}/{total}")
            success = self.run_martinize2(pdb)
            self.save_progress()
        
        # Final summary
        logger.info("\n" + "="*50)
        logger.info("PROCESSING COMPLETE")
        logger.info(f"Total processed: {len(self.progress['completed'])}")
        logger.info(f"Failed: {len(self.progress['failed'])}")
        
        if self.progress['stats'].get('total_beads'):
            avg_beads = self.progress['stats']['total_beads'] / len(self.progress['completed'])
            logger.info(f"Average beads per peptide: {avg_beads:.1f}")
        
        if self.progress['failed']:
            logger.info("\nFailed sequences:")
            for fail in self.progress['failed']:
                error_first_line = fail['error'].split('\n')[0] if fail['error'] else 'Unknown error'
                logger.info(f"  - {fail['seq_id']}: {error_first_line}")

def main():
    parser = argparse.ArgumentParser(description='Batch martinize2 processing')
    parser.add_argument('--input-dir', default='data/processed/pdb',
                        help='Directory containing PDB files')
    parser.add_argument('--output-dir', default='data/processed',
                        help='Output directory')
    parser.add_argument('--force-field', default='martini3001',
                        help='Martini force field version')
    parser.add_argument('--secondary-structure', default='auto',
                        help='Secondary structure (auto, or specific string)')
    parser.add_argument('--test-run', action='store_true',
                        help='Process only first 5 sequences for testing')
    
    args = parser.parse_args()
    
    # Create logs directory
    Path('logs').mkdir(exist_ok=True)
    
    # Initialize processor
    processor = Martinize2Batch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        force_field=args.force_field,
        secondary_structure=args.secondary_structure
    )
    
    if args.test_run:
        # For testing, only process available PDB files
        pdb_files = sorted(processor.input_dir.glob("*.pdb"))[:5]
        
        if not pdb_files:
            logger.error("No PDB files found for testing!")
            return
        
        # Create temporary list
        processor.get_pdb_files = lambda: pdb_files[:5]
    
    processor.run()

if __name__ == "__main__":
    main()
