#!/usr/bin/env python3
"""
Batch processing script for ColabFold structure predictions.
Processes all FASTA files and generates PDB structures.
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
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/colabfold_batch.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class ColabFoldBatch:
    def __init__(self, input_dir, output_dir, colabfold_path, num_models=5, 
                 num_recycle=3, batch_size=10, max_workers=1):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.colabfold_path = Path(colabfold_path)
        self.num_models = num_models
        self.num_recycle = num_recycle
        self.batch_size = batch_size
        self.max_workers = max_workers
        
        # Create output directories
        self.pdb_dir = self.output_dir / "pdb"
        self.raw_output_dir = self.output_dir / "colabfold_raw"
        self.pdb_dir.mkdir(parents=True, exist_ok=True)
        self.raw_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Track progress
        self.progress_file = self.output_dir / "progress.json"
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
                'processing': []
            }
    
    def save_progress(self):
        """Save processing progress to file."""
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def get_fasta_files(self):
        """Get list of FASTA files to process."""
        fasta_files = sorted(self.input_dir.glob("*.fasta"))
        
        # Filter out already processed files
        remaining = []
        for fasta in fasta_files:
            seq_id = fasta.stem
            if seq_id not in self.progress['completed']:
                remaining.append(fasta)
        
        return remaining
    
    def run_colabfold(self, fasta_path):
        """Run ColabFold on a single FASTA file."""
        seq_id = fasta_path.stem
        output_path = self.raw_output_dir / seq_id
        
        # Set up environment
        env = os.environ.copy()
        env['PATH'] = f"{self.colabfold_path}/bin:{env['PATH']}"
        
        # ColabFold command
        cmd = [
            'colabfold_batch',
            '--num-models', str(self.num_models),
            '--num-recycle', str(self.num_recycle),
            str(fasta_path),
            str(output_path)
        ]
        
        logger.info(f"Processing {seq_id}...")
        start_time = time.time()
        
        try:
            # Run ColabFold
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                env=env,
                timeout=600  # 10 minute timeout per sequence
            )
            
            if result.returncode != 0:
                logger.error(f"ColabFold failed for {seq_id}")
                logger.error(f"STDERR: {result.stderr}")
                self.progress['failed'].append({
                    'seq_id': seq_id,
                    'error': result.stderr,
                    'timestamp': datetime.now().isoformat()
                })
                return False
            
            # Find best PDB file and copy to final location
            pdb_files = list(output_path.glob("*_rank_001_*.pdb"))
            if pdb_files:
                best_pdb = pdb_files[0]
                final_pdb = self.pdb_dir / f"{seq_id}.pdb"
                shutil.copy2(best_pdb, final_pdb)
                
                # Extract pLDDT score from JSON
                json_files = list(output_path.glob("*_scores_rank_001_*.json"))
                if json_files:
                    with open(json_files[0], 'r') as f:
                        scores = json.load(f)
                        plddt_list = scores.get('plddt', [])
                        # Calculate average pLDDT if it's a list
                        if isinstance(plddt_list, list) and plddt_list:
                            plddt = sum(plddt_list) / len(plddt_list)
                        else:
                            plddt = plddt_list if isinstance(plddt_list, (int, float)) else 0
                else:
                    plddt = 0
                
                elapsed = time.time() - start_time
                logger.info(f"Completed {seq_id} in {elapsed:.1f}s - pLDDT: {plddt:.1f}")
                
                self.progress['completed'].append({
                    'seq_id': seq_id,
                    'plddt': plddt,
                    'time': elapsed,
                    'timestamp': datetime.now().isoformat()
                })
                return True
            else:
                logger.error(f"No PDB files generated for {seq_id}")
                self.progress['failed'].append({
                    'seq_id': seq_id,
                    'error': 'No PDB files generated',
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
            import traceback
            logger.error(traceback.format_exc())
            self.progress['failed'].append({
                'seq_id': seq_id,
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            })
            return False
    
    def process_batch(self, batch):
        """Process a batch of FASTA files."""
        results = []
        for fasta in batch:
            success = self.run_colabfold(fasta)
            results.append(success)
            self.save_progress()
            
            # Small delay to avoid overwhelming the MSA server
            time.sleep(2)
        
        return results
    
    def run(self):
        """Run batch processing."""
        fasta_files = self.get_fasta_files()
        total = len(fasta_files)
        completed = len(self.progress['completed'])
        
        logger.info(f"Found {total} FASTA files to process")
        logger.info(f"Already completed: {completed}")
        logger.info(f"Remaining: {total}")
        
        if total == 0:
            logger.info("All files already processed!")
            return
        
        # Process in batches
        for i in range(0, total, self.batch_size):
            batch = fasta_files[i:i+self.batch_size]
            batch_num = i // self.batch_size + 1
            total_batches = (total + self.batch_size - 1) // self.batch_size
            
            logger.info(f"\nProcessing batch {batch_num}/{total_batches}")
            self.process_batch(batch)
            
            # Longer delay between batches to respect MSA server limits
            if i + self.batch_size < total:
                logger.info("Waiting 30s before next batch...")
                time.sleep(30)
        
        # Final summary
        logger.info("\n" + "="*50)
        logger.info("PROCESSING COMPLETE")
        logger.info(f"Total processed: {len(self.progress['completed'])}")
        logger.info(f"Failed: {len(self.progress['failed'])}")
        
        if self.progress['failed']:
            logger.info("\nFailed sequences:")
            for fail in self.progress['failed']:
                logger.info(f"  - {fail['seq_id']}: {fail['error']}")

def main():
    parser = argparse.ArgumentParser(description='Batch ColabFold processing')
    parser.add_argument('--input-dir', default='data/raw/fasta_split',
                        help='Directory containing FASTA files')
    parser.add_argument('--output-dir', default='data/processed',
                        help='Output directory')
    parser.add_argument('--colabfold-path', 
                        default='localcolabfold/localcolabfold/colabfold-conda',
                        help='Path to ColabFold conda environment')
    parser.add_argument('--num-models', type=int, default=5,
                        help='Number of models to use')
    parser.add_argument('--num-recycle', type=int, default=3,
                        help='Number of recycling iterations')
    parser.add_argument('--batch-size', type=int, default=10,
                        help='Number of sequences per batch')
    parser.add_argument('--test-run', action='store_true',
                        help='Process only first 5 sequences for testing')
    
    args = parser.parse_args()
    
    # Create logs directory
    Path('logs').mkdir(exist_ok=True)
    
    # Initialize processor
    processor = ColabFoldBatch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        colabfold_path=args.colabfold_path,
        num_models=args.num_models,
        num_recycle=args.num_recycle,
        batch_size=args.batch_size
    )
    
    if args.test_run:
        # For testing, only process first 5 files
        processor.input_dir = Path(args.input_dir)
        fasta_files = sorted(processor.input_dir.glob("*.fasta"))[:5]
        
        # Create temporary directory with only test files
        test_dir = Path("test_fasta")
        test_dir.mkdir(exist_ok=True)
        
        for fasta in fasta_files:
            shutil.copy2(fasta, test_dir)
        
        processor.input_dir = test_dir
        processor.run()
        
        # Clean up
        shutil.rmtree(test_dir)
    else:
        processor.run()

if __name__ == "__main__":
    main()
