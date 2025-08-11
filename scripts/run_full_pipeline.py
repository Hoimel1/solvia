#!/usr/bin/env python3
"""
Run the complete SOLVIA pipeline:
1. Structure prediction (ColabFold)
2. Coarse-graining (Martinize2)
3. MD simulation setup
4. Feature extraction
5. ML model training

This script coordinates all modules for end-to-end processing.
"""

import os
import sys
import subprocess
import json
import pandas as pd
from pathlib import Path
import logging
import argparse
from datetime import datetime
import shutil

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SolviaPipeline:
    def __init__(self, config_file=None):
        """Initialize pipeline with configuration."""
        self.project_root = Path(__file__).parent.parent
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Default configuration
        self.config = {
            'peptide_list': 'data/raw/test_peptides.txt',
            'fasta_dir': 'data/raw/fasta_split',
            'pdb_dir': 'data/processed/pdb',
            'cg_pdb_dir': 'data/processed/cg_pdb',
            'topology_dir': 'data/processed/topologies',
            'simulation_dir': 'simulations/systems',
            'feature_dir': 'data/features',
            'model_dir': 'models',
            'results_dir': 'results',
            'num_workers': 4,
            'use_gpu': True,
            'test_mode': False
        }
        
        # Load custom config if provided
        if config_file and Path(config_file).exists():
            with open(config_file, 'r') as f:
                custom_config = json.load(f)
                self.config.update(custom_config)
                
        # Create directories
        for key, path in self.config.items():
            if key.endswith('_dir'):
                Path(path).mkdir(parents=True, exist_ok=True)
                
    def run_structure_prediction(self, peptide_ids):
        """Module 1: Run ColabFold for structure prediction."""
        logger.info("="*60)
        logger.info("MODULE 1: Structure Prediction (ColabFold)")
        logger.info("="*60)
        
        # Check if already done
        completed = []
        pending = []
        
        for peptide_id in peptide_ids:
            pdb_file = Path(self.config['pdb_dir']) / f"{peptide_id}.pdb"
            if pdb_file.exists():
                completed.append(peptide_id)
            else:
                pending.append(peptide_id)
                
        logger.info(f"Already completed: {len(completed)}")
        logger.info(f"Pending: {len(pending)}")
        
        if not pending:
            logger.info("All structures already predicted!")
            return True
            
        # Run batch ColabFold
        cmd = [
            "python", "scripts/batch_colabfold.py",
            "--input-dir", self.config['fasta_dir'],
            "--output-dir", self.config['pdb_dir']
        ]
        
        if self.config['test_mode']:
            cmd.extend(["--test-run"])
            
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=self.project_root)
        
        return result.returncode == 0
        
    def run_coarse_graining(self, peptide_ids):
        """Module 2: Run Martinize2 for coarse-graining."""
        logger.info("="*60)
        logger.info("MODULE 2: Coarse-Graining (Martinize2)")
        logger.info("="*60)
        
        # Check completed
        completed = []
        pending = []
        
        for peptide_id in peptide_ids:
            cg_file = Path(self.config['cg_pdb_dir']) / f"{peptide_id}_cg.pdb"
            if cg_file.exists():
                completed.append(peptide_id)
            else:
                pending.append(peptide_id)
                
        logger.info(f"Already completed: {len(completed)}")
        logger.info(f"Pending: {len(pending)}")
        
        if not pending:
            logger.info("All CG structures already generated!")
            return True
            
        # Run batch Martinize2
        cmd = [
            "python", "scripts/batch_martinize2.py",
            "--input-dir", self.config['pdb_dir'],
            "--output-dir", self.config['cg_pdb_dir']
        ]
        
        if self.config['test_mode']:
            cmd.extend(["--test-run"])
            
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=self.project_root)
        
        return result.returncode == 0
        
    def setup_simulations(self, peptide_ids):
        """Module 3: Setup MD simulations."""
        logger.info("="*60)
        logger.info("MODULE 3: MD Simulation Setup")
        logger.info("="*60)
        
        logger.info("For full RBC membrane simulations:")
        logger.info("1. Use CHARMM-GUI Martini Maker")
        logger.info("2. Upload peptide structures from data/processed/cg_pdb/")
        logger.info("3. Select Martini 3 force field")
        logger.info("4. Build asymmetric RBC membrane")
        logger.info("")
        
        # For demo, create simple test systems
        if self.config['test_mode']:
            logger.info("Creating demo simulation for SOLVIA_1...")
            
            cmd = [
                "python", "scripts/quick_martini3_demo.py",
                "--peptide-id", "SOLVIA_1",
                "--output-dir", "simulations/demo_test"
            ]
            
            result = subprocess.run(cmd, cwd=self.project_root)
            
            if result.returncode == 0:
                logger.info("Demo simulation setup complete!")
                logger.info("To run: cd simulations/demo_test && ./run_demo.sh")
            
        return True
        
    def extract_features(self, peptide_ids):
        """Module 4: Extract features from MD simulations."""
        logger.info("="*60)
        logger.info("MODULE 4: Feature Extraction")
        logger.info("="*60)
        
        # For demo, extract features even without full simulations
        cmd = [
            "python", "scripts/extract_md_features.py",
            "--batch",
            "--simulation-base-dir", self.config['simulation_dir'],
            "--output-dir", self.config['feature_dir']
        ]
        
        # Create peptide list file
        peptide_list_file = Path(self.config['feature_dir']) / "peptide_list.txt"
        with open(peptide_list_file, 'w') as f:
            for pid in peptide_ids:
                f.write(f"{pid}\n")
                
        cmd.extend(["--peptide-list", str(peptide_list_file)])
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=self.project_root)
        
        # Check if features were extracted
        feature_file = Path(self.config['feature_dir']) / "all_features.csv"
        if feature_file.exists():
            df = pd.read_csv(feature_file)
            logger.info(f"Extracted features for {len(df)} peptides")
            logger.info(f"Feature dimensions: {df.shape}")
        
        return result.returncode == 0
        
    def train_model(self):
        """Module 5: Train ML model."""
        logger.info("="*60)
        logger.info("MODULE 5: ML Model Training")
        logger.info("="*60)
        
        feature_file = Path(self.config['feature_dir']) / "all_features.csv"
        
        if not feature_file.exists():
            logger.error("No features found! Run feature extraction first.")
            return False
            
        # Train model
        cmd = [
            "python", "scripts/train_hemolysis_model.py",
            "--features", str(feature_file),
            "--output-dir", self.config['results_dir'],
            "--model-dir", self.config['model_dir'],
            "--create-labels"  # For demo, create synthetic labels
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=self.project_root)
        
        # Check results
        if result.returncode == 0:
            model_info = Path(self.config['model_dir']) / "model_info.json"
            if model_info.exists():
                with open(model_info, 'r') as f:
                    info = json.load(f)
                    
                logger.info("\nModel Performance:")
                for metric, value in info['performance'].items():
                    logger.info(f"  {metric}: {value:.3f}")
                    
        return result.returncode == 0
        
    def generate_report(self):
        """Generate pipeline execution report."""
        logger.info("="*60)
        logger.info("Generating Pipeline Report")
        logger.info("="*60)
        
        report = {
            'timestamp': self.timestamp,
            'config': self.config,
            'modules': {}
        }
        
        # Check each module's output
        # Module 1: Structure prediction
        pdb_files = list(Path(self.config['pdb_dir']).glob("*.pdb"))
        report['modules']['structure_prediction'] = {
            'completed': len(pdb_files),
            'files': [f.name for f in pdb_files[:10]]  # First 10
        }
        
        # Module 2: Coarse-graining
        cg_files = list(Path(self.config['cg_pdb_dir']).glob("*_cg.pdb"))
        report['modules']['coarse_graining'] = {
            'completed': len(cg_files),
            'files': [f.name for f in cg_files[:10]]
        }
        
        # Module 4: Features
        feature_file = Path(self.config['feature_dir']) / "all_features.csv"
        if feature_file.exists():
            df = pd.read_csv(feature_file)
            report['modules']['feature_extraction'] = {
                'n_samples': len(df),
                'n_features': len(df.columns) - 1
            }
            
        # Module 5: Model
        model_info_file = Path(self.config['model_dir']) / "model_info.json"
        if model_info_file.exists():
            with open(model_info_file, 'r') as f:
                model_info = json.load(f)
                report['modules']['ml_model'] = model_info
                
        # Save report
        report_file = Path(self.config['results_dir']) / f"pipeline_report_{self.timestamp}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
            
        logger.info(f"Report saved to {report_file}")
        
        # Print summary
        logger.info("\nPipeline Summary:")
        logger.info(f"  Structures predicted: {report['modules']['structure_prediction']['completed']}")
        logger.info(f"  CG structures: {report['modules']['coarse_graining']['completed']}")
        
        if 'feature_extraction' in report['modules']:
            logger.info(f"  Features extracted: {report['modules']['feature_extraction']['n_samples']} samples")
            
        if 'ml_model' in report['modules']:
            perf = report['modules']['ml_model'].get('performance', {})
            logger.info(f"  Model AUC: {perf.get('auc', 'N/A'):.3f}")
            
        return report
        
    def run(self):
        """Run the complete pipeline."""
        logger.info("="*60)
        logger.info("SOLVIA PIPELINE STARTED")
        logger.info("="*60)
        logger.info(f"Timestamp: {self.timestamp}")
        logger.info(f"Test mode: {self.config['test_mode']}")
        
        # Load peptide list
        if Path(self.config['peptide_list']).exists():
            with open(self.config['peptide_list'], 'r') as f:
                peptide_ids = [line.strip() for line in f if line.strip()]
        else:
            # Create test list
            peptide_ids = [f"SOLVIA_{i}" for i in range(1, 11)]
            with open(self.config['peptide_list'], 'w') as f:
                for pid in peptide_ids:
                    f.write(f"{pid}\n")
                    
        logger.info(f"Processing {len(peptide_ids)} peptides")
        
        # Run modules
        success = True
        
        # Module 1: Structure prediction
        if success and not self.config.get('skip_structure_prediction', False):
            success = self.run_structure_prediction(peptide_ids)
            
        # Module 2: Coarse-graining
        if success and not self.config.get('skip_coarse_graining', False):
            success = self.run_coarse_graining(peptide_ids)
            
        # Module 3: Simulation setup
        if success and not self.config.get('skip_simulation_setup', False):
            success = self.setup_simulations(peptide_ids)
            
        # Module 4: Feature extraction
        if success and not self.config.get('skip_feature_extraction', False):
            success = self.extract_features(peptide_ids)
            
        # Module 5: Model training
        if success and not self.config.get('skip_model_training', False):
            success = self.train_model()
            
        # Generate report
        report = self.generate_report()
        
        logger.info("="*60)
        logger.info("PIPELINE COMPLETED!")
        logger.info("="*60)
        
        return success


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run SOLVIA pipeline")
    parser.add_argument('--config', help='Configuration file')
    parser.add_argument('--test', action='store_true', help='Run in test mode')
    parser.add_argument('--skip', nargs='+', 
                       choices=['structure_prediction', 'coarse_graining', 
                               'simulation_setup', 'feature_extraction', 'model_training'],
                       help='Skip specific modules')
    
    args = parser.parse_args()
    
    # Create pipeline
    pipeline = SolviaPipeline(args.config)
    
    # Apply command line options
    if args.test:
        pipeline.config['test_mode'] = True
        
    if args.skip:
        for module in args.skip:
            pipeline.config[f'skip_{module}'] = True
            
    # Run pipeline
    success = pipeline.run()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()