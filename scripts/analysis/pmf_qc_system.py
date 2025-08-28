#!/usr/bin/env python3
"""
SOLVIA PMF Quality Control System
Comprehensive QC checks for PMF calculations
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from datetime import datetime
import json

class PMFQualityControl:
    """Complete QC system for PMF calculations"""
    
    def __init__(self, pmf_dir, config_file=None):
        self.pmf_dir = Path(pmf_dir)
        self.load_config(config_file)
        self.qc_results = {
            "timestamp": datetime.now().isoformat(),
            "pmf_directory": str(pmf_dir),
            "checks": {},
            "passed": False,
            "summary": {}
        }
    
    def load_config(self, config_file=None):
        """Load QC configuration"""
        if config_file:
            config_path = Path(config_file)
        else:
            config_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
        
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                self.qc_config = config.get('pmf', {}).get('qc', {})
        else:
            # Default QC parameters
            self.qc_config = {
                'min_neighbor_overlap': 0.10,
                'target_overlap': 0.20,
                'min_ess_frames': 200,
                'half_time_tolerance': 2.0,
                'replicate_tolerance': 2.0,
                'max_windows': 22,
                'max_time_per_window': 30
            }
    
    def check_window_files(self):
        """Check that all necessary window files exist"""
        window_dirs = list(self.pmf_dir.glob("windows/z_*"))
        
        missing_files = []
        incomplete_windows = []
        
        for window_dir in window_dirs:
            required_files = ['pullf.xvg', 'pullx.xvg', 'umbrella.tpr', 'umbrella.log']
            
            for req_file in required_files:
                if not (window_dir / req_file).exists():
                    missing_files.append(str(window_dir / req_file))
                    if window_dir not in incomplete_windows:
                        incomplete_windows.append(str(window_dir))
        
        self.qc_results["checks"]["file_completeness"] = {
            "n_windows": len(window_dirs),
            "incomplete_windows": incomplete_windows,
            "missing_files": missing_files,
            "passed": len(missing_files) == 0
        }
        
        return len(missing_files) == 0
    
    def check_window_overlap(self):
        """Check histogram overlap between adjacent windows"""
        metadata_file = self.pmf_dir / "pmf_metadata.yaml"
        
        if not metadata_file.exists():
            self.qc_results["checks"]["window_overlap"] = {
                "error": "Metadata file not found",
                "passed": False
            }
            return False
        
        with open(metadata_file, 'r') as f:
            metadata = yaml.safe_load(f)
        
        windows = metadata.get('windows', [])
        if not windows:
            return False
        
        overlap_data = []
        low_overlaps = []
        
        # Sort windows by center
        windows = [w for w in windows if w is not None]
        windows.sort(key=lambda x: x['center'], reverse=True)
        
        for i in range(len(windows) - 1):
            w1, w2 = windows[i], windows[i+1]
            
            # Load data
            data1 = self.load_pullx_data(self.pmf_dir / w1['pullx'])
            data2 = self.load_pullx_data(self.pmf_dir / w2['pullx'])
            
            if len(data1) > 0 and len(data2) > 0:
                overlap = self.calculate_overlap(data1, data2)
                
                overlap_data.append({
                    'window1': w1['center'],
                    'window2': w2['center'],
                    'overlap': overlap,
                    'passed': overlap >= self.qc_config['min_neighbor_overlap']
                })
                
                if overlap < self.qc_config['min_neighbor_overlap']:
                    low_overlaps.append({
                        'windows': [w1['center'], w2['center']],
                        'overlap': overlap,
                        'suggested_center': (w1['center'] + w2['center']) / 2
                    })
        
        mean_overlap = np.mean([d['overlap'] for d in overlap_data]) if overlap_data else 0
        min_overlap = min([d['overlap'] for d in overlap_data]) if overlap_data else 0
        
        self.qc_results["checks"]["window_overlap"] = {
            "mean_overlap": float(mean_overlap),
            "min_overlap": float(min_overlap),
            "target_overlap": self.qc_config['target_overlap'],
            "min_required": self.qc_config['min_neighbor_overlap'],
            "low_overlaps": low_overlaps,
            "n_failed": len(low_overlaps),
            "passed": len(low_overlaps) == 0
        }
        
        return len(low_overlaps) == 0
    
    def load_pullx_data(self, pullx_file):
        """Load z-coordinates from pullx file"""
        data = []
        
        if not pullx_file.exists():
            pullx_file = self.pmf_dir / pullx_file
        
        if pullx_file.exists():
            with open(pullx_file, 'r') as f:
                for line in f:
                    if not line.startswith(('#', '@')):
                        parts = line.split()
                        if len(parts) >= 2:
                            data.append(float(parts[1]))
        
        return np.array(data)
    
    def calculate_overlap(self, data1, data2):
        """Calculate histogram overlap between two data sets"""
        if len(data1) == 0 or len(data2) == 0:
            return 0.0
        
        # Create histograms
        min_val = min(data1.min(), data2.min())
        max_val = max(data1.max(), data2.max())
        bins = np.linspace(min_val, max_val, 50)
        
        hist1, _ = np.histogram(data1, bins=bins, density=True)
        hist2, _ = np.histogram(data2, bins=bins, density=True)
        
        # Calculate overlap
        overlap = np.minimum(hist1, hist2).sum() / min(hist1.sum(), hist2.sum())
        
        return overlap
    
    def check_convergence(self):
        """Check convergence by comparing first and second half of data"""
        window_dirs = list(self.pmf_dir.glob("windows/z_*"))
        
        convergence_data = []
        unconverged_windows = []
        
        for window_dir in window_dirs:
            pullx_file = window_dir / "pullx.xvg"
            if not pullx_file.exists():
                continue
            
            data = self.load_pullx_data(pullx_file)
            if len(data) < 100:
                continue
            
            # Split data in half
            mid_point = len(data) // 2
            first_half = data[:mid_point]
            second_half = data[mid_point:]
            
            # Calculate means and difference
            mean1 = np.mean(first_half)
            mean2 = np.mean(second_half)
            diff = abs(mean1 - mean2)
            
            # Extract window center from directory name
            window_center = float(window_dir.name.split('_')[1])
            
            converged = diff <= self.qc_config['half_time_tolerance'] / 100  # Convert to nm
            
            convergence_data.append({
                'window': window_center,
                'mean_first_half': mean1,
                'mean_second_half': mean2,
                'difference': diff,
                'converged': converged
            })
            
            if not converged:
                unconverged_windows.append(window_center)
        
        self.qc_results["checks"]["convergence"] = {
            "n_windows_checked": len(convergence_data),
            "n_unconverged": len(unconverged_windows),
            "unconverged_windows": unconverged_windows,
            "max_difference": max([d['difference'] for d in convergence_data]) if convergence_data else 0,
            "tolerance": self.qc_config['half_time_tolerance'] / 100,
            "passed": len(unconverged_windows) == 0
        }
        
        return len(unconverged_windows) == 0
    
    def check_effective_sample_size(self):
        """Check effective sample size for each window"""
        window_dirs = list(self.pmf_dir.glob("windows/z_*"))
        
        ess_data = []
        low_ess_windows = []
        
        for window_dir in window_dirs:
            pullx_file = window_dir / "pullx.xvg"
            if not pullx_file.exists():
                continue
            
            data = self.load_pullx_data(pullx_file)
            if len(data) == 0:
                continue
            
            # Calculate autocorrelation
            ess = self.calculate_ess(data)
            window_center = float(window_dir.name.split('_')[1])
            
            ess_data.append({
                'window': window_center,
                'ess': ess,
                'n_frames': len(data),
                'passed': ess >= self.qc_config['min_ess_frames']
            })
            
            if ess < self.qc_config['min_ess_frames']:
                low_ess_windows.append({
                    'window': window_center,
                    'ess': ess,
                    'required': self.qc_config['min_ess_frames']
                })
        
        mean_ess = np.mean([d['ess'] for d in ess_data]) if ess_data else 0
        min_ess = min([d['ess'] for d in ess_data]) if ess_data else 0
        
        self.qc_results["checks"]["effective_sample_size"] = {
            "mean_ess": float(mean_ess),
            "min_ess": float(min_ess),
            "required_ess": self.qc_config['min_ess_frames'],
            "low_ess_windows": low_ess_windows,
            "n_failed": len(low_ess_windows),
            "passed": len(low_ess_windows) == 0
        }
        
        return len(low_ess_windows) == 0
    
    def calculate_ess(self, data):
        """Calculate effective sample size using autocorrelation"""
        n = len(data)
        if n < 10:
            return n
        
        # Calculate autocorrelation
        mean = np.mean(data)
        c0 = np.var(data)
        
        if c0 == 0:
            return n
        
        # Calculate autocorrelation for different lags
        max_lag = min(n // 4, 1000)
        sum_acf = 1.0
        
        for lag in range(1, max_lag):
            c_lag = np.mean((data[:-lag] - mean) * (data[lag:] - mean))
            acf = c_lag / c0
            
            if acf < 0.05:  # Cutoff for negligible correlation
                break
            
            sum_acf += 2 * acf
        
        # Effective sample size
        ess = n / sum_acf
        
        return int(ess)
    
    def check_replicate_consistency(self):
        """Check consistency between replicates if available"""
        parent_dir = self.pmf_dir.parent
        replicate_dirs = list(parent_dir.glob("replicate_*"))
        
        if len(replicate_dirs) < 2:
            self.qc_results["checks"]["replicate_consistency"] = {
                "n_replicates": len(replicate_dirs),
                "message": "Single replicate only",
                "passed": True  # Not a failure if only one replicate
            }
            return True
        
        # Load analysis results for each replicate
        features = []
        for rep_dir in replicate_dirs:
            results_file = rep_dir / "pmf_analysis_results.yaml"
            if results_file.exists():
                with open(results_file, 'r') as f:
                    results = yaml.safe_load(f)
                    features.append(results.get('features', {}))
        
        if len(features) < 2:
            return True
        
        # Compare key features
        comparisons = []
        for key in ['delta_g_ads', 'delta_g_insert']:
            if key in features[0] and key in features[1]:
                if features[0][key] is not None and features[1][key] is not None:
                    diff = abs(features[0][key] - features[1][key])
                    comparisons.append({
                        'feature': key,
                        'rep1': features[0][key],
                        'rep2': features[1][key],
                        'difference': diff,
                        'passed': diff <= self.qc_config['replicate_tolerance']
                    })
        
        all_passed = all(c['passed'] for c in comparisons) if comparisons else True
        
        self.qc_results["checks"]["replicate_consistency"] = {
            "n_replicates": len(replicate_dirs),
            "comparisons": comparisons,
            "tolerance": self.qc_config['replicate_tolerance'],
            "passed": all_passed
        }
        
        return all_passed
    
    def check_computational_budget(self):
        """Check if computational budget limits were respected"""
        window_dirs = list(self.pmf_dir.glob("windows/z_*"))
        
        # Count windows
        n_windows = len(window_dirs)
        budget_exceeded = n_windows > self.qc_config['max_windows']
        
        # Check simulation times
        long_runs = []
        for window_dir in window_dirs:
            log_file = window_dir / "umbrella.log"
            if log_file.exists():
                # Parse simulation time from log
                sim_time = self.parse_simulation_time(log_file)
                if sim_time > self.qc_config['max_time_per_window']:
                    long_runs.append({
                        'window': window_dir.name,
                        'time_ns': sim_time
                    })
        
        self.qc_results["checks"]["computational_budget"] = {
            "n_windows": n_windows,
            "max_windows": self.qc_config['max_windows'],
            "windows_exceeded": budget_exceeded,
            "long_runs": long_runs,
            "passed": not budget_exceeded and len(long_runs) == 0
        }
        
        return not budget_exceeded and len(long_runs) == 0
    
    def parse_simulation_time(self, log_file):
        """Parse simulation time from GROMACS log file"""
        sim_time = 0
        
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    if "Time:" in line and "Core t" in line:
                        # Extract simulation time
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == "Time:":
                                sim_time = float(parts[i+1])
                                break
        except:
            pass
        
        return sim_time / 1000  # Convert ps to ns
    
    def generate_qc_report(self):
        """Generate comprehensive QC report"""
        print("\n" + "="*60)
        print("PMF QUALITY CONTROL REPORT")
        print("="*60)
        
        # Run all checks
        checks = [
            ("File Completeness", self.check_window_files),
            ("Window Overlap", self.check_window_overlap),
            ("Convergence", self.check_convergence),
            ("Effective Sample Size", self.check_effective_sample_size),
            ("Replicate Consistency", self.check_replicate_consistency),
            ("Computational Budget", self.check_computational_budget)
        ]
        
        all_passed = True
        for check_name, check_func in checks:
            print(f"\n{check_name}:")
            passed = check_func()
            all_passed = all_passed and passed
            
            if passed:
                print(f"  ✓ PASSED")
            else:
                print(f"  ✗ FAILED")
                
                # Print details of failure
                if check_name in ["Window Overlap", "Convergence", "Effective Sample Size"]:
                    details = self.qc_results["checks"][check_name.lower().replace(" ", "_")]
                    if "low_overlaps" in details:
                        for item in details["low_overlaps"]:
                            print(f"    - Windows {item['windows']}: overlap = {item['overlap']:.3f}")
                    elif "unconverged_windows" in details:
                        for window in details["unconverged_windows"]:
                            print(f"    - Window z={window:.3f} not converged")
                    elif "low_ess_windows" in details:
                        for item in details["low_ess_windows"]:
                            print(f"    - Window z={item['window']:.3f}: ESS = {item['ess']}")
        
        # Overall summary
        self.qc_results["passed"] = all_passed
        self.qc_results["summary"] = {
            "total_checks": len(checks),
            "passed_checks": sum(1 for c in self.qc_results["checks"].values() if c.get("passed", False)),
            "recommendation": "Ready for analysis" if all_passed else "Requires attention"
        }
        
        print("\n" + "-"*60)
        if all_passed:
            print("✓ ALL QUALITY CHECKS PASSED")
            print("  System ready for PMF analysis")
        else:
            print("✗ QUALITY ISSUES DETECTED")
            print("  Review failed checks and take corrective action")
        
        # Save report
        report_file = self.pmf_dir / "qc_full_report.yaml"
        with open(report_file, 'w') as f:
            yaml.dump(self.qc_results, f, default_flow_style=False)
        
        print(f"\nDetailed report saved to: {report_file}")
        
        return self.qc_results
    
    def suggest_corrections(self):
        """Suggest corrections for failed QC checks"""
        suggestions = []
        
        for check_name, check_data in self.qc_results["checks"].items():
            if not check_data.get("passed", True):
                if check_name == "window_overlap":
                    for low_overlap in check_data.get("low_overlaps", []):
                        suggestions.append({
                            "issue": "Low window overlap",
                            "location": f"Between z={low_overlap['windows'][0]:.3f} and z={low_overlap['windows'][1]:.3f}",
                            "action": f"Add window at z={low_overlap['suggested_center']:.3f}",
                            "command": f"python3 scripts/universal/run_single_window.py --center {low_overlap['suggested_center']:.3f}"
                        })
                
                elif check_name == "convergence":
                    for window in check_data.get("unconverged_windows", []):
                        suggestions.append({
                            "issue": "Poor convergence",
                            "location": f"Window z={window:.3f}",
                            "action": "Extend simulation by 10 ns",
                            "command": f"python3 scripts/universal/extend_window.py --window z_{window:+.3f} --add-ns 10"
                        })
                
                elif check_name == "effective_sample_size":
                    for item in check_data.get("low_ess_windows", []):
                        suggestions.append({
                            "issue": "Low ESS",
                            "location": f"Window z={item['window']:.3f}",
                            "action": "Extend simulation or check for sampling issues",
                            "command": f"python3 scripts/analysis/check_sampling.py --window z_{item['window']:+.3f}"
                        })
        
        if suggestions:
            print("\n" + "="*60)
            print("SUGGESTED CORRECTIONS")
            print("="*60)
            for i, suggestion in enumerate(suggestions, 1):
                print(f"\n{i}. {suggestion['issue']} at {suggestion['location']}")
                print(f"   Action: {suggestion['action']}")
                print(f"   Command: {suggestion['command']}")
        
        return suggestions

def main():
    parser = argparse.ArgumentParser(
        description="PMF Quality Control System"
    )
    parser.add_argument("pmf_dir", help="PMF directory to check")
    parser.add_argument("--config", help="Custom QC configuration file")
    parser.add_argument("--suggest", action="store_true",
                       help="Suggest corrections for failed checks")
    
    args = parser.parse_args()
    
    # Run QC
    qc = PMFQualityControl(args.pmf_dir, args.config)
    report = qc.generate_qc_report()
    
    # Suggest corrections if requested
    if args.suggest and not report["passed"]:
        qc.suggest_corrections()
    
    # Exit with error if QC failed
    if not report["passed"]:
        sys.exit(1)

if __name__ == "__main__":
    main()
