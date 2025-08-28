#!/usr/bin/env python3
"""
SOLVIA PMF Analysis Pipeline with MBAR/WHAM
Implements robust PMF analysis with bootstrap confidence intervals
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy.interpolate import PchipInterpolator
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

# Try to import pymbar
try:
    import pymbar
    from pymbar import MBAR
    HAS_MBAR = True
except ImportError:
    HAS_MBAR = False
    print("Warning: pymbar not installed. Install with: pip install pymbar")

class PMFAnalyzer:
    """Complete PMF analysis with MBAR/WHAM and feature extraction"""
    
    def __init__(self, pmf_dir, config=None):
        self.pmf_dir = Path(pmf_dir)
        
        # Try different locations for metadata file
        possible_metadata_files = []
        
        # If pmf_dir already contains a tag (like pmf_local), check directly
        if self.pmf_dir.name.startswith("pmf_"):
            possible_metadata_files.append(self.pmf_dir / "pmf_metadata.yaml")
        else:
            # Check subdirectories if we're in the pmf parent directory
            for subdir in sorted(self.pmf_dir.glob("pmf_*/")):
                if subdir.is_dir():
                    candidate = subdir / "pmf_metadata.yaml"
                    if candidate.exists():
                        possible_metadata_files.append(candidate)
            # Also check parent directory
            possible_metadata_files.append(self.pmf_dir / "pmf_metadata.yaml")
        
        # Additional fallback paths
        possible_metadata_files.extend([
            self.pmf_dir / "pmf_local" / "pmf_metadata.yaml",
            self.pmf_dir / "pmf_test" / "pmf_metadata.yaml",
            self.pmf_dir / "pmf_final" / "pmf_metadata.yaml",
            self.pmf_dir / "replicate_1" / "pmf_metadata.yaml",
        ])
        
        self.metadata_file = None
        for candidate in possible_metadata_files:
            if candidate.exists():
                self.metadata_file = candidate
                print(f"Found metadata file: {self.metadata_file}")
                break
        
        if not self.metadata_file:
            # List what was searched for debugging
            print("Searched in the following locations:")
            for path in possible_metadata_files[:5]:  # Show first 5
                print(f"  - {path}")
            raise FileNotFoundError(f"PMF metadata not found in {self.pmf_dir}")
        
        # Load metadata - handle numpy objects in YAML
        with open(self.metadata_file, 'r') as f:
            content = f.read()
            
            # First, try safe_load
            try:
                self.metadata = yaml.safe_load(content)
            except Exception:
                # If fails due to numpy objects, parse manually
                self.metadata = {}
                
                # Parse key fields manually from the YAML
                import re
                
                # Extract force constant
                match = re.search(r'^force_constant:\s*(\d+)', content, re.MULTILINE)
                if match:
                    self.metadata['force_constant'] = int(match.group(1))
                
                # Extract n_windows
                match = re.search(r'^n_windows:\s*(\d+)', content, re.MULTILINE)
                if match:
                    self.metadata['n_windows'] = int(match.group(1))
                
                # Extract production_time
                match = re.search(r'^production_time:\s*([\d.]+)', content, re.MULTILINE)
                if match:
                    self.metadata['production_time'] = float(match.group(1))
                
                # Extract window information
                self.metadata['windows'] = []
                lines = content.split('\n')
                for i, line in enumerate(lines):
                    if line == 'windows:':
                        # Parse the windows list
                        j = i + 1
                        current_window = None
                        while j < len(lines):
                            if lines[j].startswith('- center:'):
                                if current_window:
                                    self.metadata['windows'].append(current_window)
                                current_window = {}
                                center_match = re.search(r'center:\s*([-.\deE]+)', lines[j])
                                if center_match:
                                    current_window['center'] = float(center_match.group(1))
                            elif lines[j].startswith('  pullf:'):
                                if current_window:
                                    current_window['pullf'] = lines[j].split('pullf:')[1].strip()
                            elif lines[j].startswith('  pullx:'):
                                if current_window:
                                    current_window['pullx'] = lines[j].split('pullx:')[1].strip()
                            elif lines[j].startswith('  tpr:'):
                                if current_window:
                                    current_window['tpr'] = lines[j].split('tpr:')[1].strip()
                            elif not lines[j].startswith(' ') and not lines[j].startswith('-'):
                                # End of windows section
                                if current_window:
                                    self.metadata['windows'].append(current_window)
                                break
                            j += 1
                        break
        
        # Update pmf_dir to the directory containing metadata
        self.pmf_dir = self.metadata_file.parent
        
        self.config = config or self.load_analysis_config()
        self.windows = self.metadata['windows']
        self.k = self.metadata.get('force_constant', 900)  # kJ/mol/nm^2
        self.temperature = 310  # K
        self.kb = 0.008314462618  # kJ/mol/K
        self.beta = 1.0 / (self.kb * self.temperature)
        
    def load_analysis_config(self):
        """Load PMF analysis configuration"""
        config_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                return config.get('pmf', {}).get('analysis', {})
        return {}
    
    def load_umbrella_data(self):
        """Load all umbrella sampling data"""
        data = []
        centers = []
        
        for window in self.windows:
            if window is None:
                continue
            
            pullx_file = Path(window['pullx'])
            if not pullx_file.exists():
                pullx_file = self.pmf_dir / pullx_file
            
            if not pullx_file.exists():
                print(f"Warning: Missing pullx file: {pullx_file}")
                continue
            
            # Load data (skip comments)
            window_data = []
            with open(pullx_file, 'r') as f:
                for line in f:
                    if not line.startswith(('#', '@')):
                        parts = line.split()
                        if len(parts) >= 2:
                            # Time and z-coordinate
                            window_data.append(float(parts[1]))
            
            if window_data:
                data.append(np.array(window_data))
                centers.append(window['center'])
        
        return data, np.array(centers)
    
    def calculate_mbar(self, data, centers, bootstrap=False, n_bootstrap=0):
        """Calculate PMF using MBAR"""
        if not HAS_MBAR:
            print("MBAR not available, falling back to WHAM")
            return self.calculate_wham(data, centers, bootstrap, n_bootstrap)
        
        # Prepare data for MBAR
        n_windows = len(data)
        N_k = np.array([len(d) for d in data])
        N_total = N_k.sum()
        
        # Concatenate all data
        all_data = np.concatenate(data)
        
        # Compute bias matrix u_kn
        u_kn = np.zeros((n_windows, N_total))
        
        idx = 0
        for j, d in enumerate(data):
            for z_val in d:
                for k in range(n_windows):
                    u_kn[k, idx] = 0.5 * self.k * (z_val - centers[k])**2 * self.beta
                idx += 1
        
        # Initialize MBAR with error handling
        try:
            mbar = MBAR(u_kn, N_k, verbose=False, initialize='BAR')
        except Exception as e:
            print(f"Warning: BAR initialization failed ({e}), trying zeros...")
            try:
                mbar = MBAR(u_kn, N_k, verbose=False, initialize='zeros')
            except Exception as e2:
                print(f"Warning: zeros initialization failed ({e2}), using default...")
                mbar = MBAR(u_kn, N_k, verbose=False)
        
        # Calculate PMF on a grid
        z_min = min(d.min() for d in data)
        z_max = max(d.max() for d in data)
        z_grid = np.linspace(z_min - 0.2, z_max + 0.2, 200)
        n_bins = len(z_grid)
        
        # Create bias matrix for each grid point
        u_kn_grid = np.zeros((n_windows, n_bins))
        for i, z in enumerate(z_grid):
            for k in range(n_windows):
                u_kn_grid[k, i] = 0.5 * self.k * (z - centers[k])**2 * self.beta
        
        # Compute PMF using MBAR
        pmf_results = mbar.compute_pmf(u_kn_grid)
        pmf = pmf_results[0]
        
        # Convert to kJ/mol and normalize
        pmf = pmf / self.beta
        pmf = pmf - pmf[-10:].mean()  # Set bulk to zero
        
        # Bootstrap for uncertainty
        if bootstrap and n_bootstrap > 0:
            pmf_bootstrap = []
            for _ in range(n_bootstrap):
                # Resample data
                boot_data = []
                for d in data:
                    idx = np.random.choice(len(d), len(d), replace=True)
                    boot_data.append(d[idx])
                
                # Calculate PMF for bootstrap sample
                try:
                    pmf_boot = self.calculate_mbar(boot_data, centers, bootstrap=False)[1]
                    pmf_bootstrap.append(pmf_boot)
                except:
                    continue
            
            if pmf_bootstrap:
                pmf_std = np.std(pmf_bootstrap, axis=0)
            else:
                pmf_std = np.zeros_like(pmf)
        else:
            pmf_std = np.zeros_like(pmf)
        
        return z_grid, pmf, pmf_std
    
    def calculate_wham(self, data, centers, bootstrap=False, n_bootstrap=0):
        """Calculate PMF using WHAM (simplified implementation)"""
        # Create histograms
        z_min = min(d.min() for d in data)
        z_max = max(d.max() for d in data)
        n_bins = 100
        z_edges = np.linspace(z_min - 0.2, z_max + 0.2, n_bins + 1)
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2
        
        # Count observations in each bin for each window
        N_ki = np.zeros((len(data), n_bins))
        for k, d in enumerate(data):
            N_ki[k], _ = np.histogram(d, bins=z_edges)
        
        # Bias matrix
        V_ki = np.zeros((len(data), n_bins))
        for k in range(len(data)):
            for i in range(n_bins):
                V_ki[k, i] = 0.5 * self.k * (z_centers[i] - centers[k])**2
        
        # WHAM iteration
        tolerance = 1e-6
        max_iter = 1000
        
        # Initialize weights
        f_i = np.ones(n_bins)
        f_k = np.ones(len(data))
        
        for iteration in range(max_iter):
            f_i_old = f_i.copy()
            
            # Update f_i
            denominator = np.zeros(n_bins)
            for k in range(len(data)):
                N_k = len(data[k])
                denominator += N_k * f_k[k] * np.exp(-self.beta * V_ki[k])
            
            numerator = N_ki.sum(axis=0)
            f_i = numerator / (denominator + 1e-20)
            
            # Update f_k
            for k in range(len(data)):
                f_k[k] = 1.0 / np.sum(f_i * np.exp(-self.beta * V_ki[k]))
            
            # Check convergence
            if np.abs(f_i - f_i_old).max() < tolerance:
                break
        
        # Convert to PMF
        pmf = -np.log(f_i + 1e-20) / self.beta
        pmf = pmf - pmf[-10:].mean()  # Normalize
        
        # Bootstrap for uncertainty
        if bootstrap and n_bootstrap > 0:
            pmf_bootstrap = []
            for _ in range(n_bootstrap):
                # Resample data
                boot_data = []
                for d in data:
                    idx = np.random.choice(len(d), len(d), replace=True)
                    boot_data.append(d[idx])
                
                # Calculate PMF for bootstrap sample
                try:
                    pmf_boot = self.calculate_wham(boot_data, centers, bootstrap=False)[1]
                    pmf_bootstrap.append(pmf_boot)
                except:
                    continue
            
            if pmf_bootstrap:
                pmf_std = np.std(pmf_bootstrap, axis=0)
            else:
                pmf_std = np.zeros_like(pmf)
        else:
            pmf_std = np.zeros_like(pmf)
        
        return z_centers, pmf, pmf_std
    
    def extract_features(self, z_grid, pmf, pmf_std=None):
        """Extract key features from PMF profile"""
        features = {}
        
        # Find bulk region (usually z > 2.0 nm)
        bulk_mask = z_grid > 2.0
        if bulk_mask.any():
            pmf_bulk = pmf[bulk_mask].mean()
        else:
            pmf_bulk = pmf[-10:].mean()
        
        # ΔG_ads: difference between bulk and adsorbed state
        # Adsorbed state around z = 0.5-1.5 nm (just outside membrane)
        ads_mask = (z_grid > 0.5) & (z_grid < 1.5)
        if ads_mask.any():
            pmf_ads = pmf[ads_mask].min()
            features['delta_g_ads'] = pmf_ads - pmf_bulk
            features['z_ads'] = float(z_grid[ads_mask][np.argmin(pmf[ads_mask])])
        
        # ΔG_insert: difference between bulk and inserted state
        # Inserted state around z = -1.0 to 0 nm (inside membrane)
        insert_mask = (z_grid > -1.0) & (z_grid < 0)
        if insert_mask.any():
            pmf_insert = pmf[insert_mask].min()
            features['delta_g_insert'] = pmf_insert - pmf_bulk
            features['z_insert'] = float(z_grid[insert_mask][np.argmin(pmf[insert_mask])])
        
        # ΔG‡: barrier height
        # Find maximum between adsorbed and inserted states
        if 'z_ads' in features and 'z_insert' in features:
            barrier_mask = (z_grid > features['z_insert']) & (z_grid < features['z_ads'])
            if barrier_mask.any():
                pmf_barrier = pmf[barrier_mask].max()
                features['delta_g_barrier'] = pmf_barrier - pmf_ads
                features['z_barrier'] = float(z_grid[barrier_mask][np.argmax(pmf[barrier_mask])])
        
        # Add uncertainties if available
        if pmf_std is not None and len(pmf_std) == len(pmf):
            if ads_mask.any():
                features['delta_g_ads_std'] = float(pmf_std[ads_mask].mean())
            if insert_mask.any():
                features['delta_g_insert_std'] = float(pmf_std[insert_mask].mean())
        
        return features
    
    def calculate_convergence(self, data, centers):
        """Calculate PMF convergence over time"""
        convergence = []
        n_checkpoints = 5
        
        for frac in np.linspace(0.2, 1.0, n_checkpoints):
            # Use only fraction of data
            data_subset = []
            for d in data:
                n_use = int(len(d) * frac)
                if n_use > 10:  # Minimum data points
                    data_subset.append(d[:n_use])
            
            if len(data_subset) < 3:  # Need minimum windows
                continue
            
            try:
                # Calculate PMF for subset
                z_grid, pmf, _ = self.calculate_mbar(
                    data_subset, centers[:len(data_subset)], 
                    bootstrap=False
                )
                
                convergence.append({
                    'fraction': frac,
                    'n_samples': sum(len(d) for d in data_subset),
                    'z': z_grid.tolist(),
                    'pmf': pmf.tolist()
                })
            except:
                continue
        
        return convergence
    
    def calculate_overlap_matrix(self, data, centers):
        """Calculate window overlap matrix"""
        n_windows = len(centers)
        overlap_matrix = np.zeros((n_windows, n_windows))
        
        for i in range(n_windows):
            for j in range(n_windows):
                if i == j:
                    overlap_matrix[i, j] = 1.0
                else:
                    # Calculate overlap between distributions
                    # Simple approach: check if ranges overlap
                    range_i = (data[i].min(), data[i].max())
                    range_j = (data[j].min(), data[j].max())
                    
                    overlap_start = max(range_i[0], range_j[0])
                    overlap_end = min(range_i[1], range_j[1])
                    
                    if overlap_end > overlap_start:
                        # Estimate overlap fraction
                        overlap_range = overlap_end - overlap_start
                        total_range = max(range_i[1], range_j[1]) - min(range_i[0], range_j[0])
                        overlap_matrix[i, j] = overlap_range / total_range
        
        return overlap_matrix
    
    def plot_pmf(self, z, pmf, pmf_std, features, output_file):
        """Plot PMF profile with features"""
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot PMF
        ax.plot(z, pmf, 'b-', linewidth=2, label='PMF')
        if pmf_std is not None and len(pmf_std) == len(pmf):
            ax.fill_between(z, pmf - pmf_std, pmf + pmf_std, alpha=0.3)
        
        # Mark features
        if 'z_ads' in features:
            ax.axvline(features['z_ads'], color='g', linestyle='--', alpha=0.5, 
                      label=f"Adsorbed: {features['delta_g_ads']:.1f} kJ/mol")
        if 'z_insert' in features:
            ax.axvline(features['z_insert'], color='r', linestyle='--', alpha=0.5,
                      label=f"Inserted: {features['delta_g_insert']:.1f} kJ/mol")
        if 'z_barrier' in features:
            ax.axvline(features['z_barrier'], color='orange', linestyle='--', alpha=0.5,
                      label=f"Barrier: {features['delta_g_barrier']:.1f} kJ/mol")
        
        # Add membrane region
        ax.axvspan(-1.5, 1.5, alpha=0.1, color='gray', label='Membrane')
        
        ax.set_xlabel('z-coordinate (nm)', fontsize=12)
        ax.set_ylabel('PMF (kJ/mol)', fontsize=12)
        ax.set_title('Potential of Mean Force Profile', fontsize=14)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_overlap_heatmap(self, overlap_matrix, centers, output_file):
        """Plot window overlap heatmap"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        im = ax.imshow(overlap_matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Overlap', fontsize=12)
        
        # Set ticks
        n_windows = len(centers)
        tick_positions = np.arange(0, n_windows, max(1, n_windows // 10))
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels([f"{centers[i]:.2f}" for i in tick_positions], rotation=45)
        ax.set_yticklabels([f"{centers[i]:.2f}" for i in tick_positions])
        
        ax.set_xlabel('Window center (nm)', fontsize=12)
        ax.set_ylabel('Window center (nm)', fontsize=12)
        ax.set_title('Window Overlap Matrix', fontsize=14)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    def generate_analysis_report(self):
        """Generate complete PMF analysis report"""
        print("\n" + "="*60)
        print("PMF ANALYSIS - MBAR/WHAM")
        print("="*60)
        
        # Load data
        print("Loading umbrella sampling data...")
        data, centers = self.load_umbrella_data()
        print(f"Loaded {len(data)} windows")
        
        if len(data) < 3:
            print("Error: Not enough windows for analysis (minimum 3)")
            return None
        
        # Get method and bootstrap settings
        method = self.config.get('method', 'mbar')
        bootstrap_config = self.config.get('bootstrap', {})
        do_bootstrap = bootstrap_config.get('enabled', False)
        n_bootstrap = bootstrap_config.get('n_bootstrap', 0) if do_bootstrap else 0
        
        # Calculate PMF
        print(f"\nCalculating PMF using {method.upper()}...")
        if do_bootstrap and n_bootstrap > 0:
            print(f"Running {n_bootstrap} bootstrap samples for uncertainty estimation...")
        
        if method == 'mbar':
            z_grid, pmf, pmf_std = self.calculate_mbar(data, centers, do_bootstrap, n_bootstrap)
        else:
            z_grid, pmf, pmf_std = self.calculate_wham(data, centers, do_bootstrap, n_bootstrap)
        
        # Extract features
        print("\nExtracting PMF features...")
        features = self.extract_features(z_grid, pmf, pmf_std)
        
        # Calculate convergence
        print("Checking convergence...")
        convergence = self.calculate_convergence(data, centers)
        
        # Calculate overlap
        print("Calculating window overlap...")
        overlap_matrix = self.calculate_overlap_matrix(data, centers)
        
        # Generate plots
        print("\nGenerating plots...")
        plots_dir = self.pmf_dir / "analysis_plots"
        plots_dir.mkdir(exist_ok=True)
        
        # PMF plot
        self.plot_pmf(z_grid, pmf, pmf_std, features, 
                     plots_dir / "pmf_profile.png")
        
        # Overlap heatmap
        self.plot_overlap_heatmap(overlap_matrix, centers, 
                                 plots_dir / "overlap_matrix.png")
        
        # Save results
        results = {
            'method': method,
            'n_windows': len(data),
            'window_centers': centers.tolist(),
            'features': features,
            'pmf_data': {
                'z': z_grid.tolist(),
                'pmf': pmf.tolist(),
                'pmf_std': pmf_std.tolist()
            },
            'quality_metrics': {
                'mean_overlap': float(overlap_matrix[overlap_matrix < 1].mean()),
                'min_overlap': float(overlap_matrix[overlap_matrix > 0].min()),
                'convergence': convergence[-1]['pmf'][-1] if convergence else None
            },
            'bootstrap': {
                'enabled': do_bootstrap,
                'n_samples': n_bootstrap
            }
        }
        
        # Save to file
        results_file = self.pmf_dir / "pmf_analysis_results.yaml"
        with open(results_file, 'w') as f:
            yaml.dump(results, f, default_flow_style=False)
        
        # Print summary
        print("\n" + "="*60)
        print("PMF ANALYSIS COMPLETE")
        print("="*60)
        print(f"Method: {method.upper()}")
        print(f"Windows analyzed: {len(data)}")
        if do_bootstrap:
            print(f"Bootstrap samples: {n_bootstrap}")
        print(f"\nKey Features:")
        print(f"  ΔG_ads: {features.get('delta_g_ads', 'N/A'):.2f} kJ/mol")
        if features.get('delta_g_insert') is not None:
            print(f"  ΔG_insert: {features['delta_g_insert']:.2f} kJ/mol")
        if features.get('delta_g_barrier') is not None:
            print(f"  ΔG‡: {features['delta_g_barrier']:.2f} kJ/mol")
        print(f"\nQuality Metrics:")
        print(f"  Mean overlap: {overlap_matrix[overlap_matrix < 1].mean():.3f}")
        print(f"  Min overlap: {overlap_matrix[overlap_matrix > 0].min():.3f}")
        print(f"\nOutput files:")
        print(f"  Results: {results_file}")
        print(f"  Plots: {plots_dir}")
        
        return results

def main():
    parser = argparse.ArgumentParser(
        description="Analyze PMF from umbrella sampling with MBAR/WHAM"
    )
    parser.add_argument("pmf_dir", help="PMF directory with window data")
    parser.add_argument("--method", default="mbar", choices=["mbar", "wham"],
                       help="Analysis method")
    parser.add_argument("--bootstrap", type=int, default=0,
                       help="Number of bootstrap samples (0 = no bootstrap)")
    parser.add_argument("--no-bootstrap", action="store_true",
                       help="Disable bootstrap uncertainty estimation")
    parser.add_argument("--no-plots", action="store_true",
                       help="Skip plot generation")
    
    args = parser.parse_args()
    
    # Create analyzer
    try:
        analyzer = PMFAnalyzer(args.pmf_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Override config if needed
    if args.method:
        analyzer.config['method'] = args.method
    
    # Handle bootstrap configuration
    if args.no_bootstrap:
        analyzer.config['bootstrap'] = {'enabled': False, 'n_bootstrap': 0}
    elif args.bootstrap > 0:
        analyzer.config['bootstrap'] = {'enabled': True, 'n_bootstrap': args.bootstrap}
    else:
        analyzer.config['bootstrap'] = {'enabled': False, 'n_bootstrap': 0}
    
    # Run analysis
    results = analyzer.generate_analysis_report()
    
    if results:
        print("\n✓ Analysis completed successfully")
    else:
        print("\n✗ Analysis failed")
        sys.exit(1)

if __name__ == "__main__":
    main()