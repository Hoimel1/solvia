#!/usr/bin/env python3
from __future__ import annotations
"""
SOLVIA PMF Analysis Pipeline with MBAR/WHAM
Implements robust PMF analysis with bootstrap confidence intervals
"""

import sys
import yaml
import numpy as np
from typing import Optional
import matplotlib
# Force non-interactive backend for headless/CI execution
try:
    matplotlib.use('Agg', force=True)
except Exception:
    # Safe fallback; plotting will still attempt default backend
    pass
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy.interpolate import PchipInterpolator
from scipy import integrate
 
from scipy.special import logsumexp

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
        
        # Stash provided config (final config is set after metadata load)
        self.config = config if config is not None else {}
        # Try different locations for metadata file
        possible_metadata_files = []
        
        # 1. Check current directory first
        possible_metadata_files.append(self.pmf_dir / "pmf_metadata.yaml")
        
        # 2. If current dir starts with pmf_, check directly 
        if self.pmf_dir.name.startswith("pmf_"):
            possible_metadata_files.append(self.pmf_dir / "pmf_metadata.yaml")
            
        # 3. Check subdirectories of current directory
        for subdir in sorted(self.pmf_dir.glob("pmf_*/")):
            if subdir.is_dir():
                possible_metadata_files.append(subdir / "pmf_metadata.yaml")
                
        # 4. Check parent directory and its pmf_* subdirectories
        parent_dir = self.pmf_dir.parent
        if parent_dir.exists():
            # Check parent directory itself
            possible_metadata_files.append(parent_dir / "pmf_metadata.yaml")
            # Check sibling pmf_* directories 
            for sibling in sorted(parent_dir.glob("pmf_*/")):
                if sibling.is_dir() and sibling != self.pmf_dir:
                    possible_metadata_files.append(sibling / "pmf_metadata.yaml")
        
        # 5. Additional common fallback patterns
        possible_metadata_files.extend([
            self.pmf_dir / "pmf_local" / "pmf_metadata.yaml",
            self.pmf_dir / "pmf_test" / "pmf_metadata.yaml", 
            self.pmf_dir / "pmf_final" / "pmf_metadata.yaml",
            self.pmf_dir / "replicate_1" / "pmf_metadata.yaml",
            self.pmf_dir / "replicate_2" / "pmf_metadata.yaml",
        ])
        
        # 6. Go up one more level and check pmf directories there
        grandparent = self.pmf_dir.parent.parent
        if grandparent.exists() and grandparent.name != self.pmf_dir.parent.name:
            for pmf_dir in sorted(grandparent.glob("*/pmf*/pmf_metadata.yaml")):
                possible_metadata_files.append(pmf_dir)
        
        self.metadata_file = None
        for candidate in possible_metadata_files:
            if candidate.exists():
                self.metadata_file = candidate
                print(f"Found metadata file: {self.metadata_file}")
                break
        
        if not self.metadata_file:
            # List what was searched for debugging
            print("Searched in the following locations:")
            # Remove duplicates while preserving order
            seen = set()
            unique_paths = []
            for path in possible_metadata_files:
                if path not in seen:
                    seen.add(path)
                    unique_paths.append(path)
            
            for path in unique_paths[:10]:  # Show first 10 unique paths
                print(f"  - {path}")
            if len(unique_paths) > 10:
                print(f"  ... and {len(unique_paths) - 10} more locations")
            raise FileNotFoundError(f"PMF metadata not found. Please ensure pmf_metadata.yaml exists in one of the searched locations.")
        
        # Load metadata strictly via YAML
        with open(self.metadata_file, 'r') as f:
            try:
                self.metadata = yaml.safe_load(f)
            except Exception as e:
                raise RuntimeError(f"Failed to parse YAML metadata: {self.metadata_file}: {e}")
        
        # Update pmf_dir to the directory containing metadata
        self.pmf_dir = self.metadata_file.parent
        
        # Final config then derive parameters
        self.config = config or self.load_analysis_config()
        self.windows = self.metadata.get('windows', [])
        self.k = float(self.metadata.get('force_constant', 900))  # kJ/mol/nm^2
        temp_cfg = ((self.config.get('simulation') or {}).get('temperature')) if isinstance(self.config.get('simulation', {}), dict) else None
        self.temperature = float(self.metadata.get('temperature', temp_cfg if temp_cfg is not None else 310.0))
        self.R = 0.008314462618  # kJ/mol/K
        self.beta = 1.0 / (self.R * self.temperature)
        # Time filtering options
        tf = (self.config or {}).get('time_filter', {})
        self.begin_ps = tf.get('begin_ps', None)
        self.end_ps = tf.get('end_ps', None)
        self.dt_ps = tf.get('dt_ps', None)
        # Fractional discard as fallback (e.g., from QC defaults)
        self.discard_fraction = (self.config or {}).get('discard_fraction', None)
        # Bulk detection configuration
        bulk_cfg = (self.config or {}).get('bulk_detection', {})
        self.auto_bulk = bool(bulk_cfg.get('enabled', True))
        self.bulk_grad_thresh = float(bulk_cfg.get('grad_kj_per_nm', 2.0))
        self.bulk_min_bins = int(bulk_cfg.get('min_bins', 10))
        self.bulk_min_count = int(((self.config.get('plot') or {}).get('min_bin_count_for_ci') or 25))
        # WHAM weighting / grid refinement
        wham_cfg = (self.config or {}).get('wham', {})
        self.wham_use_ess_weighting = bool(wham_cfg.get('use_ess_weighting', True))
        adapt_cfg = (wham_cfg.get('adaptive_grid') or {})
        self.wham_adaptive_grid = bool(adapt_cfg.get('enabled', False))
        self.wham_adaptive_slope = float(adapt_cfg.get('slope_kj_per_nm', 50.0))
        self.wham_adaptive_passes = int(adapt_cfg.get('max_refine_passes', 1))
        # Symmetrization (optional)
        self.symmetrize = bool((self.config or {}).get('symmetrize', False))
        # Ensure grid spacing and bulk come from final config
        self.grid_spacing = (self.config.get('smoothing') or {}).get('grid_spacing', 0.01)
        self.bulk_min = self.config.get('z_bulk_min_nm', 2.4)
        # Optional bootstrap seed for reproducibility (use local RNG)
        bs_seed = ((self.config.get('bootstrap') or {}).get('seed')) if isinstance(self.config.get('bootstrap', {}), dict) else None
        seed_val = None
        if bs_seed is not None:
            try:
                seed_val = int(bs_seed)
            except Exception:
                seed_val = None
        self.rng = np.random.default_rng(seed_val)
        # QC/validation defaults
        qc_cfg = (self.config.get('qc') or {})
        self.ergodicity_block_r_max = float(qc_cfg.get('ergodicity_block_r_max', 0.2))
        
    def load_analysis_config(self):
        """Load PMF analysis configuration"""
        config_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                return config.get('pmf', {}).get('analysis', {})
        return {}
    
    def load_umbrella_data(self):
        """Load all umbrella sampling data with optional time filtering.

        Applies begin/end time cuts and optional thinning by dt_ps.
        If begin_ps is not provided but discard_fraction is set, drop the
        initial fraction per window based on its local time span.
        """
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

            # Load (time, z) pairs, skip comments
            t_list, z_list = [], []
            with open(pullx_file, 'r') as f:
                for line in f:
                    if line.startswith(('#', '@')):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            t = float(parts[0])
                            z = float(parts[1])
                        except ValueError:
                            continue
                        t_list.append(t)
                        z_list.append(z)

            if not t_list:
                continue

            t = np.asarray(t_list)
            z = np.asarray(z_list)

            # Determine local begin/end in ps
            b_ps = self.begin_ps
            e_ps = self.end_ps
            if b_ps is None and self.discard_fraction:
                frac = float(self.discard_fraction)
                if 0.0 < frac < 1.0:
                    t0, t1 = float(t[0]), float(t[-1])
                    b_ps = t0 + frac * (t1 - t0)

            # Apply begin/end filters
            mask = np.ones_like(t, dtype=bool)
            if b_ps is not None:
                mask &= (t >= b_ps)
            if e_ps is not None:
                mask &= (t <= e_ps)
            t = t[mask]
            z = z[mask]

            # Optional thinning by dt_ps (in ps), time-based to preserve independence
            if self.dt_ps and len(t) > 1:
                dt_target = float(self.dt_ps)
                kept_idx = []
                next_keep = t[0]
                for i in range(len(t)):
                    if t[i] + 1e-9 >= next_keep:
                        kept_idx.append(i)
                        next_keep = t[i] + dt_target
                t = t[kept_idx]
                z = z[kept_idx]

            if z.size:
                data.append(z)
                centers.append(window['center'])

        # Raise error if too few windows (<50% of expected)
        if len(data) < len(self.windows) / 2:
            raise RuntimeError("Too few pullx files found")
        return data, np.array(centers)

    def estimate_ess_per_window(self, data):
        """Estimate statistical inefficiency g and ESS per window using pymbar.timeseries.

        Returns two numpy arrays: g_k (>=1.0) and ess_k.
        """
        try:
            from pymbar import timeseries
        except Exception:
            timeseries = None

        g_list = []
        ess_list = []
        for d in data:
            try:
                if timeseries is not None:
                    g = float(timeseries.statistical_inefficiency(d))
                else:
                    # Simple fallback: integrated autocorrelation of mean-removed series
                    x = np.asarray(d, dtype=float)
                    x = x - np.mean(x)
                    n = len(x)
                    if n < 2:
                        g = 1.0
                    else:
                        var = np.var(x)
                        if var <= 0:
                            g = 1.0
                        else:
                            # compute until correlation drops below 0
                            max_lag = min(1000, n // 2)
                            acsum = 0.0
                            for lag in range(1, max_lag + 1):
                                c = float(np.dot(x[:-lag], x[lag:]) / (n - lag))
                                rho = c / var
                                if rho <= 0:
                                    break
                                acsum += rho
                            tau_int = 1.0 + 2.0 * acsum
                            g = max(1.0, float(tau_int))
                if not np.isfinite(g) or g < 1.0:
                    g = 1.0
            except Exception:
                g = 1.0
            g_list.append(g)
            ess_list.append(len(d) / g)
        return np.array(g_list, dtype=float), np.array(ess_list, dtype=float)

    def _detect_bulk_mask(self, z_centers: np.ndarray, pmf: np.ndarray, counts: Optional[np.ndarray] = None):
        """Detect a plateau region at large z to define bulk normalization.

        Criteria:
          - Use the highest-z contiguous segment where |dPMF/dz| < grad_thresh
          - Segment length >= min_bins
          - If counts provided: require counts >= bulk_min_count
        Falls back to z >= self.bulk_min.
        """
        z = np.asarray(z_centers)
        g = np.asarray(pmf)
        if z.size < 3:
            return z >= getattr(self, 'bulk_min', 2.4)
        # approximate gradient of PMF vs z
        dg = np.gradient(g, z)
        ok = np.abs(dg) <= self.bulk_grad_thresh
        if counts is not None:
            ok = np.logical_and(ok, np.asarray(counts) >= self.bulk_min_count)
        # Focus on right side (largest z). Find longest contiguous True block at right end
        idx_sorted = np.argsort(z)
        z_sorted = z[idx_sorted]
        ok_sorted = ok[idx_sorted]
        # rightmost contiguous block
        right_true = np.where(ok_sorted)[0]
        mask_out = np.zeros_like(ok_sorted, dtype=bool)
        if right_true.size > 0:
            # take the trailing block ending at the last index
            last_idx = len(ok_sorted) - 1
            j = last_idx
            run = 0
            while j >= 0 and ok_sorted[j]:
                mask_out[j] = True
                run += 1
                j -= 1
            if run < self.bulk_min_bins:
                mask_out[:] = False
        # Map back to original order
        mask_unsorted = np.zeros_like(mask_out)
        mask_unsorted[idx_sorted] = mask_out
        # Fallback: require both z >= bulk_min and sufficient counts if provided
        if not np.any(mask_unsorted):
            if counts is not None:
                return np.logical_and(z >= getattr(self, 'bulk_min', 2.4), np.asarray(counts) >= self.bulk_min_count)
            return z >= getattr(self, 'bulk_min', 2.4)
        return mask_unsorted

    def _bulk_mask(self, z_centers: np.ndarray, pmf: np.ndarray, counts: Optional[np.ndarray] = None):
        """Flexible bulk-region detection with config options.

        bulk_detection config options:
          - mode: 'auto' (default) -> gradient-based _detect_bulk_mask
          - mode: 'fixed_z' with 'z_min_nm' -> z >= z_min_nm
          - mode: 'percentile' with 'top_percent' -> top X% highest z as bulk
        """
        cfg = (self.config.get('bulk_detection') or {})
        mode = str(cfg.get('mode', 'auto')).lower()
        z = np.asarray(z_centers)
        if mode == 'fixed_z':
            zmin = float(cfg.get('z_min_nm', getattr(self, 'bulk_min', 2.4)))
            return z >= zmin
        if mode == 'percentile':
            perc = float(cfg.get('top_percent', 10.0))
            if not (0.0 < perc <= 100.0):
                perc = 10.0
            thresh = np.percentile(z, 100.0 - perc)
            return z >= thresh
        return self._detect_bulk_mask(z_centers, pmf, counts)

    def _symmetrize_profile(self, z: np.ndarray, y: np.ndarray):
        z = np.asarray(z)
        y = np.asarray(y)
        if z.min() >= 0:  # nothing to symmetrize
            return y.copy()
        y_sym = y.copy()
        dz = float(np.median(np.diff(z))) if len(z) > 1 else 1e-6
        tol = dz * 0.51
        for i, zi in enumerate(z):
            if zi < -tol:
                # find matching +z index
                j = np.argmin(np.abs(z - (-zi)))
                if abs(z[j] + zi) <= tol:
                    y_sym[i] = y_sym[j] = 0.5 * (y[i] + y[j])
        return y_sym
    
    def calculate_mbar(self, data, centers, bootstrap=False, n_bootstrap=0):
        """Calculate PMF using MBAR

        Notes on performance/stability:
        - MBAR scales with total sample count across all windows. To avoid
          excessive memory/CPU use (and apparent hangs), we adaptively thin
          the input data so that the total number of samples fed to MBAR does
          not exceed a safe threshold. This preserves the shape of the z
          distribution while making the solve tractable.
        - Bootstrap resamples may be additionally thinned via the
          'bootstrap_thin' config or CLI flag.
        """
        if not HAS_MBAR:
            print("MBAR not available, falling back to WHAM")
            return self.calculate_wham(data, centers, bootstrap, n_bootstrap)
        
        # Prepare data for MBAR
        n_windows = len(data)
        # Adaptive thinning: cap total samples for MBAR to avoid huge u_kn
        # matrices. Default cap ~200k samples across all windows.
        mbar_cap = int((self.config.get('mbar') or {}).get('max_samples', 200_000))
        lengths = np.array([len(d) for d in data], dtype=int)
        total = int(lengths.sum())
        thin_factor = 1
        if total > mbar_cap and total > 0:
            # Choose a global thin factor that reduces to <= cap
            thin_factor = int(np.ceil(total / float(mbar_cap)))
            print(f"[MBAR] Thinning data by factor {thin_factor} to cap {mbar_cap} samples (from {total}).")
            data = [d[::thin_factor] for d in data]
            lengths = np.array([len(d) for d in data], dtype=int)

        N_k = lengths.copy()
        N_total = int(N_k.sum())
        
        # Concatenate all data
        all_data = np.concatenate(data)
        print(f"[MBAR] Preparing energies (states={n_windows}, samples={N_total})...")
        
        # Compute bias matrix u_kn
        # Compute bias matrix u_kn (vectorized)
        # u_kn[k, n] = beta * 1/2 k_spring * (z_n - center_k)^2
        # Create a single z array for all samples and broadcast vs centers
        u_kn = 0.5 * self.k * self.beta * (centers[:, None] - all_data[None, :])**2
        
        # Initialize MBAR with error handling
        # MBAR solver configuration
        mbar_cfg = (self.config.get('mbar') or {})
        max_iter = int(mbar_cfg.get('max_iterations', 5000))
        rel_tol = float(mbar_cfg.get('relative_tolerance', 1e-7))
        # Use a conservative, broadly supported initialization
        try:
            mbar = MBAR(
                u_kn,
                N_k,
                maximum_iterations=max_iter,
                relative_tolerance=rel_tol,
                verbose=False,
                initialize='zeros',
            )
        except Exception as e:
            print(f"Warning: MBAR(init=zeros) failed ({e}), trying default...")
            try:
                mbar = MBAR(
                    u_kn,
                    N_k,
                    maximum_iterations=max_iter,
                    relative_tolerance=rel_tol,
                    verbose=False,
                )
            except Exception as e2:
                print(f"Warning: MBAR default init failed ({e2}), using minimal constructor...")
                mbar = MBAR(u_kn, N_k, verbose=False)
        print("[MBAR] Solver finished.")
        mbar_converged = getattr(mbar, 'converged', True)
        
        # Compute PMF using MBAR weights and a weighted histogram over z
        # 1) Reduced free energies of windows (dimensionless)
        f_k = mbar.f_k
        # 2) Per-sample normalization denominator using log-sum-exp for stability:
        #    denom_n = sum_k N_k * exp(f_k - u_kn)
        logs = (np.log(N_k + 1e-300)[:, None] + f_k[:, None] - u_kn)
        denom = np.exp(logsumexp(logs, axis=0))
        w_n = 1.0 / (denom + 1e-300)
        
        # 3) Probability density along z using MBAR expectations
        # Build fixed histogram grid and trim to contiguous support region
        z_min = min(d.min() for d in data)
        z_max = max(d.max() for d in data)
        n_bins = int(max(10, (z_max - z_min + 0.4) / self.grid_spacing))
        z_edges = np.linspace(z_min - 0.2, z_max + 0.2, n_bins + 1)
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2.0
        nb = len(z_edges) - 1

        # Assign samples to bins and get raw counts
        inds = np.digitize(all_data, z_edges) - 1
        inds = np.clip(inds, 0, nb - 1)
        counts_all = np.bincount(inds, minlength=nb)

        # Define a support mask based on minimum counts for stability; fallback to >0
        plot_cfg = (self.config.get('plot') or {})
        min_count_pmf = int(plot_cfg.get('pmf_min_bin_count', plot_cfg.get('min_bin_count_for_ci', 25)))
        mask_good = (counts_all >= min_count_pmf)
        if not mask_good.any():
            # Fallback: use non-zero support
            mask_good = (counts_all > 0)

        # Keep the longest contiguous True block to avoid sparse tails
        if mask_good.any():
            m = mask_good.astype(np.int8)
            diff = np.diff(np.concatenate(([0], m, [0])))
            starts = np.where(diff == 1)[0]
            ends = np.where(diff == -1)[0] - 1
            if len(starts):
                i_long = int(np.argmax(ends - starts + 1))
                keep = slice(starts[i_long], ends[i_long] + 1)
            else:
                keep = slice(0, nb)
        else:
            keep = slice(0, nb)

        # Trim grid and reindex
        z_edges = z_edges[keep.start: keep.stop + 1]
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2.0
        nb = len(z_edges) - 1
        inds = np.clip(inds - keep.start, 0, nb - 1)
        counts_all = counts_all[keep]

        # Use direct weighted histogram from MBAR weights (robust across PyMBAR versions)
        hist, _ = np.histogram(all_data, bins=z_edges, weights=w_n)
        p = hist / (hist.sum() + 1e-300)
        # Avoid log(0) by clipping
        p = np.clip(p, 1e-300, 1.0)
        # PMF (dimensionless) then convert to kJ/mol
        pmf = -np.log(p)
        pmf = pmf / self.beta
        # Normalize PMF by detected bulk plateau or fallback (use trimmed counts)
        # counts_all already computed on the trimmed grid above
        # Only normalize using bins with sufficient support
        support = counts_all >= self.bulk_min_count
        if getattr(self, 'auto_bulk', True):
            bulk_mask = self._bulk_mask(z_centers, pmf, counts_all)
        else:
            bulk_mask = z_centers >= getattr(self, 'bulk_min', 2.4)
        bulk_mask = np.logical_and(bulk_mask, support)
        if np.any(bulk_mask):
            pmf = pmf - pmf[bulk_mask].mean()
        else:
            # Fallback: use top 10% highest z among supported bins
            if np.any(support):
                z_sup = z_centers[support]
                thresh = np.percentile(z_sup, 90.0)
                # Map back to full mask
                mask_top = np.logical_and(support, z_centers >= thresh)
                pmf = pmf - pmf[mask_top].mean()
            else:
                # Last resort: subtract median to avoid huge offsets
                pmf = pmf - np.median(pmf)

        # Optional symmetrization for reporting
        if self.symmetrize:
            pmf = self._symmetrize_profile(z_centers, pmf)
        
        # Update outputs to reflect histogram grid
        z_grid = z_centers
        
        pmf_lower = None
        pmf_upper = None
        pmf_bootstrap_full = None
        # Bootstrap for uncertainty and confidence intervals
        if bootstrap and n_bootstrap > 0:
            thin = int(self.config.get('bootstrap_thin', 1) or 1)
            bs_cfg = (self.config.get('bootstrap') or {})
            bs_method = str(bs_cfg.get('method', 'block')).lower()
            pmf_bootstrap = []
            g_k, _ess = self.estimate_ess_per_window(data)
            for b in range(n_bootstrap):
                if (b + 1) % max(1, n_bootstrap // 10) == 0:
                    print(f"[MBAR] Bootstrap {b+1}/{n_bootstrap}...")
                # Resample data with optional blocking
                boot_data = []
                for k_idx, d in enumerate(data):
                    d_use = d[::thin] if thin > 1 else d
                    if bs_method in ('block', 'blocks', 'blocked'):
                        bl = int(max(2, round(2.0 * g_k[k_idx])))
                        if bl >= len(d_use):
                            idx = self.rng.choice(len(d_use), len(d_use), replace=True)
                            boot_data.append(d_use[idx])
                        else:
                            # Non-overlapping block bootstrap (NBB)
                            blocks = [d_use[i:i+bl] for i in range(0, len(d_use), bl)]
                            n_blocks = max(1, int(np.ceil(len(d_use) / bl)))
                            sel = self.rng.choice(len(blocks), n_blocks, replace=True)
                            series = np.concatenate([blocks[i] for i in sel])
                            series = series[:len(d_use)]
                            boot_data.append(series)
                    elif bs_method in ('stationary', 'stationary_bootstrap', 'sb'):
                        L_mean = max(2, int(round(g_k[k_idx])))
                        boot_series = self._stationary_bootstrap(d_use, L_mean)
                        boot_data.append(boot_series)
                    else:
                        idx = self.rng.choice(len(d_use), len(d_use), replace=True)
                        boot_data.append(d_use[idx])

                # Calculate PMF for bootstrap sample on fixed z grid
                try:
                    all_boot = np.concatenate(boot_data)
                    N_k_boot = np.array([len(d) for d in boot_data], dtype=int)
                    u_kn_boot = 0.5 * self.k * self.beta * (centers[:, None] - all_boot[None, :])**2
                    try:
                        mbar_boot = MBAR(
                            u_kn_boot,
                            N_k_boot,
                            maximum_iterations=int((self.config.get('mbar') or {}).get('max_iterations', 5000)),
                            relative_tolerance=float((self.config.get('mbar') or {}).get('relative_tolerance', 1e-7)),
                            verbose=False,
                            initialize='zeros'
                        )
                    except Exception:
                        mbar_boot = MBAR(u_kn_boot, N_k_boot, verbose=False)
                    logs_b = (np.log(N_k_boot + 1e-300)[:, None] + mbar_boot.f_k[:, None] - u_kn_boot)
                    denom_b = np.exp(logsumexp(logs_b, axis=0))
                    w_n_b = 1.0 / (denom_b + 1e-300)
                    # Weighted histogram directly via MBAR weights for bootstrap sample
                    hist_b, _ = np.histogram(all_boot, bins=z_edges, weights=w_n_b)
                    p_b = hist_b / (hist_b.sum() + 1e-300)
                    p_b = np.clip(p_b, 1e-300, 1.0)
                    pmf_b = -np.log(p_b) / self.beta
                    counts_b, _ = np.histogram(all_boot, bins=z_edges)
                    if getattr(self, 'auto_bulk', True):
                        bulk_mask_b = self._bulk_mask(z_centers, pmf_b, counts_b)
                    else:
                        bulk_mask_b = z_centers >= getattr(self, 'bulk_min', 2.4)
                    if np.any(bulk_mask_b):
                        pmf_b = pmf_b - pmf_b[bulk_mask_b].mean()
                    elif len(pmf_b) >= 10:
                        pmf_b = pmf_b - pmf_b[-10:].mean()
                    else:
                        pmf_b = pmf_b - pmf_b.mean()
                    pmf_bootstrap.append(pmf_b)
                except Exception:
                    continue
            
            if pmf_bootstrap:
                pmf_bootstrap_full = np.vstack(pmf_bootstrap)
                pmf_std = np.std(pmf_bootstrap_full, axis=0)
                pmf_lower = np.percentile(pmf_bootstrap_full, 2.5, axis=0)
                pmf_upper = np.percentile(pmf_bootstrap_full, 97.5, axis=0)
            else:
                pmf_std = np.zeros_like(pmf)
                pmf_lower = None
                pmf_upper = None
        else:
            pmf_std = np.zeros_like(pmf)
        # Fallback policy for MBAR results
        bad_numeric = not np.all(np.isfinite(pmf))
        do_fallback = False
        reason = None
        mbar_cfg = (self.config.get('mbar') or {})
        policy = str(mbar_cfg.get('fallback_policy', 'numeric_only')).lower()
        if policy == 'numeric_only':
            do_fallback = bad_numeric or (mbar_converged is False)
            reason = 'numeric' if bad_numeric else ('not_converged' if (mbar_converged is False) else None)
        elif policy == 'strict':
            max_range = float(mbar_cfg.get('pmf_range_kj', 200.0))
            min_allowed = float(mbar_cfg.get('pmf_min_kj', -100.0))
            bad_strict = ((pmf.max() - pmf.min()) > max_range) or (np.nanmin(pmf) < min_allowed)
            do_fallback = bad_numeric or bad_strict or (mbar_converged is False)
            reason = 'numeric' if bad_numeric else ('range' if bad_strict else ('not_converged' if (mbar_converged is False) else None))
        elif policy == 'never':
            do_fallback = False
        if do_fallback:
            print(f"[MBAR] Falling back to WHAM due to: {reason}")
            return self.calculate_wham(data, centers, bootstrap=False, n_bootstrap=0)
        return z_grid, pmf, pmf_std, pmf_lower, pmf_upper, pmf_bootstrap_full
    
    def calculate_wham(self, data, centers, bootstrap=False, n_bootstrap=0):
        """Calculate PMF using WHAM (simplified implementation)

        Uses a fixed histogram grid derived from the full dataset so that
        bootstrap replicates produce PMFs on the same z grid.
        """
        # Fixed histogram grid from full data
        z_min = min(d.min() for d in data)
        z_max = max(d.max() for d in data)
        n_bins = int(max(10, (z_max - z_min + 0.4) / self.grid_spacing))
        z_edges = np.linspace(z_min - 0.2, z_max + 0.2, n_bins + 1)
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2

        # Bias matrix will be built after potential trimming below

        last_iters = None
        def wham_solve(N_ki_local: np.ndarray, N_k_local: np.ndarray):
            nonlocal last_iters
            tolerance = 1e-6
            max_iter = 1000
            # Use probability per bin p_i and scale factors g_k for states
            p_i = np.ones(n_bins)
            g_k = np.ones(len(N_k_local))
            # Bin widths Δz for normalization integrals
            dz_vec = np.diff(z_edges)
            if dz_vec.size != n_bins:
                dz_vec = np.full(n_bins, float(np.median(np.diff(z_edges))))
            it = 0
            for it in range(max_iter):
                p_i_old = p_i.copy()
                denom = (N_k_local[:, None] * g_k[:, None] * exp_neg_beta_V).sum(axis=0) + 1e-20
                numerator = N_ki_local.sum(axis=0)
                p_i = numerator / denom
                # Include Δz in the state-normalization integral
                denom_k = (p_i[None, :] * exp_neg_beta_V * dz_vec[None, :]).sum(axis=1) + 1e-20
                g_k = 1.0 / denom_k
                if np.abs(p_i - p_i_old).max() < tolerance:
                    break
            last_iters = int(it + 1)
            pmf_local = -np.log(p_i + 1e-20) / self.beta
            # Normalize to bulk region
            bulk_mask = z_centers >= getattr(self, 'bulk_min', 2.4)
            if np.any(bulk_mask):
                pmf_local = pmf_local - pmf_local[bulk_mask].mean()
            elif len(pmf_local) >= 10:
                pmf_local = pmf_local - pmf_local[-10:].mean()
            else:
                pmf_local = pmf_local - pmf_local.mean()
            return pmf_local

        # Main solve on full data
        N_ki_main = np.zeros((len(data), n_bins), dtype=float)
        for k, d in enumerate(data):
            N_ki_main[k], _ = np.histogram(d, bins=z_edges)
        # Trim unsupported edges to avoid artificial vertical walls
        total_counts = N_ki_main.sum(axis=0)
        # define support as contiguous region covering all nonzero counts
        if (total_counts > 0).any():
            first = int(np.argmax(total_counts > 0))
            last = int(len(total_counts) - 1 - np.argmax((total_counts > 0)[::-1]))
        else:
            first, last = 0, len(total_counts) - 1
        # slice to support
        N_ki_main = N_ki_main[:, first:last+1]
        z_edges = z_edges[first:last+2]
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2
        n_bins = z_centers.size

        # Build bias matrix on trimmed grid
        V_ki_template = np.zeros((len(data), n_bins))
        for k in range(len(data)):
            V_ki_template[k, :] = 0.5 * self.k * (z_centers - centers[k])**2
        exp_neg_beta_V = np.exp(-self.beta * V_ki_template)

        N_k_main = np.array([len(d) for d in data], dtype=float)
        # Optional ESS weighting (Hub 2010): scale counts by 1/g_k
        if self.wham_use_ess_weighting:
            g_k, _ess = self.estimate_ess_per_window(data)
            scale = 1.0 / np.clip(g_k, 1.0, np.inf)
            N_ki_main = (N_ki_main.T * scale).T
            N_k_main = N_k_main * scale
        pmf = wham_solve(N_ki_main, N_k_main)

        # Optional adaptive grid refinement in steep regions
        if self.wham_adaptive_grid and self.wham_adaptive_passes > 0 and len(z_centers) > 5:
            eps = float(((self.config.get('wham') or {}).get('convergence_kj', 0.2)))
            min_counts = int(((self.config.get('wham') or {}).get('min_counts_refine', 20)))
            for _ in range(self.wham_adaptive_passes):
                slope = np.abs(np.gradient(pmf, z_centers))
                # counts support on current grid
                total_counts_cur = N_ki_main.sum(axis=0)
                cand = np.where((slope > self.wham_adaptive_slope) & (total_counts_cur >= min_counts))[0]
                if cand.size == 0:
                    break
                inserts = []
                for i in cand:
                    if i < len(z_centers) - 1:
                        inserts.append(0.5 * (z_centers[i] + z_centers[i+1]))
                if not inserts:
                    break
                pmf_old = pmf.copy()
                z_edges = np.sort(np.unique(np.concatenate([z_edges, inserts])))
                z_centers = (z_edges[:-1] + z_edges[1:]) / 2
                n_bins = len(z_centers)
                # Rebuild histograms and bias
                N_ki_main = np.zeros((len(data), n_bins), dtype=float)
                for k, d in enumerate(data):
                    N_ki_main[k], _ = np.histogram(d, bins=z_edges)
                if self.wham_use_ess_weighting:
                    g_k, _ess = self.estimate_ess_per_window(data)
                    scale = 1.0 / np.clip(g_k, 1.0, np.inf)
                    N_ki_main = (N_ki_main.T * scale).T
                    N_k_main = (np.array([len(d) for d in data], dtype=float) * scale)
                else:
                    N_k_main = np.array([len(d) for d in data], dtype=float)
                V_ki_template = np.zeros((len(data), n_bins))
                for k in range(len(data)):
                    V_ki_template[k, :] = 0.5 * self.k * (z_centers - centers[k])**2
                exp_neg_beta_V = np.exp(-self.beta * V_ki_template)
                pmf = wham_solve(N_ki_main, N_k_main)
                if np.max(np.abs(pmf - pmf_old)) < eps:
                    break

        # Normalize by detected bulk plateau if enabled (use only supported bins)
        if getattr(self, 'auto_bulk', True):
            counts_support = total_counts[first:last+1]
            support = counts_support >= self.bulk_min_count
            bulk_mask = self._bulk_mask(z_centers, pmf, counts_support)
            bulk_mask = np.logical_and(bulk_mask, support)
            if np.any(bulk_mask):
                pmf = pmf - pmf[bulk_mask].mean()
            else:
                if np.any(support):
                    z_sup = z_centers[support]
                    thresh = np.percentile(z_sup, 90.0)
                    mask_top = np.logical_and(support, z_centers >= thresh)
                    pmf = pmf - pmf[mask_top].mean()
                else:
                    pmf = pmf - np.median(pmf)
        # Optional symmetrization for reporting
        if self.symmetrize:
            pmf = self._symmetrize_profile(z_centers, pmf)

        # Expose convergence diagnostic
        try:
            self.wham_last_iters = int(last_iters) if last_iters is not None else None
        except Exception:
            self.wham_last_iters = None

        pmf_lower = None
        pmf_upper = None
        pmf_bootstrap_full = None
        if bootstrap and n_bootstrap > 0:
            thin = int(self.config.get('bootstrap_thin', 1) or 1)
            bs_cfg = (self.config.get('bootstrap') or {})
            bs_method = str(bs_cfg.get('method', 'standard')).lower()
            pmf_bootstrap = []
            # Precompute ESS-based block lengths if needed
            g_k, _ess = self.estimate_ess_per_window(data) if bs_method in ('block', 'blocks', 'blocked') else (None, None)
            for b in range(n_bootstrap):
                # Progress update every ~10%
                if (b + 1) % max(1, n_bootstrap // 10) == 0:
                    print(f"[WHAM] Bootstrap {b+1}/{n_bootstrap}...")
                boot_data = []
                for k, d in enumerate(data):
                    d_use = d[::thin] if thin > 1 else d
                    if bs_method in ('block', 'blocks', 'blocked') and g_k is not None:
                        # Non-overlapping block bootstrap with block length ~ 2*g_k frames
                        bl = int(max(2, round(2.0 * g_k[k])))
                        if bl >= len(d_use):
                            # fallback to simple resample
                            idx = self.rng.choice(len(d_use), len(d_use), replace=True)
                            boot_data.append(d_use[idx])
                        else:
                            # Build non-overlapping blocks
                            blocks = [d_use[i:i+bl] for i in range(0, len(d_use), bl)]
                            n_blocks = max(1, int(np.ceil(len(d_use) / bl)))
                            sel = self.rng.choice(len(blocks), n_blocks, replace=True)
                            boot_series = np.concatenate([blocks[i] for i in sel])
                            # Trim to original length
                            if boot_series.size > len(d_use):
                                boot_series = boot_series[:len(d_use)]
                            boot_data.append(boot_series)
                    elif bs_method in ('stationary', 'stationary_bootstrap', 'sb') and g_k is not None:
                        L_mean = max(2, int(round(g_k[k])))
                        boot_series = self._stationary_bootstrap(d_use, L_mean)
                        boot_data.append(boot_series)
                    else:
                        # IID bootstrap
                        idx = self.rng.choice(len(d_use), len(d_use), replace=True)
                        boot_data.append(d_use[idx])

                N_ki_boot = np.zeros((len(data), n_bins))
                for k, bdat in enumerate(boot_data):
                    N_ki_boot[k], _ = np.histogram(bdat, bins=z_edges)
                N_k_boot = np.array([len(bdat) for bdat in boot_data], dtype=float)
                try:
                    pmf_bootstrap.append(wham_solve(N_ki_boot, N_k_boot))
                except Exception:
                    continue
            if pmf_bootstrap:
                pmf_bootstrap_full = np.vstack(pmf_bootstrap)
                pmf_std = np.std(pmf_bootstrap_full, axis=0)
                pmf_lower = np.percentile(pmf_bootstrap_full, 2.5, axis=0)
                pmf_upper = np.percentile(pmf_bootstrap_full, 97.5, axis=0)
            else:
                pmf_std = np.zeros_like(pmf)
        else:
            pmf_std = np.zeros_like(pmf)
        return z_centers, pmf, pmf_std, pmf_lower, pmf_upper, pmf_bootstrap_full
    
    def extract_features(self, z_grid, pmf, pmf_std=None, counts=None):
        """Extract key features from PMF profile.

        If counts are provided, bulk determination is constrained to bins with
        sufficient support to avoid using empty, clipped tails.
        """
        features = {}
        
        # Find bulk region; constrain to supported bins if counts provided
        if getattr(self, 'auto_bulk', True):
            bulk_mask = self._bulk_mask(z_grid, pmf, counts)
        else:
            bulk_mask = (z_grid >= getattr(self, 'bulk_min', 2.4))
        if counts is not None:
            support = (np.asarray(counts) >= self.bulk_min_count)
            bulk_mask = np.logical_and(bulk_mask, support)
        if isinstance(bulk_mask, np.ndarray) and bulk_mask.any():
            pmf_bulk = float(np.mean(pmf[bulk_mask]))
        else:
            # PMF is already normalized to bulk earlier; default to 0.0
            pmf_bulk = 0.0
        
        # Adaptive feature ranges (use analysis.feature_params)
        feat_cfg = (self.config.get('feature_params') or {})
        # Adsorbed: search between z >= 0 and upper bound (default min(bulk_min, 2.0))
        ads_upper = float(feat_cfg.get('ads_max_z', min(getattr(self, 'bulk_min', 2.4), 2.0)))
        ads_region = (z_grid >= 0.0) & (z_grid <= ads_upper)
        if np.any(ads_region):
            idx = np.argmin(pmf[ads_region])
            z_ads = float(z_grid[ads_region][idx])
            pmf_ads = float(pmf[ads_region][idx])
            features['delta_g_ads'] = float(pmf_ads - pmf_bulk)
            features['z_ads'] = z_ads
        else:
            features['delta_g_ads'] = None
            features['z_ads'] = None
        
        # Inserted: search between lower bound (default -1.5 from config) and 0
        ins_lower = float(feat_cfg.get('insert_min_z', -1.5))
        insert_region = (z_grid >= ins_lower) & (z_grid <= 0.0)
        if np.any(insert_region):
            idx2 = np.argmin(pmf[insert_region])
            z_ins = float(z_grid[insert_region][idx2])
            pmf_ins = float(pmf[insert_region][idx2])
            features['delta_g_insert'] = float(pmf_ins - pmf_bulk)
            features['z_insert'] = z_ins
        else:
            features['delta_g_insert'] = None
            features['z_insert'] = None
        
        # ΔG‡: barrier height
        # Find maximum between adsorbed and inserted states
        if features.get('z_ads') is not None and features.get('z_insert') is not None:
            barrier_mask = (z_grid > features['z_insert']) & (z_grid < features['z_ads'])
            if barrier_mask.any():
                pmf_barrier = pmf[barrier_mask].max()
                features['delta_g_barrier'] = float(pmf_barrier - pmf_ads)
                features['z_barrier'] = float(z_grid[barrier_mask][np.argmax(pmf[barrier_mask])])
        
        # Add uncertainties if available
        if pmf_std is not None and len(pmf_std) == len(pmf):
            if np.any(ads_region):
                features['delta_g_ads_std'] = float(pmf_std[ads_region].mean())
            if np.any(insert_region):
                features['delta_g_insert_std'] = float(pmf_std[insert_region].mean())

        return features

    def compute_adsorption_metrics(self, z_grid, pmf, pmf_std=None, pmf_lower=None, pmf_upper=None, counts=None):
        """Compute adsorption partition integral Kp_ads and HC50 estimate.

        Uses config under analysis.adsorption with defaults:
          - window_mode: energy_band | pmf_negative | around_min | fixed
          - energy_band_kj: 3.0 (only for window_mode=energy_band)
          - z_low/z_high for fixed mode or fallback
          - area_per_lipid_nm2: 0.62
          - lp_star_range: [154, 515]
        """
        ads_cfg = (self.config or {}).get('adsorption', {})
        mode = str(ads_cfg.get('window_mode', 'fixed')).lower()
        z_low = float(ads_cfg.get('z_low', 0.5))
        z_high = float(ads_cfg.get('z_high', 1.6))
        area_per_lipid = float(ads_cfg.get('area_per_lipid_nm2', 0.62))
        lp_star = ads_cfg.get('lp_star_range', [154, 515])
        if isinstance(lp_star, (int, float)):
            lp_star = [float(lp_star), float(lp_star)]
        else:
            lp_star = [float(lp_star[0]), float(lp_star[1])] if lp_star else [154.0, 515.0]

        # Support mask (counts-based), membrane bounds and bulk side
        if counts is not None:
            support_min = int(ads_cfg.get('support_min_count', self.bulk_min_count))
            support = (np.asarray(counts) >= support_min)
        else:
            support = np.ones_like(z_grid, dtype=bool)
        mem = (ads_cfg.get('membrane') or {})
        z_lo_mem = float(mem.get('z_lo', -1.5))
        z_hi_mem = float(mem.get('z_hi', 1.5))
        margin = float(mem.get('margin_nm', 0.2))
        bulk_side_mask = (z_grid >= z_hi_mem + margin)

        # Determine window automatically if requested
        if mode == 'pmf_negative':
            raw = (pmf < 0.0) & bulk_side_mask & support
            if np.any(raw):
                m = raw.astype(np.int8)
                diff = np.diff(np.concatenate(([0], m, [0])))
                starts = np.where(diff == 1)[0]
                ends = np.where(diff == -1)[0] - 1
                if len(starts):
                    lengths = ends - starts + 1
                    i = int(np.argmax(lengths))
                    z_low = float(z_grid[starts[i]])
                    z_high = float(z_grid[ends[i]])
                else:
                    # no contiguous block
                    z_low = z_high = float(z_grid[0])
            else:
                z_low = z_high = float(z_grid[0])
        elif mode == 'around_min':
            # Choose adsorption minimum on the water side within feature range (0 .. ads_max_z)
            feat_cfg = (self.config.get('feature_params') or {})
            ads_upper = float(feat_cfg.get('ads_max_z', min(getattr(self, 'bulk_min', 2.4), 2.0)))
            cand = (z_grid >= 0.0) & (z_grid <= ads_upper)
            if counts is not None:
                cand = np.logical_and(cand, support)
            if np.any(cand):
                masked = np.where(cand, pmf, np.inf)
                i_min = int(np.argmin(masked))
            else:
                i_min = int(np.argmin(pmf))
            z0 = float(z_grid[i_min])
            width = float(ads_cfg.get('auto_width_nm', 0.6))
            zl = z0 - 0.5 * width
            zh = z0 + 0.5 * width
            # clamp to data extents
            z_low = float(max(z_grid.min(), zl))
            z_high = float(min(z_grid.max(), zh))
            # tighten to supported bins if available
            if counts is not None:
                mask_tmp = (z_grid >= z_low) & (z_grid <= z_high) & support
                idx = np.where(mask_tmp)[0]
                if idx.size:
                    z_low = float(z_grid[idx[0]])
                    z_high = float(z_grid[idx[-1]])
                else:
                    z_low = z_high = float(z_grid[0])
        elif mode in ('energy_band', 'energyband', 'band'):
            # Robust: define adsorption band by energy threshold around the minimum
            # 1) Find adsorption minimum within 0..ads_max_z on water side
            feat_cfg = (self.config.get('feature_params') or {})
            ads_upper = float(feat_cfg.get('ads_max_z', min(getattr(self, 'bulk_min', 2.4), 2.0)))
            region = (z_grid >= 0.0) & (z_grid <= ads_upper)
            if counts is not None:
                region = np.logical_and(region, support)
            if np.any(region):
                masked = np.where(region, pmf, np.inf)
                i_min = int(np.argmin(masked))
            else:
                i_min = int(np.argmin(pmf))
            w_min = float(pmf[i_min])
            band_kj = float(ads_cfg.get('energy_band_kj', 3.0))
            band_mask = (pmf <= (w_min + band_kj))
            # Restrict band to water side feature region
            band_mask = np.logical_and(band_mask, (z_grid >= 0.0) & (z_grid <= ads_upper))
            if counts is not None:
                band_mask = np.logical_and(band_mask, support)
            if np.any(band_mask):
                # Find contiguous segment around the minimum index
                m = band_mask.astype(np.int8)
                diff = np.diff(np.concatenate(([0], m, [0])))
                starts = np.where(diff == 1)[0]
                ends = np.where(diff == -1)[0] - 1
                # Select the block that contains i_min (or nearest)
                sel = None
                for s, e in zip(starts, ends):
                    if s <= i_min <= e:
                        sel = (s, e)
                        break
                if sel is None:
                    # Fallback: choose block whose center is closest to i_min
                    centers_idx = [(s + e) // 2 for s, e in zip(starts, ends)]
                    k = int(np.argmin([abs(c - i_min) for c in centers_idx])) if centers_idx else None
                    if k is not None:
                        sel = (int(starts[k]), int(ends[k]))
                if sel is not None:
                    z_low = float(z_grid[sel[0]])
                    z_high = float(z_grid[sel[1]])
                else:
                    # Last resort: small symmetric window around minimum
                    width = float(ads_cfg.get('auto_width_nm', 0.6))
                    z0 = float(z_grid[i_min])
                    z_low = float(max(z_grid.min(), z0 - 0.5 * width))
                    z_high = float(min(z_grid.max(), z0 + 0.5 * width))
            else:
                # No band after restrictions; fallback to narrow window around min
                width = float(ads_cfg.get('auto_width_nm', 0.6))
                z0 = float(z_grid[i_min])
                z_low = float(max(z_grid.min(), z0 - 0.5 * width))
                z_high = float(min(z_grid.max(), z0 + 0.5 * width))

        # Optional smoothing for integration stability
        smooth_cfg = (self.config or {}).get('smoothing', {})
        pmf_use = pmf
        if isinstance(smooth_cfg, dict) and smooth_cfg.get('method', '').lower() == 'pchip':
            try:
                interp = PchipInterpolator(z_grid, pmf, extrapolate=True)
                pmf_use = interp(z_grid)
            except Exception:
                pmf_use = pmf

        # Select adsorption region mask
        mask = (z_grid >= z_low) & (z_grid <= z_high)
        # Constrain to support if provided
        if counts is not None:
            mask = np.logical_and(mask, support)
        if not np.any(mask):
            # No supported region: Kp=0 by construction
            return {
                'z_low': z_low, 'z_high': z_high,
                'area_per_lipid_nm2': area_per_lipid,
                'lp_star_range': lp_star,
                'bilayer_factor': float(ads_cfg.get('bilayer_factor', 2.0)),
                'kp_ads_nm': 0.0,
                'kp_ads_nm_ci': [0.0, 0.0],
                'hc50_molar_range': None,
                'hc50_uM_range': None,
                'hc50_uM_ci': None,
                'dg_ads_eff': None
            }
        # Kp_ads integral (length in nm) with numerical stabilization
        # Kp = ∫(e^{x} − 1) dz, where x = −β W(z)
        x = -self.beta * pmf_use[mask]
        x0 = float(np.max(x))
        # ∫ e^{x} dz = e^{x0} ∫ e^{x−x0} dz; ∫ 1 dz computed on masked grid
        kp_scaled = float(integrate.trapezoid(np.exp(x - x0), z_grid[mask]))
        dz_len = float(integrate.trapezoid(np.ones_like(x), z_grid[mask]))
        kp_ads = float(kp_scaled * np.exp(x0) - dz_len)
        kp_ads = float(max(kp_ads, 0.0))
        kp_ci = None
        # Optional CI from pmf_lower/pmf_upper percentiles if provided
        if pmf_lower is not None and pmf_upper is not None:
            if len(pmf_lower) == len(pmf) == len(z_grid):
                # Upper PMF -> lower Kp; Lower PMF -> upper Kp (use same stabilization)
                x_low = -self.beta * pmf_upper[mask]
                x0_low = float(np.max(x_low))
                kp_scaled_low = float(integrate.trapezoid(np.exp(x_low - x0_low), z_grid[mask]))
                kp_low = float(kp_scaled_low * np.exp(x0_low) - dz_len)
                x_high = -self.beta * pmf_lower[mask]
                x0_high = float(np.max(x_high))
                kp_scaled_high = float(integrate.trapezoid(np.exp(x_high - x0_high), z_grid[mask]))
                kp_high = float(kp_scaled_high * np.exp(x0_high) - dz_len)
                kp_low = float(max(kp_low, 0.0)); kp_high = float(max(kp_high, 0.0))
                kp_ci = [kp_low, kp_high]

        # Critical surface coverage Gamma* from L/P* and area per lipid
        # HC50 = Gamma* / (0.602214 * Kp_ads)
        conv = 0.602214  # 1 M = 0.602214 nm^-3
        bilayer_factor = float(ads_cfg.get('bilayer_factor', 2.0))
        kp_eff = kp_ads * bilayer_factor if np.isfinite(kp_ads) else kp_ads
        gamma_min = 1.0 / (area_per_lipid * lp_star[1])
        gamma_max = 1.0 / (area_per_lipid * lp_star[0])
        # Guard against invalid Kp
        if (kp_eff is None) or (not np.isfinite(kp_eff)) or (kp_eff <= 0.0):
            return {
                'z_low': z_low, 'z_high': z_high,
                'area_per_lipid_nm2': area_per_lipid,
                'lp_star_range': [lp_star[0], lp_star[1]],
                'bilayer_factor': bilayer_factor,
                'kp_ads_nm': kp_ads,
                'kp_ads_nm_ci': kp_ci,
                'hc50_molar_range': None,
                'hc50_uM_range': None,
                'hc50_uM_ci': None,
                'dg_ads_eff': None
            }
        hc50_min = float(gamma_min / (conv * kp_eff))  # mol/L
        hc50_max = float(gamma_max / (conv * kp_eff))
        hc50_ci = None
        if kp_ci is not None:
            kp_low, kp_high = kp_ci
            hc_min_low = float(gamma_min / (conv * (kp_high * bilayer_factor)))
            hc_min_high = float(gamma_min / (conv * (kp_low * bilayer_factor)))
            hc_max_low = float(gamma_max / (conv * (kp_high * bilayer_factor)))
            hc_max_high = float(gamma_max / (conv * (kp_low * bilayer_factor)))
            hc50_ci = {
                'lp_star_min_uM': [hc_min_low * 1e6, hc_min_high * 1e6],
                'lp_star_max_uM': [hc_max_low * 1e6, hc_max_high * 1e6]
            }
        # Effective ΔG_ads from Kp with Δz_eff as region width
        dz_eff = max(1e-6, (z_high - z_low))
        dg_ads_eff = float(- (1.0 / self.beta) * np.log(max(1e-300, kp_eff / dz_eff)))

        return {
            'z_low': z_low, 'z_high': z_high,
            'area_per_lipid_nm2': area_per_lipid,
            'lp_star_range': [lp_star[0], lp_star[1]],
            'bilayer_factor': bilayer_factor,
            'kp_ads_nm': kp_eff,
            'kp_ads_nm_ci': kp_ci,
            'hc50_molar_range': [hc50_min, hc50_max],
            'hc50_uM_range': [hc50_min * 1e6, hc50_max * 1e6],
            'hc50_uM_ci': hc50_ci,
            'dg_ads_eff': dg_ads_eff
        }
    
    # Convergence check removed per user request: redundant and costly

    def estimate_block_stability(self, data, n_blocks: int = 5):
        """Per-window block stability metric.

        Split each window time series into n_blocks contiguous chunks and compute the
        variance of block means normalized by overall variance. Lower is better.
        Returns numpy array of R values (np.nan if insufficient data).
        """
        R = []
        for d in data:
            x = np.asarray(d, dtype=float)
            n = len(x)
            if n < n_blocks:
                R.append(np.nan)
                continue
            blocks = np.array_split(x, n_blocks)
            means = np.array([b.mean() for b in blocks])
            var_all = float(np.var(x))
            var_means = float(np.var(means))
            R.append(0.0 if var_all <= 0 else (var_means / var_all))
        return np.array(R, dtype=float)
    
    def calculate_overlap_matrix(self, data, centers):
        """Calculate window overlap matrix using histogram-based overlap (Runner style)"""
        n_windows = len(centers)
        overlap_matrix = np.zeros((n_windows, n_windows))
        for i in range(n_windows):
            for j in range(n_windows):
                if i == j:
                    overlap_matrix[i,j] = 1.0
                else:
                    zmin = min(data[i].min(), data[j].min())
                    zmax = max(data[i].max(), data[j].max())
                    # Guard against degenerate ranges and NaNs
                    if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax <= zmin:
                        overlap_matrix[i, j] = 0.0
                        continue
                    bins = np.linspace(zmin, zmax, 50)
                    h1,_ = np.histogram(data[i], bins=bins, density=True)
                    h2,_ = np.histogram(data[j], bins=bins, density=True)
                    bin_width = (zmax - zmin) / 50
                    overlap = np.minimum(h1,h2).sum() * bin_width
                    overlap_matrix[i,j] = max(0.0, min(1.0, overlap))
        return overlap_matrix

    def summarize_overlap(self, overlap_matrix, centers):
        """Summarize overlap in two ways: adjacent pairs and all pairs.

        - Adjacent: pairs of windows that are neighbors when centers are sorted.
        - All pairs: upper triangle (i<j) of the full overlap matrix.

        Returns a dict with mean/min for adjacent and all-pairs summaries.
        """
        # Indices for centers sorted by value
        order = np.argsort(centers)
        # Adjacent overlaps along the sorted order
        adjacent_values = []
        for a, b in zip(order[:-1], order[1:]):
            adjacent_values.append(float(overlap_matrix[a, b]))

        adjacent_values = np.array(adjacent_values) if adjacent_values else np.array([])

        # All pairs: use upper triangle without diagonal
        iu = np.triu_indices_from(overlap_matrix, k=1)
        all_pairs = overlap_matrix[iu]

        # Compute metrics safely
        result = {}
        if adjacent_values.size > 0:
            result['mean_overlap_adjacent'] = float(adjacent_values.mean())
            result['min_overlap_adjacent'] = float(adjacent_values.min())
        else:
            result['mean_overlap_adjacent'] = None
            result['min_overlap_adjacent'] = None

        if all_pairs.size > 0:
            # Mean over all off-diagonal pairs (can be quite small)
            result['mean_overlap_all_pairs'] = float(all_pairs.mean())
            # Min over positive off-diagonal entries (avoid exact zeros)
            positive = all_pairs[all_pairs > 0]
            result['min_overlap_all_pairs_nonzero'] = float(positive.min()) if positive.size else 0.0
        else:
            result['mean_overlap_all_pairs'] = None
            result['min_overlap_all_pairs_nonzero'] = None

        return result

    @staticmethod
    def _js_divergence(p: np.ndarray, q: np.ndarray) -> float:
        eps = 1e-12
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        p = p / (p.sum() + eps)
        q = q / (q.sum() + eps)
        m = 0.5 * (p + q)
        def _kl(a, b):
            a = np.clip(a, eps, 1.0)
            b = np.clip(b, eps, 1.0)
            return float(np.sum(a * (np.log(a) - np.log(b))))
        return 0.5 * (_kl(p, m) + _kl(q, m))

    def js_divergence_adjacent(self, data, centers, bins: int = 50):
        """Compute Jensen-Shannon divergence for adjacent window histograms.

        Returns list of {'z1','z2','js'} aligned with adjacent sorted centers and a summary dict.
        """
        centers = np.asarray(centers, dtype=float)
        order = np.argsort(centers)
        cs = centers[order]
        out = []
        for a, b in zip(order[:-1], order[1:]):
            zmin = min(data[a].min(), data[b].min())
            zmax = max(data[a].max(), data[b].max())
            if zmax <= zmin:
                continue
            edges = np.linspace(zmin, zmax, bins + 1)
            h1, _ = np.histogram(data[a], bins=edges, density=False)
            h2, _ = np.histogram(data[b], bins=edges, density=False)
            js = self._js_divergence(h1, h2)
            out.append({'z1': float(cs[np.where(order==a)[0][0]]),
                        'z2': float(cs[np.where(order==b)[0][0]]),
                        'js': float(js)})
        # summary
        js_vals = np.array([d['js'] for d in out], dtype=float) if out else np.array([], dtype=float)
        cfg = (self.config.get('qc') or {})
        js_max = float(cfg.get('js_divergence_max', 0.5))
        summary = {
            'js_max': js_max,
            'mean_js_adjacent': float(js_vals.mean()) if js_vals.size else None,
            'max_js_adjacent': float(js_vals.max()) if js_vals.size else None,
            'pass_fraction': float(np.mean(js_vals <= js_max)) if js_vals.size else None
        }
        return out, summary

    def _stationary_bootstrap(self, x: np.ndarray, L_mean: int) -> np.ndarray:
        """Stationary bootstrap (Politis-Romano) for correlated 1D series.

        - L_mean: expected block length (>=2). Uses geometric draws with p=1/L_mean.
        - Wraps indices circularly to preserve local correlation.
        """
        x = np.asarray(x, dtype=float)
        n = x.size
        if n == 0:
            return x.copy()
        L = max(2, int(L_mean))
        p = 1.0 / L
        out = np.empty(n, dtype=float)
        i = int(self.rng.integers(0, n))
        for t in range(n):
            out[t] = x[i]
            if self.rng.random() < p:
                i = int(self.rng.integers(0, n))
            else:
                i = (i + 1) % n
        return out

    def propose_densification(self, overlap_matrix, centers):
        """Suggest additional window centers to bridge low-overlap neighbors.

        Uses adjacent pairs along sorted centers. For each pair with overlap
        below the configured minimum (default 0.10), propose the midpoint.
        Returns dict with pairs_below_threshold and suggested_centers.
        """
        # Prefer local analysis config; fallback to global defaults
        qc = (self.config.get('qc') or {})
        if not qc:
            cfg = _load_global_pmf_config()
            qc = ((cfg.get('pmf') or {}).get('qc') or {})
        min_ol = float(qc.get('min_neighbor_overlap', 0.10))
        centers = np.asarray(centers, dtype=float)
        order = np.argsort(centers)
        cs = centers[order]
        M = overlap_matrix[np.ix_(order, order)]

        pairs = []
        suggestions = []
        for i in range(len(cs)-1):
            a = float(cs[i]); b = float(cs[i+1])
            ol = float(M[i, i+1])
            if not np.isfinite(ol):
                continue
            if ol < min_ol:
                pairs.append({'z1': a, 'z2': b, 'overlap': ol})
                mid = 0.5 * (a + b)
                suggestions.append(mid)
        # Deduplicate suggestions
        if suggestions:
            suggestions = sorted({float(f"{s:.3f}") for s in suggestions}, reverse=True)
        return {'min_neighbor_overlap': min_ol,
                'pairs_below_threshold': pairs,
                'suggested_centers': suggestions}
    
    def plot_pmf(self, z, pmf, pmf_std, features, output_file, pmf_lower=None, pmf_upper=None, ci_mask=None, centers_for_ticks=None, coverage_counts=None):
        """Plot PMF profile with features and optional confidence interval.

        ci_mask: optional boolean array same length as z; if provided, CI is only
                 drawn where mask is True (helps hide edge bins with poor support).
        centers_for_ticks: optional list of window centers to draw small rug ticks.
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot PMF (hide unsupported bins if ci_mask provided)
        pmf_line = pmf.copy()
        if ci_mask is not None and isinstance(ci_mask, np.ndarray) and ci_mask.shape == pmf.shape:
            pmf_line = pmf_line.copy()
            pmf_line[~ci_mask] = np.nan
        ax.plot(z, pmf_line, 'b-', linewidth=2, label='PMF')
        # CI drawing with optional mask to avoid edge artifacts
        if pmf_lower is not None and pmf_upper is not None and len(pmf_lower) == len(pmf):
            if ci_mask is None or (isinstance(ci_mask, np.ndarray) and ci_mask.dtype==bool and ci_mask.shape==pmf.shape):
                # draw in contiguous segments where mask True
                mask = ci_mask if ci_mask is not None else np.ones_like(pmf, dtype=bool)
                # robust segment detection via diff on padded mask
                m = mask.astype(np.int8)
                diff = np.diff(np.concatenate(([0], m, [0])))
                starts = np.where(diff == 1)[0]
                ends = np.where(diff == -1)[0] - 1
                for s, e in zip(starts, ends):
                    ax.fill_between(z[s:e+1], pmf_lower[s:e+1], pmf_upper[s:e+1], color='b', alpha=0.2)
                ax.plot([], [], color='b', alpha=0.2, label='95% CI')
            else:
                ax.fill_between(z, pmf_lower, pmf_upper, color='b', alpha=0.2, label='95% CI')
        elif pmf_std is not None and len(pmf_std) == len(pmf):
            ax.fill_between(z, pmf - pmf_std, pmf + pmf_std, alpha=0.3)
        
        # Mark features
        if features.get('z_ads') is not None:
            ax.axvline(features['z_ads'], color='g', linestyle='--', alpha=0.5, 
                      label=f"Adsorbed: {features['delta_g_ads']:.1f} kJ/mol")
        if features.get('z_insert') is not None:
            ax.axvline(features['z_insert'], color='r', linestyle='--', alpha=0.5,
                      label=f"Inserted: {features['delta_g_insert']:.1f} kJ/mol")
        if features.get('z_barrier') is not None:
            ax.axvline(features['z_barrier'], color='orange', linestyle='--', alpha=0.5,
                      label=f"ΔG‡ vs ads: {features['delta_g_barrier']:.1f} kJ/mol")
        
        # Add membrane region
        ax.axvspan(-1.5, 1.5, alpha=0.1, color='gray', label='Membrane')
        
        ax.set_xlabel('z-coordinate (nm)', fontsize=12)
        ax.set_ylabel('PMF (kJ/mol)', fontsize=12)
        ax.set_title('Potential of Mean Force Profile', fontsize=14)
        # Window center rug ticks at bottom
        ymin, ymax = ax.get_ylim()
        if centers_for_ticks is not None and len(centers_for_ticks) > 0:
            h = 0.04 * (ymax - ymin)
            for c in centers_for_ticks:
                ax.vlines(c, ymin, ymin + h, colors='k', alpha=0.25, linewidth=1)
        # Coverage bar under PMF (normalized to small band at bottom)
        if coverage_counts is not None and len(coverage_counts) == len(z):
            c = np.asarray(coverage_counts, dtype=float)
            if np.any(c > 0):
                c = c / np.max(c)
                band = 0.05 * (ymax - ymin)
                ax.fill_between(z, ymin, ymin + band * c, color='k', alpha=0.15, step='mid')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_pmf_with_histograms(self, z, pmf, pmf_std, features, output_file, z_edges, hist_per_window, centers,
                                 pmf_lower=None, pmf_upper=None, ci_mask=None, coverage_counts=None):
        """Two-panel figure: top PMF (with CI and coverage bar), bottom per-window z-histograms."""
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1.5])
        ax_top = fig.add_subplot(gs[0])
        # Top panel: PMF (mask unsupported bins if ci_mask provided)
        pmf_line = pmf.copy()
        if ci_mask is not None and isinstance(ci_mask, np.ndarray) and ci_mask.shape == pmf.shape:
            pmf_line = pmf_line.copy()
            pmf_line[~ci_mask] = np.nan
        ax_top.plot(z, pmf_line, 'b-', linewidth=2)
        if pmf_lower is not None and pmf_upper is not None:
            # draw masked CI
            if ci_mask is None:
                ax_top.fill_between(z, pmf_lower, pmf_upper, color='b', alpha=0.2)
            else:
                mask = ci_mask
                m = mask.astype(np.int8)
                diff = np.diff(np.concatenate(([0], m, [0])))
                starts = np.where(diff == 1)[0]
                ends = np.where(diff == -1)[0] - 1
                for s, e in zip(starts, ends):
                    ax_top.fill_between(z[s:e+1], pmf_lower[s:e+1], pmf_upper[s:e+1], color='b', alpha=0.2)
        # Features
        if features.get('z_ads') is not None:
            ax_top.axvline(features['z_ads'], color='g', linestyle='--', alpha=0.5,
                           label=f"Adsorbed: {features['delta_g_ads']:.1f} kJ/mol")
        if features.get('z_insert') is not None:
            ax_top.axvline(features['z_insert'], color='r', linestyle='--', alpha=0.5,
                           label=f"Inserted: {features['delta_g_insert']:.1f} kJ/mol")
        if features.get('z_barrier') is not None:
            ax_top.axvline(features['z_barrier'], color='orange', linestyle='--', alpha=0.5,
                           label=f"ΔG‡ vs ads: {features['delta_g_barrier']:.1f} kJ/mol")
        ax_top.axvspan(-1.5, 1.5, alpha=0.1, color='gray', label='Membrane')
        ax_top.set_ylabel('PMF (kJ/mol)')
        ax_top.set_title('PMF with coverage and window positions')
        # Coverage bar
        ymin, ymax = ax_top.get_ylim()
        if coverage_counts is not None and len(coverage_counts) == len(z):
            c = np.asarray(coverage_counts, dtype=float)
            if np.any(c > 0):
                c = c / np.max(c)
                band = 0.05 * (ymax - ymin)
                ax_top.fill_between(z, ymin, ymin + band * c, color='k', alpha=0.15, step='mid')
        # Rug ticks
        if centers is not None and len(centers) > 0:
            h = 0.04 * (ymax - ymin)
            for cc in centers:
                ax_top.vlines(cc, ymin, ymin + h, colors='k', alpha=0.25, linewidth=1)
        ax_top.legend(loc='best')
        ax_top.grid(True, alpha=0.3)

        # Bottom panel: per-window histograms
        ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2.0
        # Normalize each histogram by total counts (per-bin probabilities)
        for h, _c in zip(hist_per_window, centers):
            total = float(np.sum(h))
            y = (h / total) if total > 0 else h
            ax_bottom.plot(z_centers, y, color='tab:gray', alpha=0.35, linewidth=1)
        ax_bottom.set_xlabel('z-coordinate (nm)')
        ax_bottom.set_ylabel('Density (a.u.)')
        ax_bottom.set_title('Per-window z-histograms (normalized)')
        ax_bottom.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    def plot_overlap_heatmap(self, overlap_matrix, centers, output_file):
        """Plot window overlap heatmap"""
        fig, ax = plt.subplots(figsize=(10, 8))
        # Sort by center for a monotonic diagonal
        centers = np.asarray(centers)
        order = np.argsort(centers)
        M = overlap_matrix[np.ix_(order, order)]
        centers_sorted = centers[order]
        im = ax.imshow(M, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Overlap', fontsize=12)
        
        # Set ticks
        n_windows = len(centers_sorted)
        tick_positions = np.arange(0, n_windows, max(1, n_windows // 10))
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)
        ax.set_xticklabels([f"{centers_sorted[i]:.2f}" for i in tick_positions], rotation=45)
        ax.set_yticklabels([f"{centers_sorted[i]:.2f}" for i in tick_positions])
        
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
        # ESS diagnostics
        print("\nEstimating ESS per window...")
        g_k, ess_k = self.estimate_ess_per_window(data)
        print(f"ESS (mean/min): {ess_k.mean():.1f} / {ess_k.min():.1f}")
        
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
            z_grid, pmf, pmf_std, pmf_lower, pmf_upper, pmf_bootstrap_full = self.calculate_mbar(data, centers, do_bootstrap, n_bootstrap)
        else:
            z_grid, pmf, pmf_std, pmf_lower, pmf_upper, pmf_bootstrap_full = self.calculate_wham(data, centers, do_bootstrap, n_bootstrap)
        
        # Extract features (bulk using supported bins only)
        print("\nExtracting PMF features...")
        # Build counts on the same z grid for support-aware bulk detection
        try:
            dz_feat = float(np.median(np.diff(z_grid))) if len(z_grid) > 1 else self.grid_spacing
            z_edges_feat = np.concatenate(([z_grid[0] - dz_feat/2], (z_grid[:-1] + z_grid[1:]) / 2, [z_grid[-1] + dz_feat/2]))
            counts_feat = np.zeros_like(z_grid, dtype=int)
            for d in data:
                ctmp, _ = np.histogram(d, bins=z_edges_feat)
                counts_feat += ctmp
        except Exception:
            counts_feat = None
        features = self.extract_features(z_grid, pmf, pmf_std, counts=counts_feat)
        
        # Convergence check removed (redundant and expensive)
        convergence_value = None
        
        # Calculate overlap
        print("Calculating window overlap...")
        overlap_matrix = self.calculate_overlap_matrix(data, centers)
        # Ergodicity block stability
        block_R = self.estimate_block_stability(data, n_blocks=5)
        erg_pass = np.isfinite(block_R) & (block_R <= self.ergodicity_block_r_max)
        erg_pass_fraction = float(np.mean(erg_pass)) if erg_pass.size else None
        # JS divergence across adjacent windows
        js_pairs, js_summary = self.js_divergence_adjacent(data, centers, bins=50)
        
        # Generate plots (can be disabled via config)
        plots_dir = self.pmf_dir / "analysis_plots"
        if not self.config.get('no_plots', False):
            print("\nGenerating plots...")
            plots_dir.mkdir(exist_ok=True)
            ci_mask = None  # ensure defined for fallback path
            try:
                # Build coverage mask for plotting CIs only where supported
                dz = float(np.median(np.diff(z_grid))) if len(z_grid) > 1 else self.grid_spacing
                z_edges_plot = np.concatenate(([z_grid[0]-dz/2], (z_grid[:-1]+z_grid[1:])/2, [z_grid[-1]+dz/2]))
                counts = np.zeros_like(z_grid, dtype=int)
                h_per_window = []
                for d in data:
                    c, _ = np.histogram(d, bins=z_edges_plot)
                    counts += c
                    h_per_window.append(c.astype(float) / max(c.sum(), 1))
                min_count = int(((self.config.get('plot') or {}).get('min_bin_count_for_ci') or 25))
                ci_mask = counts >= min_count
                # Optional: refine CI mask using bootstrap coverage/width if available
                if pmf_bootstrap_full is not None and np.size(pmf_bootstrap_full) > 0:
                    finite_frac = np.mean(np.isfinite(pmf_bootstrap_full), axis=0)
                    ci_width = (np.nanpercentile(pmf_bootstrap_full, 97.5, axis=0) -
                                np.nanpercentile(pmf_bootstrap_full,  2.5, axis=0))
                    plot_cfg = (self.config.get('plot') or {})
                    min_frac = float(plot_cfg.get('ci_min_bootstrap_fraction', 0.90))
                    max_band = float(plot_cfg.get('ci_max_band_kj', 200.0))
                    ci_mask_extra = (finite_frac >= min_frac) & np.isfinite(ci_width) & (ci_width <= max_band)
                    if ci_mask is None:
                        ci_mask = ci_mask_extra
                    else:
                        ci_mask = ci_mask & ci_mask_extra
                # PMF plot with coverage + hist panel
                self.plot_pmf_with_histograms(
                    z_grid, pmf, pmf_std, features,
                    plots_dir / "pmf_profile.png",
                    z_edges_plot, h_per_window, centers,
                    pmf_lower=pmf_lower, pmf_upper=pmf_upper, ci_mask=ci_mask,
                    coverage_counts=counts
                )
            except Exception:
                # Fallback: simple PMF plot
                self.plot_pmf(
                    z_grid, pmf, pmf_std, features,
                    plots_dir / "pmf_profile.png",
                    pmf_lower=pmf_lower, pmf_upper=pmf_upper, ci_mask=ci_mask
                )
            # Overlap heatmap always attempted
            self.plot_overlap_heatmap(overlap_matrix, centers, plots_dir / "overlap_matrix.png")
        
        # Summarize overlap
        overlap_summary = self.summarize_overlap(overlap_matrix, centers)
        # Densification suggestions
        densify = self.propose_densification(overlap_matrix, centers)

        # Strict consistency violations: high JS and low overlap for adjacent pairs
        try:
            cfg_global = _load_global_pmf_config()
            qc_global = ((cfg_global.get('pmf') or {}).get('qc') or {})
            target_ol = float(qc_global.get('target_overlap', 0.20))
            js_max = float((self.config.get('qc') or {}).get('js_divergence_max', 0.5))
            centers_arr = np.asarray(centers, dtype=float)
            order = np.argsort(centers_arr)
            consistency_violations = []
            for idx_adj, (i, j) in enumerate(zip(order[:-1], order[1:])):
                js_val = js_pairs[idx_adj]['js'] if idx_adj < len(js_pairs) else None
                ol_val = float(overlap_matrix[i, j])
                if js_val is not None and (js_val > js_max) and (ol_val < target_ol):
                    consistency_violations.append({'z1': float(centers_arr[i]), 'z2': float(centers_arr[j]), 'js': float(js_val), 'overlap': float(ol_val)})
        except Exception:
            consistency_violations = []

        # Adsorption metrics and HC50 estimate
        # If we have bootstrap PMFs, prefer CI from Kp distribution over per-bin envelope
        adsorption = self.compute_adsorption_metrics(z_grid, pmf, pmf_std, pmf_lower, pmf_upper, counts=counts_feat)
        if pmf_bootstrap_full is not None:
            mask = (z_grid >= adsorption['z_low']) & (z_grid <= adsorption['z_high'])
            if np.any(mask):
                # Compute Kp for each bootstrap replicate with stabilization
                x = -self.beta * pmf_bootstrap_full[:, mask]
                x0 = np.max(x, axis=1, keepdims=True)
                kp_scaled = integrate.trapezoid(np.exp(x - x0), z_grid[mask], axis=1)
                dz_len = float(integrate.trapezoid(np.ones_like(z_grid[mask]), z_grid[mask]))
                kp_samples = kp_scaled * np.exp(x0.flatten()) - dz_len
                # Apply bilayer factor if configured
                bilayer_factor = float(((self.config or {}).get('adsorption') or {}).get('bilayer_factor', 2.0))
                kp_samples = kp_samples * bilayer_factor
                kp_samples = kp_samples[np.isfinite(kp_samples) & (kp_samples > 0)]
                if kp_samples.size:
                    kp_lo, kp_hi = float(np.percentile(kp_samples, 2.5)), float(np.percentile(kp_samples, 97.5))
                    adsorption['kp_ads_nm_ci'] = [kp_lo, kp_hi]
                    # Convert to HC50 CI bands for both L/P* extremes
                    area_per_lipid = adsorption['area_per_lipid_nm2']
                    lp0, lp1 = adsorption['lp_star_range']
                    conv = 0.602214
                    gamma_min = 1.0 / (area_per_lipid * lp1)
                    gamma_max = 1.0 / (area_per_lipid * lp0)
                    hc_min_low = float(gamma_min / (conv * kp_hi))
                    hc_min_high = float(gamma_min / (conv * kp_lo))
                    hc_max_low = float(gamma_max / (conv * kp_hi))
                    hc_max_high = float(gamma_max / (conv * kp_lo))
                    adsorption['hc50_uM_ci'] = {
                        'lp_star_min_uM': [hc_min_low * 1e6, hc_min_high * 1e6],
                        'lp_star_max_uM': [hc_max_low * 1e6, hc_max_high * 1e6]
                    }

        # Save results
        pmf_data = {
            'z': z_grid.tolist(),
            'pmf': pmf.tolist(),
            'pmf_std': pmf_std.tolist()
        }
        # Add coverage counts for transparency (total counts per z bin)
        try:
            dz = float(np.median(np.diff(z_grid))) if len(z_grid) > 1 else self.grid_spacing
            z_edges_plot = np.concatenate(([z_grid[0]-dz/2], (z_grid[:-1]+z_grid[1:])/2, [z_grid[-1]+dz/2]))
            counts = np.zeros_like(z_grid, dtype=int)
            for d in data:
                c, _ = np.histogram(d, bins=z_edges_plot)
                counts += c
            pmf_data['bin_counts'] = counts.tolist()
        except Exception:
            pass
        if pmf_lower is not None and pmf_upper is not None:
            pmf_data['pmf_lower'] = pmf_lower.tolist()
            pmf_data['pmf_upper'] = pmf_upper.tolist()
        
        # Compose ESS per window list for YAML
        ess_list_yaml = [
            {'center': float(c), 'g': float(g), 'ess': float(e)}
            for c, g, e in zip(centers, g_k, ess_k)
        ]

        results = {
            'method': method,
            'n_windows': len(data),
            'window_centers': centers.tolist(),
            'features': features,
            'pmf_data': pmf_data,
            'adsorption': adsorption,
            'quality_metrics': {
                # Keep original fields for backward compatibility (all-pairs aggregates)
                'mean_overlap': float(overlap_matrix[np.triu_indices_from(overlap_matrix, k=1)].mean()),
                'min_overlap': float((overlap_matrix[np.triu_indices_from(overlap_matrix, k=1)][overlap_matrix[np.triu_indices_from(overlap_matrix, k=1)] > 0]).min()) if (overlap_matrix[np.triu_indices_from(overlap_matrix, k=1)] > 0).any() else 0.0,
                # Add adjacent-pair metrics to match QC report semantics
                'mean_overlap_adjacent': overlap_summary.get('mean_overlap_adjacent'),
                'min_overlap_adjacent': overlap_summary.get('min_overlap_adjacent'),
                'mean_overlap_all_pairs': overlap_summary.get('mean_overlap_all_pairs'),
                'min_overlap_all_pairs_nonzero': overlap_summary.get('min_overlap_all_pairs_nonzero'),
                'convergence': convergence_value,
                'ess_mean': float(ess_k.mean()),
                'ess_min': float(ess_k.min()),
                'ess_by_window': ess_list_yaml,
                'densify': densify
            },
            'bootstrap': {
                'enabled': do_bootstrap,
                'n_samples': n_bootstrap
            }
        }
        # Add extra QC metrics
        results['quality_metrics']['overlap_summary'] = overlap_summary
        results['quality_metrics']['block_stability_R'] = block_R.tolist()
        results['quality_metrics']['ergodicity_pass_fraction'] = erg_pass_fraction
        if method == 'wham':
            results['quality_metrics']['wham_last_iterations'] = getattr(self, 'wham_last_iters', None)
        results['quality_metrics']['js_divergence_adjacent'] = js_pairs
        results['quality_metrics']['js_divergence_summary'] = js_summary
        results['quality_metrics']['consistency_violations'] = consistency_violations

        # QC pass/fail and strict mode behavior
        qc_cfg = (self.config.get('qc') or {})
        min_erg_frac = float(qc_cfg.get('min_ergodicity_pass_fraction', 0.8))
        min_js_frac = float(qc_cfg.get('min_js_pass_fraction', 0.8))
        # Overlap pass fraction for adjacent pairs
        cfg_global = _load_global_pmf_config()
        qc_global = ((cfg_global.get('pmf') or {}).get('qc') or {})
        min_neighbor_ol = float(qc_global.get('min_neighbor_overlap', 0.10))
        order = np.argsort(np.asarray(centers, dtype=float))
        adj_ols = []
        for i, j in zip(order[:-1], order[1:]):
            adj_ols.append(float(overlap_matrix[i, j]))
        adj_ols = np.array(adj_ols, dtype=float) if adj_ols else np.array([], dtype=float)
        ol_pass_frac = float(np.mean(adj_ols >= min_neighbor_ol)) if adj_ols.size else None
        results['quality_metrics']['overlap_pass_fraction'] = ol_pass_frac
        min_ol_frac = float(qc_cfg.get('min_overlap_pass_fraction', 0.8))
        js_pass_frac = js_summary.get('pass_fraction') if isinstance(js_summary, dict) else None
        qc_pass = True
        reasons = []
        if erg_pass_fraction is not None and erg_pass_fraction < min_erg_frac:
            qc_pass = False
            reasons.append(f"ergodicity_pass_fraction {erg_pass_fraction:.2f} < {min_erg_frac:.2f}")
        if js_pass_frac is not None and js_pass_frac < min_js_frac:
            qc_pass = False
            reasons.append(f"js_pass_fraction {js_pass_frac:.2f} < {min_js_frac:.2f}")
        if ol_pass_frac is not None and ol_pass_frac < min_ol_frac:
            qc_pass = False
            reasons.append(f"overlap_pass_fraction {ol_pass_frac:.2f} < {min_ol_frac:.2f}")
        if consistency_violations and bool(qc_cfg.get('strict_mode', False)):
            qc_pass = False
            reasons.append(f"consistency_violations {len(consistency_violations)}")
        results['quality_metrics']['qc_status'] = {'passed': bool(qc_pass), 'reasons': reasons}
        
        # Save to file
        results_file = self.pmf_dir / "pmf_analysis_results.yaml"
        with open(results_file, 'w') as f:
            yaml.dump(results, f, default_flow_style=False)
        # Also write densification suggestions separately for the runner
        densify_file = self.pmf_dir / "densify_suggestions.yaml"
        with open(densify_file, 'w') as f:
            yaml.dump(densify, f, default_flow_style=False)
        
        # Print summary
        print("\n" + "="*60)
        print("PMF ANALYSIS COMPLETE")
        print("="*60)
        print(f"Method: {method.upper()}")
        print(f"Windows analyzed: {len(data)}")
        if do_bootstrap:
            print(f"Bootstrap samples: {n_bootstrap}")
        print(f"\nKey Features:")
        dg_ads_val = features.get('delta_g_ads')
        if dg_ads_val is None:
            print("  ΔG_ads: N/A")
        else:
            print(f"  ΔG_ads: {dg_ads_val:.2f} kJ/mol")
        if features.get('delta_g_insert') is not None:
            print(f"  ΔG_insert: {features['delta_g_insert']:.2f} kJ/mol")
        if features.get('delta_g_barrier') is not None:
            print(f"  ΔG‡: {features['delta_g_barrier']:.2f} kJ/mol")
        # Print adsorption/HC50 summary
        if adsorption and adsorption.get('kp_ads_nm') is not None:
            z_low = adsorption['z_low']; z_high = adsorption['z_high']
            hc50_uM_rng = adsorption.get('hc50_uM_range')
            print(f"  Kp_ads (z {z_low:.2f}-{z_high:.2f} nm): {adsorption['kp_ads_nm']:.2f} nm")
            if hc50_uM_rng:
                print(f"  HC50 (L/P*={int(adsorption['lp_star_range'][0])}-{int(adsorption['lp_star_range'][1])}): {hc50_uM_rng[0]:.2f}–{hc50_uM_rng[1]:.2f} µM")
        print(f"\nQuality Metrics:")
        # Adjacent neighbors (matches QC report)
        if overlap_summary['mean_overlap_adjacent'] is not None:
            print(f"  Mean overlap (adjacent): {overlap_summary['mean_overlap_adjacent']:.3f}")
            print(f"  Min overlap (adjacent): {overlap_summary['min_overlap_adjacent']:.3f}")
        else:
            print("  Overlap (adjacent): N/A")
        # All-pairs for completeness
        print(f"  Mean overlap (all pairs): {overlap_summary['mean_overlap_all_pairs']:.3f}")
        print(f"  Min overlap >0 (all pairs): {overlap_summary['min_overlap_all_pairs_nonzero']:.3f}")
        print(f"\nOutput files:")
        print(f"  Results: {results_file}")
        if not self.config.get('no_plots', False):
            print(f"  Plots: {plots_dir}")
        # QC summary line and top-5 consistency violations
        qc_status = results['quality_metrics'].get('qc_status', {})
        passed = qc_status.get('passed')
        reasons = qc_status.get('reasons') or []
        print("\nQC Summary:")
        print(f"  Status: {'PASS' if passed else 'FAIL'}")
        ol_frac = results['quality_metrics'].get('overlap_pass_fraction')
        js_frac = results['quality_metrics'].get('js_divergence_summary', {}).get('pass_fraction')
        erg_frac = results['quality_metrics'].get('ergodicity_pass_fraction')
        if ol_frac is not None:
            print(f"  Overlap pass fraction: {ol_frac:.2f}")
        if js_frac is not None:
            print(f"  JS pass fraction: {js_frac:.2f}")
        if erg_frac is not None:
            print(f"  Ergodicity pass fraction: {erg_frac:.2f}")
        if reasons:
            print("  Reasons:")
            for r in reasons:
                print(f"   - {r}")
        viol = results['quality_metrics'].get('consistency_violations') or []
        if viol:
            print("\nTop consistency violations (up to 5):")
            for v in viol[:5]:
                print(f"  Pair [{v['z1']:.2f}, {v['z2']:.2f}]  JS={v['js']:.3f}, overlap={v['overlap']:.3f}")
        
        return results

def _load_global_pmf_config():
    cfg_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
    if cfg_path.exists():
        with open(cfg_path, 'r') as f:
            return yaml.safe_load(f) or {}
    return {}

def _find_replicate_dirs(base: Path):
    """Find replicate/tag dirs containing pmf_metadata.yaml under base.

    Returns a list of Path objects. If base itself has pmf_metadata.yaml, returns [base].
    Otherwise scans children matching 'replicate_*' and 'pmf_*'.
    """
    if (base / "pmf_metadata.yaml").exists():
        return [base]
    dirs = []
    # Prioritize replicate_* then pmf_*
    for pat in ("replicate_*", "pmf_*"):
        for d in sorted(base.glob(pat)):
            if d.is_dir() and (d / "pmf_metadata.yaml").exists():
                dirs.append(d)
    return dirs

def _classify(results: dict):
    """Toxicity classification using HC50 when available; fallback to legacy ΔG rule.

    Rules (configurable in future):
      - If adsorption.hc50_uM_range exists: toxic = (max_uM <= 50)
        and borderline = (50 < max_uM <= 200)
      - Else fallback to legacy: delta_g_ads < -8 and (delta_g_insert <= -3 or barrier < 12)
    """
    # Preferred: HC50-based
    ads = (results or {}).get('adsorption') or {}
    rng = ads.get('hc50_uM_range')
    classification = {"toxic": None, "basis": None}
    if isinstance(rng, (list, tuple)) and len(rng) == 2 and all(isinstance(x, (int, float)) for x in rng):
        max_uM = float(max(rng))
        classification["basis"] = "hc50"
        classification["toxic"] = bool(max_uM <= 50.0)
        classification["severity"] = (
            "toxic" if max_uM <= 50.0 else ("borderline" if max_uM <= 200.0 else "non_toxic")
        )
        classification["hc50_uM"] = rng
        return classification

    # Fallback legacy rule
    features = (results or {}).get('features') or {}
    thr = {"ads": -8.0, "insert": -3.0, "barrier": 12.0}
    dg_ads = features.get('delta_g_ads')
    dg_ins = features.get('delta_g_insert')
    dg_bar = features.get('delta_g_barrier')
    if dg_ads is None:
        return {"toxic": None, "basis": "insufficient_features"}
    cond = (dg_ads < thr["ads"]) and ((dg_ins is not None and dg_ins <= thr["insert"]) or (dg_bar is not None and dg_bar < thr["barrier"]))
    return {"toxic": bool(cond), "basis": "legacy_features"}

def _aggregate_features(rep_features: list):
    """Aggregate features across replicates with simple mean.

    rep_features: list of {features: {...}, path: Path}
    Returns dict of aggregated features with per-key mean.
    """
    keys = set()
    for rf in rep_features:
        keys.update(rf['features'].keys())
    agg = {}
    for k in sorted(keys):
        vals = [rf['features'][k] for rf in rep_features if rf['features'].get(k) is not None]
        if vals:
            agg[k] = float(np.mean(vals))
        else:
            agg[k] = None
    return agg

def main():
    parser = argparse.ArgumentParser(
        description="Analyze PMF from umbrella sampling with MBAR/WHAM"
    )
    parser.add_argument("pmf_dir", help="PMF directory (replicate/tag dir or pmf root)")
    # Default to WHAM for robust, dependency-light analysis
    parser.add_argument("--method", default="wham", choices=["mbar", "wham"], help="Analysis method")
    parser.add_argument("--bootstrap", type=int, default=0, help="Number of bootstrap samples (0 = no bootstrap)")
    parser.add_argument("--no-bootstrap", action="store_true", help="Disable bootstrap uncertainty estimation")
    parser.add_argument("--bootstrap-thin", type=int, default=1, help="Thinning factor applied only to bootstrap resamples (speed-up)")
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    parser.add_argument("--aggregate", action="store_true", help="If multiple replicates found under pmf_dir, aggregate features")
    parser.add_argument("--no-aggregate", action="store_true", help="Force single-replicate mode even if multiples found")
    parser.add_argument("--bootstrap-method", type=str, default=None, choices=["standard","block","stationary"], help="Bootstrap resampling: standard (IID), block (NBB via ESS) or stationary (Politis-Romano)")
    parser.add_argument("--strict-qc", action="store_true", help="Exit non-zero if JS/ergodicity thresholds or strict consistency violations are present")
    # Time filtering (analysis): discard early frames, restrict, thin by time
    parser.add_argument("--begin-ps", type=float, default=None, help="Discard data before this time (ps) in each window")
    parser.add_argument("--end-ps", type=float, default=None, help="Discard data after this time (ps) in each window")
    parser.add_argument("--dt-ps", type=float, default=None, help="Time-based thinning step (ps) applied during analysis")
    parser.add_argument("--discard-fraction", type=float, default=None, help="Discard this fraction from start of each window if --begin-ps not set (0<frac<1)")
    # MBAR tuning knobs
    parser.add_argument("--mbar-max-samples", type=int, default=None, help="Cap total samples for MBAR by global thinning (e.g., 200000)")
    # Removed: --mbar-solver (no effect in current MBAR init)
    parser.add_argument("--mbar-max-iter", type=int, default=None, help="MBAR maximum iterations (e.g., 5000)")
    parser.add_argument("--mbar-rel-tol", type=float, default=None, help="MBAR relative tolerance (e.g., 1e-7)")
    args = parser.parse_args()

    base = Path(args.pmf_dir)
    reps = _find_replicate_dirs(base)

    # If base looks like run_dir/pmf, also try its subdirs
    if not reps and (base / "pmf").exists():
        reps = _find_replicate_dirs(base / "pmf")
        base = base / "pmf"

    if not reps:
        # As a last resort, try richer search up to depth 2
        candidates = list(base.glob("**/pmf_metadata.yaml"))
        reps = [p.parent for p in candidates]

    if not reps:
        print(f"Error: No pmf_metadata.yaml found under {args.pmf_dir}")
        sys.exit(1)

    multi = len(reps) > 1
    do_agg = (args.aggregate or (multi and not args.no_aggregate))

    rep_results = []
    any_qc_fail = False
    for rep_dir in reps:
        print("\n=== Analyzing replicate ===")
        print(rep_dir)
        try:
            analyzer = PMFAnalyzer(rep_dir)
        except FileNotFoundError as e:
            print(f"Skip: {e}")
            continue
        # Override per-run config
        if args.method:
            analyzer.config['method'] = args.method
        if args.no_bootstrap:
            analyzer.config['bootstrap'] = {'enabled': False, 'n_bootstrap': 0}
        elif args.bootstrap > 0:
            analyzer.config['bootstrap'] = {'enabled': True, 'n_bootstrap': args.bootstrap}
        else:
            analyzer.config['bootstrap'] = analyzer.config.get('bootstrap', {'enabled': False, 'n_bootstrap': 0})
        if args.bootstrap_method is not None:
            analyzer.config['bootstrap'] = analyzer.config.get('bootstrap', {}) or {}
            analyzer.config['bootstrap']['method'] = args.bootstrap_method
        # Pass bootstrap thinning
        analyzer.config['bootstrap_thin'] = max(1, int(args.bootstrap_thin))
        # Apply analysis time filters
        tf = analyzer.config.get('time_filter', {}) if isinstance(analyzer.config.get('time_filter', {}), dict) else {}
        if args.begin_ps is not None:
            tf['begin_ps'] = float(args.begin_ps)
        if args.end_ps is not None:
            tf['end_ps'] = float(args.end_ps)
        if args.dt_ps is not None:
            tf['dt_ps'] = float(args.dt_ps)
        if tf:
            analyzer.config['time_filter'] = tf
        if args.discard_fraction is not None:
            analyzer.config['discard_fraction'] = float(args.discard_fraction)
        # Sync CLI time filters into analyzer attributes to ensure they take effect
        tf_eff = analyzer.config.get('time_filter', {}) or {}
        analyzer.begin_ps = tf_eff.get('begin_ps', analyzer.begin_ps)
        analyzer.end_ps   = tf_eff.get('end_ps',   analyzer.end_ps)
        analyzer.dt_ps    = tf_eff.get('dt_ps',    analyzer.dt_ps)
        if 'discard_fraction' in analyzer.config:
            analyzer.discard_fraction = analyzer.config['discard_fraction']
        # Apply MBAR tuning from CLI if provided
        if args.mbar_max_samples or args.mbar_max_iter or args.mbar_rel_tol:
            analyzer.config['mbar'] = analyzer.config.get('mbar', {})
            if args.mbar_max_samples:
                analyzer.config['mbar']['max_samples'] = int(args.mbar_max_samples)
            if args.mbar_max_iter:
                analyzer.config['mbar']['max_iterations'] = int(args.mbar_max_iter)
            if args.mbar_rel_tol:
                analyzer.config['mbar']['relative_tolerance'] = float(args.mbar_rel_tol)

        res = analyzer.generate_analysis_report()
        if not res:
            continue
        # Compute classification (HC50-based if available)
        clf = _classify(res)
        res['classification'] = clf
        # Augment and resave replicate results file
        rep_out = Path(rep_dir) / "pmf_analysis_results.yaml"
        with open(rep_out, 'w') as f:
            yaml.dump(res, f, default_flow_style=False)
        rep_results.append({'path': Path(rep_dir), 'results': res, 'features': res.get('features', {})})
        # Check strict QC status per replicate
        qc_status = ((res.get('quality_metrics') or {}).get('qc_status') or {})
        if qc_status.get('passed') is False:
            any_qc_fail = True

    if not rep_results:
        print("No replicates analyzed successfully.")
        sys.exit(1)

    if do_agg and len(rep_results) >= 2:
        print("\n=== Aggregating replicates ===")
        agg_features = _aggregate_features(rep_results)
        # Replicate consistency check using tolerance from config
        cfg = _load_global_pmf_config()
        tol = ((cfg.get('pmf') or {}).get('qc') or {}).get('replicate_tolerance', 2.0)
        consistency = {}
        for key in ['delta_g_ads', 'delta_g_insert']:
            vals = [rf['features'].get(key) for rf in rep_results if rf['features'].get(key) is not None]
            if len(vals) >= 2:
                spread = float(np.max(vals) - np.min(vals))
                consistency[key] = {'spread': spread, 'tolerance': tol, 'passed': spread <= tol}
        # Aggregate adsorption metrics if available
        hc_ranges = []
        for rr in rep_results:
            ads = (rr.get('results') or {}).get('adsorption') or {}
            rng = ads.get('hc50_uM_range')
            if isinstance(rng, (list, tuple)) and len(rng) == 2 and all(isinstance(x, (int, float)) for x in rng):
                hc_ranges.append([float(rng[0]), float(rng[1])])
        adsorption_agg = {}
        if hc_ranges:
            arr = np.array(hc_ranges, dtype=float)
            # Mean range across replicates
            adsorption_agg['hc50_uM_range'] = [float(arr[:,0].mean()), float(arr[:,1].mean())]
        # Classification on aggregate: prefer HC50 if available
        res_like = {'features': agg_features, 'adsorption': adsorption_agg}
        clf_agg = _classify(res_like)
        aggregate = {
            'n_replicates': len(rep_results),
            'replicates': [str(r['path'].name) for r in rep_results],
            'features_aggregate': agg_features,
            'replicate_consistency': consistency,
            'adsorption_aggregate': adsorption_agg if adsorption_agg else None,
            'classification': clf_agg
        }
        # Save aggregate at base
        out_dir = base
        out_file = out_dir / "pmf_analysis_aggregate.yaml"
        with open(out_file, 'w') as f:
            yaml.dump(aggregate, f, default_flow_style=False)
        print(f"Aggregate results: {out_file}")

    # Strict QC: exit non-zero if requested by CLI or config
    # Config toggle under pmf.analysis.qc.strict_mode
    cfg_global = _load_global_pmf_config()
    strict_cfg = (((cfg_global.get('pmf') or {}).get('analysis') or {}).get('qc') or {})
    strict_mode_cfg = bool(strict_cfg.get('strict_mode', False))
    if args.strict_qc or strict_mode_cfg:
        if any_qc_fail:
            print("\nQC: Strict mode failure (see qc_status.reasons in results). Exiting with code 2.")
            sys.exit(2)
    print("\n✓ Analysis completed successfully")

if __name__ == "__main__":
    main()
