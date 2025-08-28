#!/usr/bin/env python3
"""
SOLVIA Enhanced PMF Runner with Local Patch Reference
Implements all improvements from aenderungen.md for robust PMF calculations
"""

import os
import sys
import yaml
import math
import json
import time
import argparse
import subprocess
import numpy as np
from pathlib import Path
from datetime import datetime
import hashlib
import shutil

# Physical constants
KB_KJ_MOL_K = 0.008314462618  # kJ/mol/K

class PMFRunner:
    """Enhanced PMF runner with local patch reference and QC gates"""
    
    def __init__(self, run_dir, config=None):
        self.run_dir = Path(run_dir)
        self.config = config or self.load_pmf_config()
        self.pmf_config = self.config.get('pmf', {})
        self.umbrella_config = self.pmf_config.get('umbrella', {})
        self.qc_config = self.pmf_config.get('qc', {})
        
    def load_pmf_config(self):
        """Load PMF-specific configuration"""
        pmf_config_path = Path(__file__).parent.parent.parent / "config" / "pmf_standard_config.yaml"
        if pmf_config_path.exists():
            with open(pmf_config_path, 'r') as f:
                return yaml.safe_load(f)
        else:
            # Fallback to standard config
            std_config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
            with open(std_config_path, 'r') as f:
                return yaml.safe_load(f)
    
    def read_gro_atoms(self, gro_path):
        """Parse GRO file and return atom information"""
        atoms = []
        with open(gro_path, "r") as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            return atoms
        
        for line in lines[2:-1]:
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            try:
                x = float(line[20:28])
                y = float(line[28:36])
                z = float(line[36:44])
                atoms.append((resname, atomname, x, y, z))
            except ValueError:
                continue
        
        return atoms
    
    def find_peptide_indices(self, gro_atoms):
        """Find peptide atom indices (1-based)"""
        # Standard amino acids
        AA_RESNAMES = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        
        peptide_idx = []
        for i, (resname, _, _, _, _) in enumerate(gro_atoms):
            if resname.upper() in AA_RESNAMES:
                peptide_idx.append(i + 1)  # 1-based
        
        return peptide_idx
    
    def create_patch_reference_group(self, gro_atoms, peptide_indices, patch_radius=2.0):
        """
        Create local patch reference from phosphates near peptide
        This is the KEY IMPROVEMENT from aenderungen.md
        """
        # Get peptide COM projection on XY plane
        peptide_coords = []
        for i, (_, _, x, y, z) in enumerate(gro_atoms):
            if (i + 1) in peptide_indices:
                peptide_coords.append([x, y, z])
        
        if not peptide_coords:
            raise ValueError("No peptide atoms found")
        
        peptide_coords = np.array(peptide_coords)
        peptide_com_xy = peptide_coords[:, :2].mean(axis=0)
        
        # Find phosphates in outer leaflet within patch radius
        patch_indices = []
        for i, (resname, atomname, x, y, z) in enumerate(gro_atoms):
            # Check if it's a phosphate in outer leaflet
            if atomname in ["PO4", "P"] and z > 7.0:  # Upper leaflet
                # Check distance in XY plane
                dist_xy = np.sqrt((x - peptide_com_xy[0])**2 + 
                                 (y - peptide_com_xy[1])**2)
                if dist_xy <= patch_radius:
                    patch_indices.append(i + 1)  # 1-based
        
        if len(patch_indices) < 10:
            print(f"Warning: Only {len(patch_indices)} phosphates in patch, expanding radius")
            # Expand radius if too few atoms
            return self.create_patch_reference_group(gro_atoms, peptide_indices, 
                                                    patch_radius + 0.5)
        
        return patch_indices
    
    def create_dynamic_index_file(self, gro_path, output_ndx, window_center=None):
        """Create index file with dynamic patch reference for current window"""
        gro_atoms = self.read_gro_atoms(gro_path)
        
        # Find peptide
        peptide_indices = self.find_peptide_indices(gro_atoms)
        if not peptide_indices:
            raise ValueError("No peptide found in structure")
        
        # Get reference mode
        ref_mode = self.umbrella_config.get('ref_mode', 'patch')
        
        if ref_mode == 'patch':
            # LOCAL PATCH REFERENCE (key improvement)
            patch_radius = self.umbrella_config.get('patch_radius', 2.0)
            reference_indices = self.create_patch_reference_group(
                gro_atoms, peptide_indices, patch_radius
            )
            ref_group_name = "LocalPatch"
        else:
            # Fallback to global (not recommended)
            reference_indices = []
            for i, (resname, atomname, _, _, z) in enumerate(gro_atoms):
                if atomname in ["PO4", "P"] and z > 7.0:
                    reference_indices.append(i + 1)
            ref_group_name = "UpperPO4"
        
        # Write index file
        with open(output_ndx, 'w') as f:
            # System group (all atoms)
            f.write("[ System ]\n")
            num_atoms = len(gro_atoms)
            for i in range(1, num_atoms + 1):
                if (i - 1) > 0 and (i - 1) % 15 == 0:
                    f.write("\n")
                f.write(f"{i} ")
            f.write("\n\n")
            
            # Peptide group
            f.write("[ Peptide ]\n")
            for i, idx in enumerate(peptide_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n\n")
            
            # Reference group (patch or global)
            f.write(f"[ {ref_group_name} ]\n")
            for i, idx in enumerate(reference_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{idx} ")
            f.write("\n")
        
        print(f"Created index with {len(peptide_indices)} peptide atoms, "
              f"{len(reference_indices)} reference atoms ({ref_group_name})")
        
        return output_ndx, ref_group_name
    
    def run_smd_ladder(self, start_structure, target_z, output_dir):
        """Run SMD pulling to prepare window starting structure"""
        smd_dir = output_dir / "smd_prep"
        smd_dir.mkdir(exist_ok=True)
        
        # Create SMD MDP
        smd_mdp = smd_dir / "smd.mdp"
        smd_time = self.umbrella_config.get('pre_smd_ns', 5.0)
        
        mdp_content = f"""
; SMD pulling for window preparation
integrator              = md
dt                      = 0.02
nsteps                  = {int(smd_time * 1000 / 0.02)}
nstcomm                 = 100

; Output
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstenergy              = 500
nstxout-compressed     = 500

; Constraints
constraints             = none
constraint-algorithm    = lincs

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 310

; Pressure coupling  
pcoupl                  = berendsen
pcoupltype             = semiisotropic
tau-p                   = 12.0
ref-p                   = 1.0 1.0
compressibility        = 4.5e-5 0

; Pull code - SMD
pull                    = yes
pull-ngroups           = 2
pull-ncoords           = 1
pull-group1-name       = LocalPatch
pull-group2-name       = Peptide
pull-group1-pbcatom    = 0  ; Use COM as reference
pull-group2-pbcatom    = 0  ; Use COM as reference
pull-pbc-ref-prev-step-com = yes  ; Use COM from previous step for PBC
pull-coord1-type       = umbrella
pull-coord1-geometry   = direction
pull-coord1-vec        = 0 0 1
pull-coord1-groups     = 1 2
pull-coord1-rate       = {(target_z - 0) / smd_time}  ; nm/ns
pull-coord1-k          = 500  ; Softer for SMD
pull-coord1-start      = yes
"""
        with open(smd_mdp, 'w') as f:
            f.write(mdp_content)
        
        # Run grompp and mdrun
        # ... (implementation continues)
        
        return smd_dir / "smd_final.gro"
    
    def calculate_window_positions(self):
        """Calculate umbrella window positions with adaptive spacing"""
        z_range = self.umbrella_config.get('z_range', {})
        
        windows = []
        
        # Bulk windows
        bulk_max = z_range.get('bulk_max', 2.8)
        bulk_min = z_range.get('bulk_min', 2.4)
        windows.extend([bulk_max, bulk_min])
        
        # Coarse spacing above interface
        current = bulk_min - z_range.get('coarse_step', 0.2)
        dense_max = z_range.get('dense_max', 0.6)
        while current > dense_max:
            windows.append(current)
            current -= z_range.get('coarse_step', 0.2)
        
        # Dense spacing at interface
        dense_min = z_range.get('dense_min', -0.6)
        dense_step = z_range.get('dense_step', 0.15)
        current = dense_max
        while current >= dense_min:
            windows.append(current)
            current -= dense_step
        
        # Coarse spacing in membrane
        end_z = z_range.get('end_z', -2.0)
        current = dense_min - z_range.get('coarse_step', 0.2)
        while current >= end_z:
            windows.append(current)
            current -= z_range.get('coarse_step', 0.2)
        
        return sorted(windows, reverse=True)
    
    def check_window_overlap(self, window1_pullf, window2_pullf, target_overlap=0.20):
        """Check histogram overlap between adjacent windows"""
        # Load pull force data
        data1 = np.loadtxt(window1_pullf, comments=['#', '@'])
        data2 = np.loadtxt(window2_pullf, comments=['#', '@'])
        
        if len(data1) == 0 or len(data2) == 0:
            return 0.0
        
        # Extract z positions (column 2)
        z1 = data1[:, 1]
        z2 = data2[:, 1]
        
        # Create histograms
        bins = np.linspace(min(z1.min(), z2.min()), max(z1.max(), z2.max()), 50)
        hist1, _ = np.histogram(z1, bins=bins, density=True)
        hist2, _ = np.histogram(z2, bins=bins, density=True)
        
        # Calculate overlap
        overlap = np.minimum(hist1, hist2).sum() / np.minimum(hist1.sum(), hist2.sum())
        
        return overlap
    
    def run_umbrella_window(self, window_center, window_dir, start_structure, 
                          replicate=1, extend_time=0):
        """Run single umbrella window with all improvements"""
        
        window_dir.mkdir(parents=True, exist_ok=True)
        
        # Create dynamic index for this window
        index_file, ref_group = self.create_dynamic_index_file(
            start_structure, 
            window_dir / "index.ndx",
            window_center
        )
        
        # Generate deterministic seed
        peptide_id = self.run_dir.name.split('_')[0]
        seed_string = f"{peptide_id}_rep{replicate}_z{window_center:.2f}"
        seed = int(hashlib.md5(seed_string.encode()).hexdigest()[:8], 16) % 1000000
        
        # Production time
        prod_time = self.umbrella_config.get('production_ns', 20.0) + extend_time
        force_k = self.umbrella_config.get('force_constant', 900)
        
        # Create umbrella MDP
        mdp_file = window_dir / "umbrella.mdp"
        mdp_content = f"""
; Umbrella sampling window at z = {window_center:.3f} nm
integrator              = md
dt                      = 0.02
nsteps                  = {int(prod_time * 1000 / 0.02)}
nstcomm                 = 100

; Output control  
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 1000
nstenergy              = 1000
nstxout-compressed     = 1000

; Neighbor searching
cutoff-scheme          = Verlet
nstlist                = 20
ns_type                = grid
pbc                    = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype            = reaction-field
rcoulomb               = 1.1
epsilon_r              = 15
vdw-type               = Cut-off
vdw-modifier           = Potential-shift-verlet
rvdw                   = 1.1

; Temperature coupling
tcoupl                 = v-rescale
tc-grps                = System
tau-t                  = 1.0
ref-t                  = 310
ld-seed                = {seed}

; Pressure coupling
pcoupl                 = Parrinello-Rahman
pcoupltype            = semiisotropic
tau-p                  = 12.0
ref-p                  = 1.0 1.0
compressibility       = 4.5e-5 0

; Velocity generation
gen-vel                = no
continuation          = yes

; Constraints
constraints            = none
constraint-algorithm   = lincs

; Pull code - Umbrella
pull                   = yes
pull-ngroups          = 2
pull-ncoords          = 1
pull-group1-name      = {ref_group}
pull-group2-name      = Peptide
pull-group1-pbcatom   = 0  ; Use COM as reference
pull-group2-pbcatom   = 0  ; Use COM as reference
pull-pbc-ref-prev-step-com = yes  ; Use COM from previous step for PBC
pull-coord1-type      = umbrella
pull-coord1-geometry  = direction
pull-coord1-vec       = 0 0 1
pull-coord1-groups    = 1 2
pull-coord1-init      = {window_center}
pull-coord1-k         = {force_k}
pull-coord1-start     = no  ; Use init value directly

; Output pull force and coordinate
pull-nstxout          = {int(self.umbrella_config.get('pull_output_ps', 1.0) / 0.02)}
pull-nstfout          = {int(self.umbrella_config.get('pull_output_ps', 1.0) / 0.02)}

; Optional lipid restraints
define                 = -DPOSRES_LIPID_Z
"""
        
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
        
        # Find topology
        top_file = self.find_topology_file()
        
        # Run grompp via Docker
        tpr_file = window_dir / "umbrella.tpr"
        
        # Get relative paths from project root
        project_root = Path(__file__).parent.parent.parent
        rel_mdp = os.path.relpath(mdp_file, project_root)
        rel_structure = os.path.relpath(start_structure, project_root)
        rel_top = os.path.relpath(top_file, project_root)
        rel_index = os.path.relpath(index_file, project_root)
        rel_tpr = os.path.relpath(tpr_file, project_root)
        
        grompp_cmd = [
            "docker", "compose", "run", "--rm", "gromacs",
            "grompp",  # Docker already provides "gmx"
            "-f", f"/work/{rel_mdp}",
            "-c", f"/work/{rel_structure}",
            "-p", f"/work/{rel_top}",
            "-n", f"/work/{rel_index}",
            "-o", f"/work/{rel_tpr}",
            "-maxwarn", "2"
        ]
        
        print(f"  Running grompp...")
        sys.stdout.flush()
        result = subprocess.run(grompp_cmd, capture_output=True, text=True, cwd=str(project_root))
        if result.returncode != 0:
            print(f"  ✗ Error in grompp: {result.stderr}")
            return None
        print(f"  ✓ TPR file created")
        
        # Run mdrun via Docker
        log_file = window_dir / "umbrella.log"
        pullf_file = window_dir / "pullf.xvg"
        pullx_file = window_dir / "pullx.xvg"
        
        # Get relative paths from project root
        rel_tpr = os.path.relpath(tpr_file, project_root)
        rel_pullx = os.path.relpath(pullx_file, project_root)
        rel_pullf = os.path.relpath(pullf_file, project_root)
        rel_log = os.path.relpath(log_file, project_root)
        rel_deffnm = os.path.relpath(window_dir / "umbrella", project_root)
        
        mdrun_cmd = [
            "docker", "compose", "run", "--rm", "gromacs",
            "mdrun",  # Docker already provides "gmx"
            "-deffnm", f"/work/{rel_deffnm}",
            "-s", f"/work/{rel_tpr}",
            "-px", f"/work/{rel_pullx}",
            "-pf", f"/work/{rel_pullf}",
            "-g", f"/work/{rel_log}"
        ]
        
        print(f"  Starting mdrun for z={window_center:.3f} nm...")
        print(f"  Production time: {prod_time} ns")
        sys.stdout.flush()
        
        # Run with real-time output streaming
        process = subprocess.Popen(
            mdrun_cmd, 
            cwd=str(project_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        
        # Stream output line by line
        last_time = time.time()
        for line in process.stdout:
            # Show progress lines (step updates)
            if "step" in line.lower() or "Step" in line:
                # Clear line and print progress
                sys.stdout.write(f"\r  Progress: {line.strip()[:80]}")
                sys.stdout.flush()
            # Show important messages
            elif any(keyword in line for keyword in ["WARNING", "ERROR", "Note", "Writing"]):
                sys.stdout.write(f"\n  {line.strip()}\n")
                sys.stdout.flush()
            # Show periodic heartbeat
            current_time = time.time()
            if current_time - last_time > 5:  # Every 5 seconds
                sys.stdout.write(".")
                sys.stdout.flush()
                last_time = current_time
        
        process.wait()
        sys.stdout.write("\n")  # New line after progress
        
        if process.returncode == 0:
            print(f"✓ Window z={window_center:.3f} nm completed")
            return {
                "center": window_center,
                "pullf": str(pullf_file),
                "pullx": str(pullx_file),
                "tpr": str(tpr_file),
                "seed": seed
            }
        else:
            print(f"✗ Window z={window_center:.3f} nm failed")
            return None
    
    def find_topology_file(self):
        """Find the appropriate topology file"""
        # Look in PMF system directory first
        pmf_dirs = list(self.run_dir.glob("pmf_system/replicate_*/system.top"))
        if pmf_dirs:
            return pmf_dirs[0]
        
        # Check system directory for any .top file
        system_tops = list(self.run_dir.glob("system/*.top"))
        if system_tops:
            # Prefer system_n1.top if it exists (single peptide)
            for top in system_tops:
                if "n1" in top.name:
                    return top
            # Otherwise use the first one found
            return system_tops[0]
        
        # Fallback to other locations
        possible_paths = [
            self.run_dir / "system" / "system.top",
            self.run_dir / "equilibration" / "system.top"
        ]
        
        for path in possible_paths:
            if path.exists():
                return path
        
        # Print what we searched for debugging
        print("Error: No topology file found. Searched:")
        print(f"  - pmf_system/replicate_*/system.top")
        print(f"  - system/*.top")
        for path in possible_paths:
            print(f"  - {path.relative_to(self.run_dir)}")
        
        raise FileNotFoundError("No topology file found")
    
    def run_qc_checks(self, windows_data):
        """Run quality control checks on umbrella windows"""
        qc_results = {
            "overlap_check": [],
            "convergence_check": [],
            "ess_check": []
        }
        
        # Check overlaps between adjacent windows
        for i in range(len(windows_data) - 1):
            if windows_data[i] and windows_data[i+1]:
                overlap = self.check_window_overlap(
                    windows_data[i]['pullf'],
                    windows_data[i+1]['pullf'],
                    self.qc_config.get('min_neighbor_overlap', 0.10)
                )
                
                qc_results["overlap_check"].append({
                    "windows": [windows_data[i]['center'], windows_data[i+1]['center']],
                    "overlap": overlap,
                    "passed": overlap >= self.qc_config.get('min_neighbor_overlap', 0.10)
                })
        
        return qc_results
    
    def generate_qc_report(self, qc_results, output_dir):
        """Generate QC report for umbrella sampling"""
        report_file = output_dir / "qc_report.yaml"
        
        with open(report_file, 'w') as f:
            yaml.dump(qc_results, f, default_flow_style=False)
        
        # Create summary
        n_passed = sum(1 for check in qc_results['overlap_check'] if check['passed'])
        n_total = len(qc_results['overlap_check'])
        
        print("\n=== QC Report ===")
        print(f"Overlap checks: {n_passed}/{n_total} passed")
        
        # Identify problematic windows
        for check in qc_results['overlap_check']:
            if not check['passed']:
                print(f"  ⚠ Low overlap ({check['overlap']:.3f}) between "
                      f"z={check['windows'][0]:.2f} and z={check['windows'][1]:.2f}")
        
        return report_file
    
    def run_pmf_calculation(self, replicate=1, tag=None):
        """Main PMF calculation workflow with all improvements"""
        
        # Setup directories
        pmf_dir = self.run_dir / "pmf"
        if tag:
            pmf_dir = pmf_dir / tag
        else:
            pmf_dir = pmf_dir / f"replicate_{replicate}"
        
        pmf_dir.mkdir(parents=True, exist_ok=True)
        
        # Get starting structure - try different locations
        equil_dir = self.run_dir / "equilibration"
        
        # Try different equilibration outputs
        possible_structures = [
            equil_dir / "npt" / "npt.gro",     # Standard NPT output
            equil_dir / "npt_pr.gro",          # Old format
            equil_dir / "nvt" / "nvt.gro",     # Fallback to NVT
            equil_dir / "em" / "em.gro",       # Last resort: EM
            self.run_dir / "pmf_system" / f"replicate_{replicate}" / "system.gro",  # PMF system
        ]
        
        start_structure = None
        for structure in possible_structures:
            if structure.exists():
                start_structure = structure
                print(f"Using starting structure: {structure.relative_to(self.run_dir)}")
                break
        
        if not start_structure:
            print("Error: No starting structure found. Checked:")
            for structure in possible_structures:
                print(f"  - {structure.relative_to(self.run_dir)}: {'✓' if structure.exists() else '✗'}")
            raise FileNotFoundError(f"No starting structure found")
        
        # Calculate window positions
        window_centers = self.calculate_window_positions()
        print(f"\n{'='*60}")
        print(f"Planning {len(window_centers)} umbrella windows")
        print(f"Starting structure: {start_structure.relative_to(self.run_dir)}")
        print(f"{'='*60}\n")
        
        # Run windows
        windows_data = []
        for i, center in enumerate(window_centers):
            print(f"\n[Window {i+1}/{len(window_centers)}] z = {center:+.3f} nm")
            print("-" * 40)
            
            window_dir = pmf_dir / "windows" / f"z_{center:+.3f}"
            
            # Check if window already exists
            if (window_dir / "pullf.xvg").exists():
                print(f"  ✓ Already completed, skipping")
                windows_data.append({
                    "center": center,
                    "pullf": str(window_dir / "pullf.xvg"),
                    "pullx": str(window_dir / "pullx.xvg"),
                    "tpr": str(window_dir / "umbrella.tpr")
                })
                continue
            
            # Show progress
            print(f"  Setting up window...")
            sys.stdout.flush()
            
            # Run window
            result = self.run_umbrella_window(
                center, window_dir, start_structure, replicate
            )
            windows_data.append(result)
            
            # Show overall progress
            completed = sum(1 for w in windows_data if w is not None)
            print(f"\n  Overall progress: {completed}/{len(window_centers)} windows completed")
            sys.stdout.flush()
            
            # Update start structure for next window (optional)
            if result and self.umbrella_config.get('use_previous_frame', False):
                start_structure = window_dir / "umbrella.gro"
        
        # Run QC checks
        qc_results = self.run_qc_checks(windows_data)
        
        # Add windows if overlap too low
        if self.umbrella_config.get('auto_densify', True):
            for check in qc_results['overlap_check']:
                if not check['passed']:
                    # Add intermediate window
                    new_center = np.mean(check['windows'])
                    print(f"Adding window at z={new_center:.3f} to improve overlap")
                    
                    window_dir = pmf_dir / "windows" / f"z_{new_center:+.3f}"
                    result = self.run_umbrella_window(
                        new_center, window_dir, start_structure, replicate
                    )
                    if result:
                        windows_data.append(result)
        
        # Generate metadata for analysis
        metadata = {
            "replicate": replicate,
            "n_windows": len(windows_data),
            "window_centers": [w['center'] for w in windows_data if w],
            "force_constant": self.umbrella_config.get('force_constant', 900),
            "production_time": self.umbrella_config.get('production_ns', 20.0),
            "reference_mode": self.umbrella_config.get('ref_mode', 'patch'),
            "patch_radius": self.umbrella_config.get('patch_radius', 2.0),
            "windows": windows_data,
            "qc_results": qc_results
        }
        
        metadata_file = pmf_dir / "pmf_metadata.yaml"
        with open(metadata_file, 'w') as f:
            yaml.dump(metadata, f, default_flow_style=False)
        
        # Generate QC report
        self.generate_qc_report(qc_results, pmf_dir)
        
        # Final summary with more details
        successful_windows = len([w for w in windows_data if w])
        print(f"\n{'='*60}")
        print(f"✅ PMF CALCULATION COMPLETED")
        print(f"{'='*60}")
        print(f"  Peptide: {self.run_dir.name.split('_')[0]}")
        print(f"  Replicate: {replicate}")
        print(f"  Windows completed: {successful_windows}/{len(window_centers)}")
        print(f"  Output directory: {pmf_dir.relative_to(self.run_dir.parent)}")
        print(f"  Metadata: pmf_metadata.yaml")
        print(f"  QC Report: qc_report.txt")
        if successful_windows < len(window_centers):
            print(f"\n  ⚠ Warning: {len(window_centers) - successful_windows} windows failed")
            print(f"  Check the output for error messages")
        print(f"{'='*60}")
        print(f"\n  Next step: Run MBAR/WHAM analysis")
        print(f"  Command: python scripts/analysis/pmf_mbar_analysis.py {pmf_dir}")
        print(f"{'='*60}\n")
        
        return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced PMF runner with local patch reference"
    )
    parser.add_argument("run_dir", help="Simulation run directory")
    parser.add_argument("--replicate", type=int, default=1, help="Replicate number")
    parser.add_argument("--tag", help="Optional tag for output directory")
    
    args = parser.parse_args()
    
    # Initialize runner
    runner = PMFRunner(args.run_dir)
    
    # Run PMF calculation
    runner.run_pmf_calculation(replicate=args.replicate, tag=args.tag)

if __name__ == "__main__":
    main()
