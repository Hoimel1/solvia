#!/usr/bin/env python3
"""
Build physiologically accurate RBC membrane using insane.

RBC membrane composition:
- Outer leaflet: ~45% POPC, 10% PSM, 45% Cholesterol
- Inner leaflet: ~45% POPE, 15% POPS, 40% Cholesterol

Using insane directly as a Python module.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class InsaneRBCMembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use the CG PDB array
        self.peptide_array = Path(f"simulations/membrane_16x/{peptide_id}_16x.pdb")
        
        # Check if we already have the 16x array
        if not self.peptide_array.exists():
            logger.info(f"Creating peptide array first...")
            subprocess.run([
                "python", "scripts/setup_16peptide_membrane.py",
                "--peptide-id", peptide_id,
                "--num-peptides", "16"
            ])
            
        # Copy necessary files
        self.prepare_files()
        
    def prepare_files(self):
        """Copy all necessary files to working directory."""
        # Copy peptide array  
        shutil.copy(self.peptide_array, self.output_dir / "peptides_16x.pdb")
        
        # Copy force field files
        ff_dir = Path("force_fields/martini3")
        for ff_file in ff_dir.glob("*.itp"):
            shutil.copy(ff_file, self.output_dir / ff_file.name)
        
        # Copy peptide topology with correct name
        top_src = Path(f"data/processed/topologies/{self.peptide_id}/{self.peptide_id}.itp")
        if top_src.exists():
            shutil.copy(top_src, self.output_dir / f"{self.peptide_id}.itp")
    
    def build_rbc_membrane(self):
        """Build the asymmetric RBC membrane using insane."""
        logger.info("Building RBC membrane with insane")
        
        # Calculate lipid numbers for 12x12 nm membrane
        # ~216 lipids per leaflet (1.5 lipids/nm²)
        
        # Outer leaflet: 45% POPC, 10% PSM, 45% CHOL
        n_popc_outer = 97
        n_psm_outer = 22
        n_chol_outer = 97
        
        # Inner leaflet: 45% POPE, 15% POPS, 40% CHOL  
        n_pope_inner = 97
        n_pops_inner = 32
        n_chol_inner = 86
        
        # Total lipids
        total_popc = n_popc_outer  # Only in outer
        total_psm = n_psm_outer    # Only in outer
        total_pope = n_pope_inner  # Only in inner
        total_pops = n_pops_inner  # Only in inner
        total_chol = n_chol_outer + n_chol_inner
        
        # Build the insane command for Martini 3
        cmd = [
            "insane",
            "-f", "peptides_16x.pdb",           # Input peptides
            "-o", "system.gro",                 # Output coordinates
            "-p", "system_insane.top",          # Output topology
            "-x", "12",                         # Box X dimension
            "-y", "12",                         # Box Y dimension  
            "-z", "15",                         # Box Z dimension
            "-center",                          # Center the system
            "-sol", "W",                        # Use Martini water
            "-salt", "0.15",                    # 150 mM salt
            # Lipid composition - use M3 prefix for Martini 3
            "-l", f"M3.POPC:{total_popc}",
            "-l", f"M3.POPE:{total_pope}",
            "-l", f"M3.POPS:{total_pops}",
            "-l", f"M3.PSM:{total_psm}",       # PSM is available!
            "-l", f"M3.CHOL:{total_chol}",
            "-pbc", "rectangular",              # PBC type
            # Try to specify asymmetry - upper leaflet composition
            "-u", "M3.POPC:45",                   
            "-u", "M3.PSM:10",
            "-u", "M3.CHOL:45",
            # Membrane parameters
            "-dm", "5"                          # Membrane thickness ~5nm
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
        # Run insane
        result = subprocess.run(
            cmd,
            cwd=self.output_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Insane failed with error:")
            logger.error(f"STDERR: {result.stderr}")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"Return code: {result.returncode}")
            
            # Try without asymmetry flags if it fails
            logger.info("Trying simpler approach...")
            cmd_simple = [
                "insane",
                "-f", "peptides_16x.pdb",
                "-o", "system.gro",
                "-p", "system_insane.top",
                "-x", "12", "-y", "12", "-z", "15",
                "-center",
                "-sol", "W",
                "-salt", "0.15",
                "-l", f"M3.POPC:{total_popc}",
                "-l", f"M3.POPE:{total_pope}",
                "-l", f"M3.POPS:{total_pops}",
                "-l", f"M3.PSM:{total_psm}",
                "-l", f"M3.CHOL:{total_chol}",
                "-pbc", "rectangular",
                "-dm", "5"
            ]
            
            result = subprocess.run(
                cmd_simple,
                cwd=self.output_dir,
                capture_output=True,
                text=True
            )
        
        if result.returncode == 0:
            logger.info("Insane completed successfully!")
            logger.info("Output summary:")
            # Extract key information from output
            for line in result.stdout.split('\n'):
                if 'Adding' in line or 'beads' in line or 'lipids' in line:
                    logger.info(f"  {line.strip()}")
        else:
            logger.error("Insane failed completely")
            logger.error(f"Final STDERR: {result.stderr}")
            logger.error(f"Final STDOUT: {result.stdout}")
            return False
            
        # Fix the topology file
        self.fix_topology()
        
        return True
    
    def fix_topology(self):
        """Fix the topology file generated by insane."""
        logger.info("Fixing topology file for Martini 3")
        
        # Read insane-generated topology
        insane_top = self.output_dir / "system_insane.top"
        if not insane_top.exists():
            logger.error("Insane topology not found!")
            return
            
        with open(insane_top, 'r') as f:
            insane_content = f.read()
        
        # Create proper Martini 3 topology
        topology = f"""; RBC membrane system with 16x {self.peptide_id}
; Generated with insane for Martini 3

; Include force field parameters
#include "martini_v3.0.0.itp"

; Include lipid topologies
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0.0_phospholipids_PC_v2.itp"
#include "martini_v3.0.0_phospholipids_PE_v2.itp" 
#include "martini_v3.0.0_phospholipids_SM_v2.itp"
#include "martini_v3.0_sterols_v1.0.itp"

; Include solvent and ion topologies
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "{self.peptide_id}.itp"

[ system ]
; Name
RBC membrane with 16x {self.peptide_id}

"""
        
        # Extract molecule counts from insane topology
        in_molecules = False
        molecule_lines = []
        
        for line in insane_content.split('\n'):
            if "[ molecules ]" in line:
                in_molecules = True
                continue
            elif in_molecules and line.strip() and not line.startswith(';'):
                # Clean up molecule names
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    mol_count = parts[1]
                    
                    # Fix peptide name
                    if mol_name in ['Protein', 'protein', 'molecule_0']:
                        mol_name = self.peptide_id
                    
                    molecule_lines.append(f"{mol_name:<16} {mol_count}")
        
        # Add molecules section
        topology += "[ molecules ]\n"
        topology += "; Compound        #mols\n"
        for line in molecule_lines:
            topology += line + "\n"
        
        # Write fixed topology
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("Topology fixed and saved as system.top")
        
        # Show summary
        logger.info("System composition:")
        for line in molecule_lines:
            logger.info(f"  {line}")
    
    def create_mdp_files(self):
        """Copy MDP files from RBC membrane directory."""
        src_dir = Path("simulations/rbc_membrane")
        
        for mdp in ['em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
            src = src_dir / mdp
            if src.exists():
                shutil.copy(src, self.output_dir / mdp)
        
        logger.info("MDP files copied")
    
    def create_run_script(self):
        """Create a script to run the complete simulation."""
        
        run_script = f"""#!/bin/bash
#SBATCH --job-name=RBC_{self.peptide_id}
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# Run RBC membrane simulation for {self.peptide_id}

echo "Starting RBC membrane simulation pipeline..."
echo "System: 16x {self.peptide_id} on asymmetric RBC membrane"

# Load modules if needed
# module load gromacs/2021.4-gpu

# 1. Energy minimization
echo "========================================="
echo "Step 1: Energy minimization"
echo "========================================="
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8

# Check if minimization converged
tail -n 20 em.log | grep -E "Steepest|converged"

# 2. NVT equilibration (1 ns)
echo "========================================="
echo "Step 2: NVT equilibration (1 ns)"
echo "========================================="
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -maxwarn 10
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8

# 3. NPT equilibration (5 ns)
echo "========================================="
echo "Step 3: NPT equilibration (5 ns)"  
echo "========================================="
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8

# 4. Production run (100 ns for testing, increase to 500 ns for production)
echo "========================================="
echo "Step 4: Production run (100 ns)"
echo "========================================="
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt -maxwarn 10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu

echo "========================================="
echo "Simulation complete!"
echo "========================================="

# 5. Basic analysis
echo "Step 5: Basic analysis"

# Peptide-membrane distance
echo "1 2" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi -n index.ndx 2>/dev/null || \\
echo "Protein Membrane" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi

# Create visualization trajectory (every 1 ns)
echo "0" | gmx trjconv -f md.xtc -s md.tpr -o viz_full.xtc -pbc mol -ur compact -dt 1000

# Final structure in PDB format
echo "0" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

# System info
echo ""
echo "Analysis files created:"
echo "- peptide_membrane_dist.xvg: Distance between peptides and membrane"
echo "- viz_full.xtc: Trajectory for visualization (every 1 ns)"
echo "- final_structure.pdb: Final structure"
echo ""
echo "To visualize in VMD:"
echo "vmd final_structure.pdb viz_full.xtc"
"""
        
        script_path = self.output_dir / "run_simulation.sh"
        with open(script_path, 'w') as f:
            f.write(run_script)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created run script: {script_path}")
        
        # Also create a quick test script
        test_script = f"""#!/bin/bash
# Quick test run for {self.peptide_id}

echo "Running quick test simulation (only minimization)..."

# Energy minimization only
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
gmx mdrun -deffnm em -ntmpi 1 -ntomp 4 -nsteps 1000

echo "Test complete! Check em.log for results."
"""
        
        test_path = self.output_dir / "test_simulation.sh"
        with open(test_path, 'w') as f:
            f.write(test_script)
        os.chmod(test_path, 0o755)
    
    def run(self):
        """Run the complete pipeline."""
        logger.info(f"Building RBC membrane for {self.peptide_id} using insane")
        logger.info("Membrane composition:")
        logger.info("  Outer: 45% POPC, 10% PSM, 45% CHOL")
        logger.info("  Inner: 45% POPE, 15% POPS, 40% CHOL")
        
        # Build membrane
        if not self.build_rbc_membrane():
            logger.error("Failed to build membrane")
            return False
            
        # Copy MDP files
        self.create_mdp_files()
        
        # Create run script
        self.create_run_script()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"RBC membrane setup complete in {self.output_dir}")
        logger.info(f"{'='*60}")
        logger.info("\nNext steps:")
        logger.info(f"1. cd {self.output_dir}")
        logger.info("2. ./test_simulation.sh  # Quick test")
        logger.info("3. ./run_simulation.sh   # Full simulation")
        logger.info(f"{'='*60}\n")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build RBC membrane with insane")
    parser.add_argument('--peptide-id', default='SOLVIA_1', help='Peptide ID')
    parser.add_argument('--output-dir', default='simulations/rbc_insane', help='Output directory')
    args = parser.parse_args()
    
    builder = InsaneRBCMembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
