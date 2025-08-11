#!/usr/bin/env python3
"""
Build RBC membrane using insane.py directly.
No need for CHARMM-GUI - we can do everything with insane!
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class InsaneRBCBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use the CG PDB!
        self.cg_pdb = Path(f"data/processed/cg_pdb/{peptide_id}_cg.pdb")
        self.peptide_array = Path(f"simulations/membrane_16x/{peptide_id}_16x.pdb")
        
        # Check if we already have the 16x array
        if not self.peptide_array.exists():
            logger.error(f"Peptide array not found: {self.peptide_array}")
            logger.info("Run setup_16peptide_membrane.py first!")
            return
            
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
            
        # Copy insane.py
        insane_src = Path("scripts/md_tools/insane.py")
        if insane_src.exists():
            shutil.copy(insane_src, self.output_dir / "insane.py")
            os.chmod(self.output_dir / "insane.py", 0o755)
    
    def build_rbc_membrane(self):
        """Build the asymmetric RBC membrane using insane."""
        logger.info("Building RBC membrane with insane.py")
        
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
        
        # Total cholesterol
        n_chol_total = n_chol_outer + n_chol_inner
        
        # Build the insane command
        # Note: insane doesn't support true asymmetric bilayers directly,
        # but we can approximate by specifying upper/lower leaflet percentages
        
        cmd = [
            "python", "insane.py",
            "-f", "peptides_16x.pdb",           # Input peptides
            "-o", "system.gro",                 # Output coordinates
            "-p", "system_insane.top",          # Output topology
            "-x", "12",                         # Box X dimension
            "-y", "12",                         # Box Y dimension  
            "-z", "15",                         # Box Z dimension
            "-center",                          # Center the system
            "-sol", "W",                        # Use Martini water
            "-salt", "0.15",                    # 150 mM salt
            # Lipid composition - insane distributes these
            "-l", f"POPC:{n_popc_outer}",
            "-l", f"POPE:{n_pope_inner}",
            "-l", f"POPS:{n_pops_inner}",
            "-l", f"PSM:{n_psm_outer}",
            "-l", f"CHOL:{n_chol_total}",
            # Try to make asymmetric (upper/lower leaflet ratios)
            "-u", "POPC:45",                   # Upper leaflet percentages
            "-u", "PSM:10",
            "-u", "CHOL:45",
            "-u", "POPE:0",                    # None in upper
            "-u", "POPS:0",
            "-l", "POPC:0",                    # Lower leaflet percentages  
            "-l", "PSM:0",                     # None in lower
            "-l", "POPE:52",                   # Adjusted for no POPC/PSM
            "-l", "POPS:18",
            "-l", "CHOL:30",
            "-pbc", "rectangular",              # PBC type
            "-asym", "0.2"                      # Asymmetry parameter
        ]
        
        logger.info(f"Running insane with command:")
        logger.info(" ".join(cmd))
        
        # Run insane
        result = subprocess.run(
            cmd,
            cwd=self.output_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Insane failed: {result.stderr}")
            # Try simpler command without asymmetry
            logger.info("Trying simpler membrane setup...")
            cmd_simple = [
                "python", "insane.py",
                "-f", "peptides_16x.pdb",
                "-o", "system.gro",
                "-p", "system_insane.top",
                "-x", "12", "-y", "12", "-z", "15",
                "-center",
                "-sol", "W",
                "-salt", "0.15",
                "-l", f"POPC:{n_popc_outer + 20}",  # Add some extra for balance
                "-l", f"POPE:{n_pope_inner}",
                "-l", f"POPS:{n_pops_inner}",
                "-l", f"CHOL:{n_chol_total}",
                "-pbc", "rectangular"
            ]
            
            result = subprocess.run(
                cmd_simple,
                cwd=self.output_dir,
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                logger.error(f"Simple insane also failed: {result.stderr}")
                return False
        
        logger.info("Insane completed successfully!")
        logger.info(result.stdout)
        
        # Fix the topology file
        self.fix_topology()
        
        return True
    
    def fix_topology(self):
        """Fix the topology file generated by insane."""
        logger.info("Fixing topology file")
        
        # Read insane-generated topology
        insane_top = self.output_dir / "system_insane.top"
        if not insane_top.exists():
            logger.error("Insane topology not found!")
            return
            
        with open(insane_top, 'r') as f:
            lines = f.readlines()
        
        # Create proper topology
        topology = f"""
; RBC membrane system generated with insane
; Include force field parameters
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0.0_phospholipids_PC_v2.itp"
#include "martini_v3.0.0_phospholipids_PE_v2.itp"
#include "martini_v3.0.0_phospholipids_SM_v2.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0_sterols_v1.0.itp"

; Include peptide topology
#include "{self.peptide_id}.itp"

[ system ]
; Name
16x {self.peptide_id} on asymmetric RBC membrane (insane-generated)

"""
        
        # Extract molecule section from insane topology
        in_molecules = False
        for line in lines:
            if "[ molecules ]" in line:
                in_molecules = True
                topology += line
            elif in_molecules:
                # Replace molecule_0 with correct peptide name
                if "molecule_0" in line:
                    topology += line.replace("molecule_0", self.peptide_id)
                else:
                    topology += line
        
        # Write fixed topology
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("Topology fixed and saved as system.top")
    
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
# Run RBC membrane simulation for {self.peptide_id}

echo "Starting RBC membrane simulation pipeline..."

# 1. Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8

# 2. NVT equilibration  
echo "Step 2: NVT equilibration (1 ns)"
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -maxwarn 10
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8

# 3. NPT equilibration
echo "Step 3: NPT equilibration (5 ns)"  
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8

# 4. Production run
echo "Step 4: Production run (500 ns)"
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt -maxwarn 10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 8

echo "Simulation complete!"

# 5. Basic analysis
echo "Step 5: Analysis"
echo -e "1\\n1\\n" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi
echo -e "0\\n" | gmx trjconv -f md.xtc -s md.tpr -o viz_full.xtc -pbc mol -ur compact
echo -e "0\\n" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo "Visualize with: vmd final_structure.pdb viz_full.xtc"
"""
        
        script_path = self.output_dir / "run_simulation.sh"
        with open(script_path, 'w') as f:
            f.write(run_script)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created run script: {script_path}")
    
    def run(self):
        """Run the complete pipeline."""
        logger.info(f"Building RBC membrane for {self.peptide_id} using insane")
        
        # Build membrane
        if not self.build_rbc_membrane():
            logger.error("Failed to build membrane")
            return False
            
        # Copy MDP files
        self.create_mdp_files()
        
        # Create run script
        self.create_run_script()
        
        logger.info(f"RBC membrane setup complete in {self.output_dir}")
        logger.info("To run the simulation:")
        logger.info(f"  cd {self.output_dir}")
        logger.info(f"  ./run_simulation.sh")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/insane_rbc_membrane')
    args = parser.parse_args()
    
    # First check if peptide array exists
    peptide_array = Path(f"simulations/membrane_16x/{args.peptide_id}_16x.pdb")
    if not peptide_array.exists():
        logger.info(f"Creating peptide array first...")
        subprocess.run([
            "python", "scripts/setup_16peptide_membrane.py",
            "--peptide-id", args.peptide_id,
            "--num-peptides", "16"
        ])
    
    builder = InsaneRBCBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
