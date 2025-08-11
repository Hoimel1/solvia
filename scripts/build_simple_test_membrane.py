#!/usr/bin/env python3
"""
Build a simple test membrane without cholesterol for quick testing.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SimpleTestMembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use the CG PDB array
        self.peptide_array = Path(f"simulations/membrane_16x/{peptide_id}_16x.pdb")
        
        if not self.peptide_array.exists():
            logger.error(f"Peptide array not found: {self.peptide_array}")
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
        
        # Copy peptide topology
        top_src = Path("simulations/rbc_simple/SOLVIA_1_16x.itp")
        if top_src.exists():
            shutil.copy(top_src, self.output_dir / "SOLVIA_1_16x.itp")
    
    def build_test_membrane(self):
        """Build simple membrane without cholesterol."""
        logger.info("Building simple test membrane (no cholesterol)")
        
        # Simple lipid composition
        n_popc = 150
        n_pope = 100
        n_pops = 30
        
        cmd = [
            "insane",
            "-f", "peptides_16x.pdb",
            "-o", "system.gro",
            "-p", "topol.top",
            "-x", "10",  # Smaller box
            "-y", "10", 
            "-z", "12",
            "-center",
            "-orient",
            "-sol", "W",
            "-salt", "0.15",
            # Only phospholipids, no cholesterol
            "-l", f"POPC:{n_popc}",
            "-l", f"POPE:{n_pope}",
            "-l", f"POPS:{n_pops}",
            "-pbc", "rectangular",
            "-dm", "3"
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            cwd=self.output_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            logger.info("Insane completed successfully!")
            # Save output
            with open(self.output_dir / "insane_output.log", 'w') as f:
                f.write(result.stdout + "\n" + result.stderr)
        else:
            logger.error(f"Insane failed: {result.stderr}")
            return False
            
        # Fix topology
        self.fix_topology()
        
        return True
    
    def fix_topology(self):
        """Create proper Martini 3 topology."""
        logger.info("Creating Martini 3 topology")
        
        # Read insane output
        insane_top = self.output_dir / "topol.top"
        if not insane_top.exists():
            logger.error("Insane topology not found!")
            return
            
        with open(insane_top, 'r') as f:
            insane_lines = f.readlines()
        
        # Extract molecules section
        molecules = []
        in_molecules = False
        for line in insane_lines:
            if "[ molecules ]" in line:
                in_molecules = True
            elif in_molecules and line.strip() and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    mol_count = parts[1]
                    
                    # Fix names
                    if mol_name == "Protein":
                        mol_name = "SOLVIA_1"
                    elif mol_name == "NA+":
                        mol_name = "NA"
                    elif mol_name == "CL-":
                        mol_name = "CL"
                    
                    molecules.append((mol_name, mol_count))
        
        # Create new topology
        topology = f"""; Simple test membrane with 16x {self.peptide_id}
; No cholesterol for easier testing

; Include force field parameters
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "SOLVIA_1_16x.itp"

[ system ]
; Name
Simple test membrane with 16x {self.peptide_id}

[ molecules ]
; Compound        #mols
"""
        
        for mol_name, mol_count in molecules:
            if mol_name != "CHOL":  # Skip cholesterol
                topology += f"{mol_name:<16} {mol_count}\n"
        
        # Write topology
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("Topology created")
        
        # Show summary
        logger.info("System composition:")
        for mol_name, mol_count in molecules:
            if mol_name != "CHOL":
                logger.info(f"  {mol_name}: {mol_count}")
    
    def create_mdp_files(self):
        """Create simple MDP files for testing."""
        
        # Simple energy minimization
        em_mdp = """
; Simple energy minimization
integrator              = steep
nsteps                  = 5000
emtol                   = 100.0

; Output
nstlog                  = 100
nstenergy               = 100

; Cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz

; Electrostatics
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
rvdw                    = 1.1
"""
        
        with open(self.output_dir / "em.mdp", 'w') as f:
            f.write(em_mdp)
            
        # Quick test MD
        md_mdp = """
; Quick test MD
integrator              = md
dt                      = 0.02
nsteps                  = 5000   ; 100 ps

; Output
nstlog                  = 100
nstenergy               = 100
nstxout-compressed      = 100

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 310

; Pressure coupling
pcoupl                  = no

; Cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz

; Electrostatics
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
rvdw                    = 1.1

; Velocity generation
gen_vel                 = yes
gen_temp                = 310
"""
        
        with open(self.output_dir / "md_test.mdp", 'w') as f:
            f.write(md_mdp)
        
        logger.info("MDP files created")
    
    def create_test_script(self):
        """Create a quick test script."""
        
        script = f"""#!/bin/bash
# Quick test of the membrane system

echo "Testing membrane system..."

# 1. Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
if [ $? -eq 0 ]; then
    gmx mdrun -deffnm em -ntmpi 1 -ntomp 4
    echo "Energy minimization complete!"
else
    echo "Error in grompp! Check system.top and system.gro"
    exit 1
fi

# 2. Quick MD test
echo "Step 2: Quick MD test (100 ps)"
gmx grompp -f md_test.mdp -c em.gro -p system.top -o md.tpr -maxwarn 10
if [ $? -eq 0 ]; then
    gmx mdrun -deffnm md -ntmpi 1 -ntomp 4
    echo "Test MD complete!"
else
    echo "Error in MD setup!"
    exit 1
fi

# 3. Visualize
echo "Creating PDB for visualization..."
echo "0" | gmx trjconv -f md.gro -s md.tpr -o final.pdb -pbc mol -ur compact

echo ""
echo "Test complete! Check final.pdb in VMD"
echo "vmd final.pdb"
"""
        
        script_path = self.output_dir / "test_system.sh"
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created test script: {script_path}")
    
    def run(self):
        """Run the complete pipeline."""
        logger.info(f"Building simple test membrane for {self.peptide_id}")
        
        # Build membrane
        if not self.build_test_membrane():
            logger.error("Failed to build membrane")
            return False
            
        # Create MDP files
        self.create_mdp_files()
        
        # Create test script
        self.create_test_script()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Test membrane setup complete in {self.output_dir}")
        logger.info(f"{'='*60}")
        logger.info("\nTo test:")
        logger.info(f"cd {self.output_dir}")
        logger.info("./test_system.sh")
        logger.info(f"{'='*60}\n")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/test_membrane')
    args = parser.parse_args()
    
    builder = SimpleTestMembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
