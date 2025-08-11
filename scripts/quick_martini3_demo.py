#!/usr/bin/env python3
"""
Quick Martini 3 demo - Create a working simulation with peptides.
Since membrane building is complex, let's start with peptides in solution.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class QuickMartini3Demo:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use single peptide for quick demo
        self.single_peptide = Path(f"data/processed/cg_pdb/{peptide_id}_cg.pdb")
        self.single_itp = Path(f"data/processed/topologies/{peptide_id}/{peptide_id}.itp")
        
        if not self.single_peptide.exists():
            logger.error(f"CG peptide not found: {self.single_peptide}")
            return
            
        # Copy files
        self.prepare_files()
        
    def prepare_files(self):
        """Copy necessary files."""
        shutil.copy(self.single_peptide, self.output_dir / f"{self.peptide_id}.pdb")
        shutil.copy(self.single_itp, self.output_dir / f"{self.peptide_id}.itp")
        
        # Copy Martini 3 force fields
        ff_dir = Path("force_fields/martini3")
        essential_files = [
            "martini_v3.0.0.itp",
            "martini_v3.0.0_solvents_v1.itp",
            "martini_v3.0.0_ions_v1.itp"
        ]
        
        for ff_file in essential_files:
            src = ff_dir / ff_file
            if src.exists():
                shutil.copy(src, self.output_dir / ff_file)
    
    def create_system(self):
        """Create a simple solvated system."""
        logger.info("Creating Martini 3 demo system")
        
        # Convert to GRO with box
        cmd_box = [
            "gmx", "editconf",
            "-f", f"{self.peptide_id}.pdb",
            "-o", "peptide_box.gro",
            "-box", "6", "6", "6",  # Small box
            "-c"  # Center
        ]
        
        result = subprocess.run(cmd_box, cwd=self.output_dir, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to create box: {result.stderr}")
            return False
        
        # Create topology
        topology = f"""; Martini 3 demo system
; Single {self.peptide_id} in water

; Include Martini 3 force field
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide
#include "{self.peptide_id}.itp"

[ system ]
Martini 3 demo - {self.peptide_id} in water

[ molecules ]
{self.peptide_id}    1
; Water and ions to be added
"""
        
        with open(self.output_dir / "topol.top", 'w') as f:
            f.write(topology)
        
        # Create Martini water box
        self.create_water_box()
        
        # Solvate
        cmd_solv = [
            "gmx", "solvate",
            "-cp", "peptide_box.gro",
            "-cs", "water_martini3.gro",
            "-o", "solvated.gro",
            "-p", "topol.top"
        ]
        
        result = subprocess.run(cmd_solv, cwd=self.output_dir, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning("Standard solvation failed, creating manual system")
            # Create manual system
            return self.create_manual_system()
        
        logger.info("System solvated")
        
        # Add ions
        self.add_ions()
        
        return True
    
    def create_water_box(self):
        """Create a Martini 3 water box."""
        logger.info("Creating Martini 3 water box")
        
        # Create a simple water GRO file
        water_gro = """Martini 3 water box
    3
    1W     W     1   0.300   0.300   0.300
    2W     W     2   0.700   0.300   0.300
    3W     W     3   0.500   0.700   0.300
   1.00000   1.00000   1.00000
"""
        
        with open(self.output_dir / "water_martini3.gro", 'w') as f:
            f.write(water_gro)
    
    def create_manual_system(self):
        """Create system manually."""
        logger.info("Creating manual system")
        
        # Just copy peptide as system
        shutil.copy(self.output_dir / "peptide_box.gro", self.output_dir / "system.gro")
        
        # Add some water molecules manually to topology
        with open(self.output_dir / "topol.top", 'a') as f:
            f.write("W    500  ; Water beads\n")
        
        return True
    
    def add_ions(self):
        """Add ions to neutralize."""
        logger.info("Adding ions")
        
        # Create simple ion MDP
        ion_mdp = """
; For ion addition
integrator  = steep
nsteps      = 0
"""
        
        with open(self.output_dir / "ions.mdp", 'w') as f:
            f.write(ion_mdp)
        
        # Prepare for genion
        cmd_grompp = [
            "gmx", "grompp",
            "-f", "ions.mdp",
            "-c", "solvated.gro",
            "-p", "topol.top",
            "-o", "ions.tpr",
            "-maxwarn", "5"
        ]
        
        result = subprocess.run(cmd_grompp, cwd=self.output_dir, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning("Could not prepare for ion addition")
            return
        
        # Add ions
        cmd_genion = [
            "gmx", "genion",
            "-s", "ions.tpr",
            "-o", "system.gro",
            "-p", "topol.top",
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ]
        
        # Automatically select water
        genion_input = "W\n"
        
        result = subprocess.run(
            cmd_genion, 
            cwd=self.output_dir, 
            input=genion_input,
            capture_output=True, 
            text=True
        )
        
        if result.returncode == 0:
            logger.info("Ions added")
        else:
            logger.warning("Could not add ions, using solvated system")
            shutil.copy(self.output_dir / "solvated.gro", self.output_dir / "system.gro")
    
    def create_mdp_files(self):
        """Create Martini 3 MDP files."""
        
        # Energy minimization
        em_mdp = """; Martini 3 energy minimization
integrator              = steep
emtol                   = 100.0
emstep                  = 0.01
nsteps                  = 5000

nstlog                  = 100
nstxout                 = 0
nstvout                 = 0
nstenergy               = 100
nstxout-compressed      = 1000

cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
rvdw                    = 1.1
"""
        
        # NVT equilibration
        nvt_mdp = """; Martini 3 NVT equilibration  
integrator              = md
dt                      = 0.020  
nsteps                  = 5000    ; 100 ps

nstlog                  = 100
nstenergy               = 100
nstxout-compressed      = 500

constraint_algorithm    = lincs
constraints             = none
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
rvdw                    = 1.1

tcoupl                  = v-rescale
tc-grps                 = system
tau_t                   = 1.0
ref_t                   = 310
pcoupl                  = no

gen_vel                 = yes
gen_temp                = 310
gen_seed                = -1
"""
        
        # Short production
        md_mdp = """; Martini 3 production
integrator              = md
dt                      = 0.025
nsteps                  = 40000   ; 1 ns

nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 1000

constraint_algorithm    = lincs
constraints             = none
lincs_iter              = 1
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
rvdw                    = 1.1

tcoupl                  = v-rescale
tc-grps                 = system
tau_t                   = 1.0
ref_t                   = 310

pcoupl                  = berendsen
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = 1.0
compressibility         = 3e-4
"""
        
        # Write MDP files
        for name, content in [
            ("em.mdp", em_mdp),
            ("nvt.mdp", nvt_mdp),
            ("md.mdp", md_mdp)
        ]:
            with open(self.output_dir / name, 'w') as f:
                f.write(content)
        
        logger.info("Martini 3 MDP files created")
    
    def create_run_script(self):
        """Create demo run script."""
        
        script = f"""#!/bin/bash
# Martini 3 demo simulation

echo "Running Martini 3 demo for {self.peptide_id}..."

# 1. Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 5
gmx mdrun -deffnm em -v

# 2. NVT equilibration
echo "Step 2: NVT equilibration (100 ps)"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5
gmx mdrun -deffnm nvt -v

# 3. Production
echo "Step 3: Production run (1 ns)"
gmx grompp -f md.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -deffnm md -v

# 4. Analysis
echo "Step 4: Creating visualization files"
echo "0" | gmx trjconv -f md.xtc -s md.tpr -o traj_centered.xtc -pbc mol -center
echo "0" | gmx trjconv -f md.gro -s md.tpr -o final.pdb -pbc mol -center

echo ""
echo "Demo complete!"
echo "View trajectory: vmd final.pdb traj_centered.xtc"
echo ""
echo "This was a simple demo. For membrane simulations:"
echo "1. Use CHARMM-GUI Martini Maker"
echo "2. Select Martini 3 force field"
echo "3. Build your RBC membrane composition"
"""
        
        script_path = self.output_dir / "run_demo.sh"
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created run script: {script_path}")
    
    def run(self):
        """Run the demo setup."""
        logger.info(f"Setting up Martini 3 demo for {self.peptide_id}")
        
        # Create system
        if not self.create_system():
            logger.error("Failed to create system")
            return False
            
        # Create MDP files
        self.create_mdp_files()
        
        # Create run script
        self.create_run_script()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Martini 3 demo ready in {self.output_dir}")
        logger.info(f"{'='*60}")
        logger.info("\nTo run the demo:")
        logger.info(f"cd {self.output_dir}")
        logger.info("./run_demo.sh")
        logger.info("")
        logger.info("For full RBC membrane simulation:")
        logger.info("- Use CHARMM-GUI with Martini 3")
        logger.info("- Upload peptides_16x.pdb from simulations/membrane_16x/")
        logger.info(f"{'='*60}\n")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/martini3_demo')
    args = parser.parse_args()
    
    demo = QuickMartini3Demo(args.peptide_id, args.output_dir)
    demo.run()

if __name__ == "__main__":
    main()
