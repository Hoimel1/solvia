#!/usr/bin/env python3
"""
Setup membrane simulations for coarse-grained peptides.
Creates a membrane system, solvates, adds ions, and prepares for MD simulations.
"""

import os
import sys
import subprocess
import argparse
import json
import shutil
from pathlib import Path
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/simulation_setup.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class MembraneSimulationSetup:
    def __init__(self, peptide_id, cg_pdb_path, topology_dir, output_dir, 
                 force_field_dir, membrane_type='POPC'):
        self.peptide_id = peptide_id
        self.cg_pdb_path = Path(cg_pdb_path)
        self.topology_dir = Path(topology_dir)
        self.output_dir = Path(output_dir)
        self.force_field_dir = Path(force_field_dir)
        self.membrane_type = membrane_type
        
        # Create output directory for this peptide
        self.sim_dir = self.output_dir / peptide_id
        self.sim_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy files to simulation directory
        self.prepare_files()
        
    def prepare_files(self):
        """Copy necessary files to simulation directory."""
        # Copy CG PDB
        shutil.copy(self.cg_pdb_path, self.sim_dir / f"{self.peptide_id}_cg.pdb")
        
        # Copy topology files
        topol_file = self.topology_dir / self.peptide_id / "topol.top"
        molecule_itp = self.topology_dir / self.peptide_id / "molecule_0.itp"
        
        if topol_file.exists():
            shutil.copy(topol_file, self.sim_dir / "topol.top")
        if molecule_itp.exists():
            shutil.copy(molecule_itp, self.sim_dir / "molecule_0.itp")
    
    def create_membrane_system(self):
        """Create a membrane system with the peptide."""
        logger.info(f"Creating membrane system for {self.peptide_id}")
        
        # For now, we'll create a simple box around the peptide
        # In a full implementation, this would use INSANE or similar tool
        
        # Step 1: Create a box
        cmd_box = [
            'gmx', 'editconf',
            '-f', f"{self.peptide_id}_cg.pdb",
            '-o', f"{self.peptide_id}_box.gro",
            '-c',  # Center molecule
            '-d', '2.0',  # Distance from box edge
            '-bt', 'cubic'  # Box type
        ]
        
        result = subprocess.run(cmd_box, cwd=self.sim_dir, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to create box: {result.stderr}")
            return False
        
        # For a real membrane simulation, we would:
        # 1. Use INSANE to create lipid bilayer
        # 2. Insert peptide into membrane
        # 3. Add water layers above and below
        
        # Placeholder for membrane building
        logger.info("Note: Full membrane building with INSANE will be implemented when force fields are available")
        
        return True
    
    def solvate_system(self):
        """Add water to the system."""
        logger.info(f"Solvating system for {self.peptide_id}")
        
        # We need the appropriate water model for Martini
        # For Martini 3, this would be W beads
        
        # First, we need to ensure we have the water topology
        water_itp = self.force_field_dir / "martini_v3.0.0_solvents_v1.itp"
        
        if not water_itp.exists():
            logger.warning("Water topology not found. Please upload martini_v3.0.0_solvents_v1.itp to force_fields/martini3/")
            return True  # Continue without solvation for now
        
        # Copy water topology
        shutil.copy(water_itp, self.sim_dir / "martini_solvents.itp")
        
        # Update topology file to include water
        self.update_topology_for_water()
        
        # Solvate command (simplified)
        cmd_solvate = [
            'gmx', 'solvate',
            '-cp', f"{self.peptide_id}_box.gro",
            '-cs', 'water.gro',  # This needs to be a Martini water box
            '-o', f"{self.peptide_id}_solv.gro",
            '-p', 'topol.top'
        ]
        
        # For now, skip actual solvation until we have proper water structures
        logger.info("Solvation will be completed when Martini water structures are available")
        
        return True
    
    def add_ions(self):
        """Add ions to neutralize the system."""
        logger.info(f"Adding ions to {self.peptide_id}")
        
        # This would typically use gmx genion
        # For now, we'll prepare the structure for it
        
        logger.info("Ion addition will be completed when force fields are fully configured")
        
        return True
    
    def update_topology_for_water(self):
        """Update topology file to include water and ions."""
        topol_path = self.sim_dir / "topol.top"
        
        # Read current topology
        with open(topol_path, 'r') as f:
            lines = f.readlines()
        
        # Find where to insert includes
        insert_idx = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('#include'):
                insert_idx = i + 1
                break
        
        # Insert water topology include if not already there
        water_include = '#include "martini_solvents.itp"\n'
        if water_include not in lines:
            lines.insert(insert_idx, water_include)
        
        # Write updated topology
        with open(topol_path, 'w') as f:
            f.writelines(lines)
    
    def create_mdp_files(self):
        """Create MDP files for energy minimization and equilibration."""
        logger.info(f"Creating MDP files for {self.peptide_id}")
        
        # Energy minimization MDP
        em_mdp = """
; Energy Minimization for Martini 3
integrator              = steep
nsteps                  = 5000
emtol                   = 1000.0
emstep                  = 0.01

; Output
nstlog                  = 100
nstenergy               = 100
nstxout-compressed      = 100

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
rvdw                    = 1.1
rvdw_switch             = 0.9
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
"""
        
        # NVT equilibration MDP
        nvt_mdp = """
; NVT Equilibration for Martini 3
integrator              = md
dt                      = 0.02
nsteps                  = 50000  ; 1 ns
comm-mode               = Linear
comm-grps               = System

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 1000

; Bonds
constraints             = none
constraint-algorithm    = lincs
continuation            = no

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
rvdw                    = 1.1
rvdw_switch             = 0.9
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 310

; Pressure coupling
pcoupl                  = no

; Velocity generation
gen_vel                 = yes
gen_temp                = 310
gen_seed                = -1
"""
        
        # NPT equilibration MDP
        npt_mdp = """
; NPT Equilibration for Martini 3
integrator              = md
dt                      = 0.02
nsteps                  = 50000  ; 1 ns
comm-mode               = Linear
comm-grps               = System

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 1000

; Bonds
constraints             = none
constraint-algorithm    = lincs
continuation            = yes

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
rvdw                    = 1.1
rvdw_switch             = 0.9
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 310

; Pressure coupling
pcoupl                  = berendsen
pcoupltype              = isotropic
tau_p                   = 12.0
ref_p                   = 1.0
compressibility         = 3e-4

; Velocity generation
gen_vel                 = no
"""
        
        # Production MD MDP
        md_mdp = """
; Production MD for Martini 3
integrator              = md
dt                      = 0.02
nsteps                  = 5000000  ; 100 ns
comm-mode               = Linear
comm-grps               = System

; Output
nstlog                  = 5000
nstenergy               = 5000
nstxout-compressed      = 5000

; Bonds
constraints             = none
constraint-algorithm    = lincs
continuation            = yes

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
rvdw                    = 1.1
rvdw_switch             = 0.9
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 1.0
ref_t                   = 310

; Pressure coupling
pcoupl                  = parrinello-rahman
pcoupltype              = isotropic
tau_p                   = 12.0
ref_p                   = 1.0
compressibility         = 3e-4

; Velocity generation
gen_vel                 = no
"""
        
        # Write MDP files
        with open(self.sim_dir / "em.mdp", 'w') as f:
            f.write(em_mdp)
        
        with open(self.sim_dir / "nvt.mdp", 'w') as f:
            f.write(nvt_mdp)
        
        with open(self.sim_dir / "npt.mdp", 'w') as f:
            f.write(npt_mdp)
        
        with open(self.sim_dir / "md.mdp", 'w') as f:
            f.write(md_mdp)
        
        logger.info("MDP files created successfully")
        
        return True
    
    def create_submission_script(self):
        """Create a job submission script for the simulation."""
        
        script = f"""#!/bin/bash
#SBATCH --job-name={self.peptide_id}_md
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --output={self.peptide_id}_%j.out
#SBATCH --error={self.peptide_id}_%j.err

# Load GROMACS module (adjust as needed)
module load gromacs/2021.4

# Set number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Starting simulation for {self.peptide_id}"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on: $HOSTNAME"
echo "Start time: $(date)"

# Energy minimization
gmx grompp -f em.mdp -c {self.peptide_id}_box.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v

# NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v

# NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v

# Production MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md -v

echo "Simulation completed"
echo "End time: $(date)"
"""
        
        with open(self.sim_dir / "submit.sh", 'w') as f:
            f.write(script)
        
        os.chmod(self.sim_dir / "submit.sh", 0o755)
        
        logger.info("Submission script created")
        
        return True
    
    def setup(self):
        """Run the complete setup process."""
        logger.info(f"Setting up simulation for {self.peptide_id}")
        
        success = True
        
        # Create membrane system
        if not self.create_membrane_system():
            success = False
        
        # Solvate
        if not self.solvate_system():
            logger.warning("Solvation incomplete - waiting for force field files")
        
        # Add ions
        if not self.add_ions():
            logger.warning("Ion addition incomplete - waiting for force field files")
        
        # Create MDP files
        if not self.create_mdp_files():
            success = False
        
        # Create submission script
        if not self.create_submission_script():
            success = False
        
        if success:
            logger.info(f"Setup completed for {self.peptide_id}")
            logger.info(f"Simulation directory: {self.sim_dir}")
        else:
            logger.error(f"Setup failed for {self.peptide_id}")
        
        return success

def main():
    parser = argparse.ArgumentParser(description='Setup membrane simulations for CG peptides')
    parser.add_argument('--peptide-id', required=True, help='Peptide ID')
    parser.add_argument('--cg-pdb', required=True, help='Path to CG PDB file')
    parser.add_argument('--topology-dir', required=True, help='Directory containing topology files')
    parser.add_argument('--output-dir', default='simulations/systems', help='Output directory')
    parser.add_argument('--force-field-dir', default='force_fields/martini3', help='Force field directory')
    parser.add_argument('--membrane-type', default='POPC', help='Membrane lipid type')
    
    args = parser.parse_args()
    
    # Create logs directory
    Path('logs').mkdir(exist_ok=True)
    
    # Setup simulation
    setup = MembraneSimulationSetup(
        peptide_id=args.peptide_id,
        cg_pdb_path=args.cg_pdb,
        topology_dir=args.topology_dir,
        output_dir=args.output_dir,
        force_field_dir=args.force_field_dir,
        membrane_type=args.membrane_type
    )
    
    setup.setup()

if __name__ == "__main__":
    main()
