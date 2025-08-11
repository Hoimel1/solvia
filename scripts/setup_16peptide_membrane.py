#!/usr/bin/env python3
"""
Setup a membrane simulation with 16 copies of a peptide.
Creates a POPC bilayer with peptides distributed on the surface.
"""

import os
import sys
import subprocess
import numpy as np
import shutil
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MultiPeptideMembrane:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Paths
        self.cg_pdb = Path(f"data/processed/cg_pdb/{peptide_id}_cg.pdb")
        self.topology_dir = Path(f"data/processed/topologies/{peptide_id}")
        self.force_field_dir = Path("force_fields/martini3")
        
    def create_peptide_array(self, n_peptides=16):
        """Create 16 copies of the peptide arranged in a 4x4 grid."""
        logger.info(f"Creating array of {n_peptides} peptides")
        
        # Read original peptide
        with open(self.cg_pdb, 'r') as f:
            lines = f.readlines()
        
        # Extract atom lines
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        
        # Calculate peptide dimensions
        coords = []
        for line in atom_lines:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
        
        coords = np.array(coords)
        peptide_size = coords.max(axis=0) - coords.min(axis=0)
        peptide_center = coords.mean(axis=0)
        
        # Create 4x4 grid with 3nm spacing
        spacing = 3.0  # nm
        grid_positions = []
        for i in range(4):
            for j in range(4):
                x_offset = (i - 1.5) * spacing
                y_offset = (j - 1.5) * spacing
                z_offset = 0
                grid_positions.append([x_offset, y_offset, z_offset])
        
        # Write combined PDB
        combined_pdb = self.output_dir / f"{self.peptide_id}_16x.pdb"
        with open(combined_pdb, 'w') as out:
            out.write("TITLE     16 copies of peptide arranged in 4x4 grid\n")
            
            atom_counter = 0
            residue_offset = 0
            
            for idx, offset in enumerate(grid_positions):
                # Get the last residue number from previous peptide
                if idx > 0:
                    residue_offset += len(set([l[22:26].strip() for l in atom_lines]))
                
                for line in atom_lines:
                    if line.startswith('ATOM'):
                        atom_counter += 1
                        
                        # Parse original line
                        atom_name = line[12:16]
                        res_name = line[17:20]
                        chain = chr(65 + idx) if idx < 26 else 'A'  # A-Z
                        res_num = int(line[22:26]) + residue_offset
                        
                        # Update coordinates
                        x = float(line[30:38]) - peptide_center[0] + offset[0] + 6.0
                        y = float(line[38:46]) - peptide_center[1] + offset[1] + 6.0
                        z = float(line[46:54]) - peptide_center[2] + offset[2] + 6.0
                        
                        # Write new line
                        new_line = f"ATOM  {atom_counter:5d} {atom_name} {res_name} {chain}{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
                        out.write(new_line)
            
            out.write("END\n")
        
        logger.info(f"Created {combined_pdb}")
        return combined_pdb
    
    def create_membrane_box(self, peptide_pdb):
        """Create a membrane box with peptides positioned above."""
        logger.info("Creating membrane box")
        
        # Create a large box for membrane + peptides + water
        cmd_box = [
            'gmx', 'editconf',
            '-f', str(peptide_pdb),
            '-o', str(self.output_dir / 'peptides_box.gro'),
            '-box', '12', '12', '15',  # 12x12x15 nm box
            '-c'
        ]
        
        result = subprocess.run(cmd_box, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Box creation failed: {result.stderr}")
            return False
        
        # Create membrane using simple approach
        # In production, would use INSANE or CHARMM-GUI
        self.create_simple_membrane()
        
        return True
    
    def create_simple_membrane(self):
        """Create a simple POPC bilayer."""
        logger.info("Creating POPC bilayer")
        
        # Parameters for POPC bilayer
        box_x, box_y = 12.0, 12.0  # nm
        lipids_per_nm2 = 1.5  # approximate
        n_lipids_per_leaflet = int(box_x * box_y * lipids_per_nm2)
        
        # Create membrane topology entry
        membrane_top = f"""
; Include force field
#include "{self.force_field_dir}/martini_v3.0.0.itp"
#include "{self.force_field_dir}/martini_v3.0.0_phospholipids_v1.itp"
#include "{self.force_field_dir}/martini_v3.0.0_solvents_v1.itp"
#include "{self.force_field_dir}/martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "molecule_0.itp"

[ system ]
16 peptides on POPC membrane

[ molecules ]
molecule_0    16
POPC         {n_lipids_per_leaflet * 2}
W            10000
NA+          50
CL-          50
"""
        
        with open(self.output_dir / 'system.top', 'w') as f:
            f.write(membrane_top)
        
        # Copy peptide topology
        shutil.copy(
            self.topology_dir / 'molecule_0.itp',
            self.output_dir / 'molecule_0.itp'
        )
        
        # Note: In a real setup, we would:
        # 1. Use INSANE to generate actual lipid coordinates
        # 2. Combine peptides with membrane
        # 3. Solvate properly
        # For now, we'll use a simplified approach
        
        logger.info("Membrane topology created")
        logger.info("Note: For production, use INSANE or CHARMM-GUI for proper membrane building")
        
        return True
    
    def create_mdp_files(self):
        """Create MDP files optimized for membrane+peptide system."""
        
        # Energy minimization - gentle for membrane
        em_mdp = """
; Energy minimization for membrane+peptides
integrator              = steep
nsteps                  = 10000
emtol                   = 500.0
emstep                  = 0.01

; Output
nstlog                  = 500
nstenergy               = 500

; Cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
rvdw                    = 1.1
"""
        
        # NVT equilibration - restrain membrane
        nvt_mdp = """
; NVT equilibration for membrane+peptides
integrator              = md
dt                      = 0.02
nsteps                  = 50000  ; 1 ns

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 1000

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein POPC W_ION
tau_t                   = 1.0 1.0 1.0
ref_t                   = 310 310 310

; Pressure
pcoupl                  = no

; Constraints
constraints             = none
constraint-algorithm    = lincs

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
        
        # NPT equilibration - semi-isotropic for membrane
        npt_mdp = """
; NPT equilibration for membrane+peptides
integrator              = md
dt                      = 0.02
nsteps                  = 250000  ; 5 ns

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 1000

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein POPC W_ION
tau_t                   = 1.0 1.0 1.0
ref_t                   = 310 310 310

; Pressure coupling - semi-isotropic for membrane
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 12.0
ref_p                   = 1.0 1.0
compressibility         = 3e-4 3e-4

; Constraints
constraints             = none
constraint-algorithm    = lincs

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
        
        # Production run
        md_mdp = """
; Production MD for membrane+peptides
integrator              = md
dt                      = 0.02
nsteps                  = 5000000  ; 100 ns

; Output
nstlog                  = 5000
nstenergy               = 1000
nstxout-compressed      = 5000
compressed-x-grps       = Protein POPC

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein POPC W_ION
tau_t                   = 1.0 1.0 1.0
ref_t                   = 310 310 310

; Pressure coupling - semi-isotropic
pcoupl                  = parrinello-rahman
pcoupltype              = semiisotropic
tau_p                   = 12.0
ref_p                   = 1.0 1.0
compressibility         = 3e-4 3e-4

; Constraints
constraints             = none
constraint-algorithm    = lincs

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
        
        # Write MDP files
        for name, content in [
            ('em.mdp', em_mdp),
            ('nvt.mdp', nvt_mdp),
            ('npt.mdp', npt_mdp),
            ('md.mdp', md_mdp)
        ]:
            with open(self.output_dir / name, 'w') as f:
                f.write(content)
        
        logger.info("MDP files created")
        return True
    
    def create_run_script(self):
        """Create a script to run the simulation."""
        
        script = f"""#!/bin/bash
# Run membrane simulation with 16 peptides

echo "Starting membrane simulation with 16x {self.peptide_id}"
date

# Note: Full pipeline requires:
# 1. INSANE or packmol to build proper membrane
# 2. Combining peptides with membrane
# 3. Solvation with Martini water
# 4. Ion addition

# For now, we have a simplified setup
echo "This is a simplified setup. For production:"
echo "1. Use INSANE to build POPC bilayer"
echo "2. Position peptides using gmx insert-molecules"
echo "3. Solvate with Martini W beads"
echo "4. Add ions with gmx genion"

# Example commands (to be run after proper setup):
# gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr
# gmx mdrun -v -deffnm em

# gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro
# gmx mdrun -v -deffnm nvt

# gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt
# gmx mdrun -v -deffnm npt

# gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -r npt.gro -t npt.cpt
# gmx mdrun -v -deffnm md

echo "Setup complete. Ready for manual membrane building."
"""
        
        script_path = self.output_dir / 'run_membrane_sim.sh'
        with open(script_path, 'w') as f:
            f.write(script)
        
        os.chmod(script_path, 0o755)
        logger.info(f"Run script created: {script_path}")
        
        return True
    
    def setup(self):
        """Run the complete setup."""
        logger.info(f"Setting up 16x {self.peptide_id} on membrane")
        
        # Create peptide array
        peptide_array = self.create_peptide_array(16)
        
        # Create box
        self.create_membrane_box(peptide_array)
        
        # Create MDP files
        self.create_mdp_files()
        
        # Create run script
        self.create_run_script()
        
        logger.info(f"Setup complete in {self.output_dir}")
        logger.info("Next steps:")
        logger.info("1. Use INSANE or CHARMM-GUI to build proper membrane")
        logger.info("2. Combine with peptides")
        logger.info("3. Run simulations")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/membrane_16x')
    args = parser.parse_args()
    
    sim = MultiPeptideMembrane(args.peptide_id, args.output_dir)
    sim.setup()

if __name__ == "__main__":
    main()
