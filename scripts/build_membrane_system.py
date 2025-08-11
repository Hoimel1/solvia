#!/usr/bin/env python3
"""
Build a complete membrane system with 16 peptides using GROMACS tools.
This replaces INSANE functionality with native GROMACS commands.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy necessary files
        self.prepare_files()
        
    def prepare_files(self):
        """Copy all necessary files to working directory."""
        # Copy peptide array
        src = Path(f"simulations/membrane_16x/{self.peptide_id}_16x.pdb")
        if src.exists():
            shutil.copy(src, self.output_dir / "peptides_16x.pdb")
        
        # Copy force field files
        ff_dir = Path("force_fields/martini3")
        for ff_file in ff_dir.glob("*.itp"):
            shutil.copy(ff_file, self.output_dir / ff_file.name)
        
        # Copy peptide topology
        top_src = Path(f"data/processed/topologies/{self.peptide_id}/molecule_0.itp")
        if top_src.exists():
            shutil.copy(top_src, self.output_dir / "molecule_0.itp")
    
    def convert_to_gro(self):
        """Convert PDB to GRO format."""
        logger.info("Converting peptides to GRO format")
        
        cmd = [
            'gmx', 'editconf',
            '-f', 'peptides_16x.pdb',
            '-o', 'peptides.gro',
            '-box', '12', '12', '15'
        ]
        
        result = subprocess.run(cmd, cwd=self.output_dir, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Conversion failed: {result.stderr}")
            return False
        
        return True
    
    def create_water_box(self):
        """Create a Martini water box."""
        logger.info("Creating Martini water box")
        
        # Create a water box structure (simplified)
        water_gro = """Martini water box
 1000
    1W       W    1   0.300   0.300   0.300
    2W       W    2   0.600   0.300   0.300
    3W       W    3   0.900   0.300   0.300
    4W       W    4   0.300   0.600   0.300
    5W       W    5   0.600   0.600   0.300
  12.00000  12.00000  15.00000
"""
        
        with open(self.output_dir / "water_box.gro", 'w') as f:
            f.write(water_gro)
        
        logger.info("Note: For production, use pre-equilibrated Martini water box")
        return True
    
    def create_topology(self):
        """Create the system topology."""
        logger.info("Creating system topology")
        
        topology = f"""
; Include force field parameters
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "molecule_0.itp"

[ system ]
; Name
16x {self.peptide_id} on POPC membrane in water

[ molecules ]
; Compound        #mols
molecule_0          16
; Will add POPC, W, NA, CL after building
"""
        
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        return True
    
    def create_simple_membrane_coordinates(self):
        """Create simplified membrane coordinates for testing."""
        logger.info("Creating simplified membrane structure")
        
        # This is a placeholder - in production use CHARMM-GUI or INSANE
        
        membrane_info = """
For a complete membrane system, please use one of these methods:

1. CHARMM-GUI Martini Maker (recommended):
   - Go to http://www.charmm-gui.org/
   - Select "Martini Maker"
   - Upload your peptide structure
   - Choose POPC bilayer
   - Set 16 copies of peptide
   - Download and convert to GROMACS format

2. Manual construction:
   - Build POPC bilayer using packmol or other tools
   - Use gmx insert-molecules to add peptides
   - Solvate with Martini water beads
   - Add ions with gmx genion

3. Use MemProtMD or similar services

The current setup provides the framework, but needs proper membrane coordinates.
"""
        
        with open(self.output_dir / "MEMBRANE_BUILDING_INSTRUCTIONS.txt", 'w') as f:
            f.write(membrane_info)
        
        logger.info("Created membrane building instructions")
        return True
    
    def create_production_scripts(self):
        """Create scripts for running the production simulation."""
        
        # Create minimization script
        em_script = """#!/bin/bash
# Energy minimization
echo "Running energy minimization..."
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -nt 8
"""
        
        # Create equilibration script
        eq_script = """#!/bin/bash
# Equilibration protocol for membrane system

echo "Running NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -maxwarn 1
gmx mdrun -v -deffnm nvt -nt 8

echo "Running NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 1
gmx mdrun -v -deffnm npt -nt 8
"""
        
        # Create production script
        prod_script = """#!/bin/bash
# Production MD
echo "Running production MD..."
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -r npt.gro -t npt.cpt -maxwarn 1
gmx mdrun -v -deffnm md -nt 8
"""
        
        # Create analysis script
        analysis_script = """#!/bin/bash
# Basic analysis

echo "Extracting trajectory..."
# Extract protein and membrane for visualization
echo -e "1\\n13\\n" | gmx trjconv -f md.xtc -s md.tpr -o md_protein_membrane.xtc -pbc mol -ur compact

echo "Calculating RMSD..."
echo -e "4\\n4\\n" | gmx rms -f md.xtc -s md.tpr -o rmsd.xvg

echo "Calculating RMSF..."
echo -e "4\\n" | gmx rmsf -f md.xtc -s md.tpr -o rmsf.xvg -res

echo "Creating PDB for VMD..."
echo -e "0\\n" | gmx trjconv -f md.gro -s md.tpr -o final_frame.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo "Load final_frame.pdb and md_protein_membrane.xtc in VMD for visualization"
"""
        
        # Write all scripts
        scripts = {
            "01_minimize.sh": em_script,
            "02_equilibrate.sh": eq_script,
            "03_production.sh": prod_script,
            "04_analyze.sh": analysis_script
        }
        
        for name, content in scripts.items():
            script_path = self.output_dir / name
            with open(script_path, 'w') as f:
                f.write(content)
            os.chmod(script_path, 0o755)
        
        # Create master run script
        master_script = """#!/bin/bash
# Master script to run complete simulation

echo "SOLVIA Membrane Simulation Pipeline"
echo "==================================="
echo "Peptide: """ + self.peptide_id + """
echo "16 copies on POPC membrane"
echo ""
echo "NOTE: This requires a properly built membrane system!"
echo "See MEMBRANE_BUILDING_INSTRUCTIONS.txt"
echo ""
echo "Once you have a complete system.gro file with membrane+peptides+water+ions:"
echo ""
echo "./01_minimize.sh"
echo "./02_equilibrate.sh"
echo "./03_production.sh"
echo "./04_analyze.sh"
"""
        
        with open(self.output_dir / "RUN_PIPELINE.sh", 'w') as f:
            f.write(master_script)
        os.chmod(self.output_dir / "RUN_PIPELINE.sh", 0o755)
        
        logger.info("Created all simulation scripts")
        return True
    
    def copy_mdp_files(self):
        """Copy MDP files from the membrane setup."""
        src_dir = Path(f"simulations/membrane_16x")
        
        for mdp in ['em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
            src = src_dir / mdp
            if src.exists():
                shutil.copy(src, self.output_dir / mdp)
        
        return True
    
    def run(self):
        """Run the complete membrane building pipeline."""
        logger.info(f"Building membrane system for 16x {self.peptide_id}")
        
        # Convert to GRO
        if not self.convert_to_gro():
            return False
        
        # Create water box template
        self.create_water_box()
        
        # Create topology
        self.create_topology()
        
        # Create membrane building instructions
        self.create_simple_membrane_coordinates()
        
        # Copy MDP files
        self.copy_mdp_files()
        
        # Create run scripts
        self.create_production_scripts()
        
        logger.info(f"Setup complete in {self.output_dir}")
        logger.info("Next steps:")
        logger.info("1. Build complete membrane system using CHARMM-GUI or similar")
        logger.info("2. Run the simulation pipeline with RUN_PIPELINE.sh")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/solvia1_16x_membrane')
    args = parser.parse_args()
    
    builder = MembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
