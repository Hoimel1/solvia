#!/usr/bin/env python3
"""
Build a Martini 3 compatible membrane system using GROMACS tools directly.
Since insane generates Martini 2 format, we'll use a different approach.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Martini3MembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Paths
        self.peptide_array = Path(f"simulations/membrane_16x/{peptide_id}_16x.pdb")
        self.peptide_itp = Path("simulations/rbc_simple/SOLVIA_1_16x.itp")
        
        if not self.peptide_array.exists():
            logger.error(f"Peptide array not found: {self.peptide_array}")
            return
            
        # Copy necessary files
        self.prepare_files()
        
    def prepare_files(self):
        """Copy all necessary files to working directory."""
        # Copy peptide array  
        shutil.copy(self.peptide_array, self.output_dir / "peptides_16x.pdb")
        
        # Copy peptide ITP
        if self.peptide_itp.exists():
            shutil.copy(self.peptide_itp, self.output_dir / "SOLVIA_1_16x.itp")
        
        # Copy force field files
        ff_dir = Path("force_fields/martini3")
        for ff_file in ff_dir.glob("*.itp"):
            shutil.copy(ff_file, self.output_dir / ff_file.name)
    
    def convert_pdb_to_gro(self):
        """Convert peptide PDB to GRO format."""
        logger.info("Converting peptide PDB to GRO")
        
        cmd = [
            "gmx", "editconf",
            "-f", "peptides_16x.pdb",
            "-o", "peptides.gro",
            "-box", "12", "12", "15",
            "-c"  # Center in box
        ]
        
        result = subprocess.run(
            cmd,
            cwd=self.output_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Failed to convert PDB: {result.stderr}")
            return False
            
        logger.info("Peptides converted and centered in box")
        return True
    
    def create_membrane_coordinates(self):
        """Create simplified membrane coordinates for Martini 3."""
        logger.info("Creating Martini 3 membrane coordinates")
        
        # Box dimensions
        box_x = 12.0
        box_y = 12.0
        box_z = 15.0
        
        # Membrane parameters
        lipids_per_nm2 = 1.5
        membrane_area = box_x * box_y
        n_lipids_per_leaflet = int(membrane_area * lipids_per_nm2 / 2)
        
        # RBC-like composition (simplified without cholesterol for now)
        # Outer leaflet
        n_popc_outer = int(0.55 * n_lipids_per_leaflet)  # 55% instead of 45% (no CHOL)
        n_psm_outer = int(0.12 * n_lipids_per_leaflet)   # 12% instead of 10%
        n_dppc_outer = int(0.33 * n_lipids_per_leaflet)  # Use DPPC as PSM substitute
        
        # Inner leaflet
        n_pope_inner = int(0.60 * n_lipids_per_leaflet)  # 60% instead of 45%
        n_pops_inner = int(0.20 * n_lipids_per_leaflet)  # 20% instead of 15%
        n_popc_inner = int(0.20 * n_lipids_per_leaflet)  # Some POPC in inner too
        
        logger.info(f"Membrane composition (no cholesterol):")
        logger.info(f"  Outer: {n_popc_outer} POPC, {n_dppc_outer} DPPC (PSM substitute)")
        logger.info(f"  Inner: {n_pope_inner} POPE, {n_pops_inner} POPS, {n_popc_inner} POPC")
        
        # Create topology
        self.create_topology(
            n_popc_outer, n_dppc_outer,
            n_pope_inner, n_pops_inner, n_popc_inner
        )
        
        return True
    
    def create_topology(self, n_popc_outer, n_dppc_outer, 
                       n_pope_inner, n_pops_inner, n_popc_inner):
        """Create Martini 3 topology."""
        
        topology = f"""; Martini 3 RBC-like membrane with 16x {self.peptide_id}
; Simplified composition without cholesterol

; Include Martini 3 force field
#include "martini_v3.0.0.itp"

; Include Martini 3 lipids
#include "martini_v3.0.0_phospholipids_v1.itp"

; Include solvents and ions
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "SOLVIA_1_16x.itp"

[ system ]
; Name
Martini 3 RBC membrane with 16x {self.peptide_id}

[ molecules ]
; Compound        #mols
SOLVIA_1         1
; Membrane composition (placeholder - adjust after building)
; Outer leaflet
POPC             {n_popc_outer}
DPPC             {n_dppc_outer}
; Inner leaflet  
POPE             {n_pope_inner}
POPS             {n_pops_inner}
POPC             {n_popc_inner}
; To be added: water and ions
"""
        
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("Martini 3 topology created")
    
    def create_simple_system(self):
        """Create a simple solvated system for testing."""
        logger.info("Creating simple solvated system")
        
        # Just add water around peptides
        cmd_solvate = [
            "gmx", "solvate",
            "-cp", "peptides.gro",
            "-cs", "water.gro",  # Will need Martini water
            "-o", "solvated.gro",
            "-p", "system.top"
        ]
        
        # For now, create a simple box with peptides
        logger.info("System prepared for manual membrane insertion")
        
        # Create instructions
        instructions = """
MEMBRANE BUILDING INSTRUCTIONS FOR MARTINI 3
============================================

Since insane.py generates Martini 2 format, we need alternative approaches:

Option 1: Use CHARMM-GUI Martini Maker
---------------------------------------
1. Go to http://www.charmm-gui.org/
2. Select "Martini Maker"
3. Upload peptides_16x.pdb
4. Select Martini 3 force field
5. Build asymmetric membrane:
   - Outer: 55% POPC, 12% DPPC, 33% CHOL
   - Inner: 45% POPE, 20% POPS, 35% CHOL
6. Download and use with our topology

Option 2: Use MemGen from MemProtMD
------------------------------------
1. Visit http://memprotmd.bioch.ox.ac.uk/
2. Use their membrane builder for Martini 3

Option 3: Manual Construction
-----------------------------
1. Download pre-equilibrated Martini 3 membrane
2. Use gmx commands to insert peptides:
   gmx editconf -f membrane.gro -o membrane_box.gro -c -box 12 12 15
   gmx insert-molecules -f membrane_box.gro -ci peptides.gro -o system.gro

Option 4: Use Polyply (if available)
------------------------------------
polyply gen_coords -p system.top -o system.gro -name POPC POPE POPS

Current Status:
- Peptides are prepared in peptides.gro
- Topology template is in system.top
- Add membrane coordinates manually
"""
        
        with open(self.output_dir / "MEMBRANE_INSTRUCTIONS.txt", 'w') as f:
            f.write(instructions)
        
        return True
    
    def create_test_system_without_membrane(self):
        """Create a test system with just peptides in water."""
        logger.info("Creating test system without membrane")
        
        # Create water box topology
        topology = f"""; Martini 3 test system - peptides in water
; For testing before membrane insertion

; Include Martini 3 force field
#include "martini_v3.0.0.itp"

; Include solvents and ions
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "SOLVIA_1_16x.itp"

[ system ]
; Name
16x {self.peptide_id} in water (test system)

[ molecules ]
; Compound        #mols
SOLVIA_1         1
; Water will be added by gmx solvate
"""
        
        with open(self.output_dir / "test_nowat.top", 'w') as f:
            f.write(topology)
        
        # Create simple MDP for testing
        em_mdp = """
; Martini 3 energy minimization
integrator              = steep
nsteps                  = 5000
emtol                   = 100.0

; Output
nstlog                  = 100
nstenergy               = 100

; Cutoffs - Martini 3 settings
cutoff-scheme           = Verlet
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics - Martini 3
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15    ; Martini 3 
epsilon_rf              = 0
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
rvdw                    = 1.1
"""
        
        with open(self.output_dir / "em_martini3.mdp", 'w') as f:
            f.write(em_mdp)
        
        # Create test script
        script = f"""#!/bin/bash
# Test Martini 3 system

echo "Testing Martini 3 peptide system..."

# 1. Test without water first
echo "Step 1: Testing peptide topology..."
gmx grompp -f em_martini3.mdp -c peptides.gro -p test_nowat.top -o test.tpr -maxwarn 5
if [ $? -eq 0 ]; then
    echo "Topology test passed!"
else
    echo "Topology error - check SOLVIA_1_16x.itp"
    exit 1
fi

echo ""
echo "Peptide system is ready for membrane insertion."
echo "See MEMBRANE_INSTRUCTIONS.txt for next steps."
"""
        
        script_path = self.output_dir / "test_martini3.sh"
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        
        logger.info("Test system created")
        return True
    
    def run(self):
        """Run the pipeline."""
        logger.info(f"Building Martini 3 compatible system for {self.peptide_id}")
        
        # Convert PDB to GRO
        if not self.convert_pdb_to_gro():
            return False
            
        # Create membrane coordinates
        self.create_membrane_coordinates()
        
        # Create simple system
        self.create_simple_system()
        
        # Create test system
        self.create_test_system_without_membrane()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Martini 3 system prepared in {self.output_dir}")
        logger.info(f"{'='*60}")
        logger.info("\nNext steps:")
        logger.info("1. Test peptide topology:")
        logger.info(f"   cd {self.output_dir}")
        logger.info("   ./test_martini3.sh")
        logger.info("")
        logger.info("2. Add membrane using one of these methods:")
        logger.info("   - CHARMM-GUI Martini Maker (recommended)")
        logger.info("   - Download pre-equilibrated Martini 3 membrane")
        logger.info("   - See MEMBRANE_INSTRUCTIONS.txt")
        logger.info(f"{'='*60}\n")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/martini3_membrane')
    args = parser.parse_args()
    
    builder = Martini3MembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
