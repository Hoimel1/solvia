#!/usr/bin/env python3
"""
Build a physiologically accurate RBC membrane with asymmetric composition.

Outer leaflet: ~45% POPC, 10% PSM (Sphingomyelin), 45% Cholesterol
Inner leaflet: ~45% POPE, 15% POPS, 40% Cholesterol
Total cholesterol: 40-50 mol%
Thickness: 4.5-5.0 nm

References: Plisson et al., 2020; Bennett et al., 2020
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RBCMembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Membrane composition (approximate numbers for 12x12 nm membrane)
        self.box_x = 12.0  # nm
        self.box_y = 12.0  # nm
        self.area = self.box_x * self.box_y  # nm²
        self.lipids_per_nm2 = 1.5  # approximate for Martini
        self.total_lipids_per_leaflet = int(self.area * self.lipids_per_nm2)
        
        # Calculate lipid numbers based on percentages
        # Outer leaflet
        self.n_popc_outer = int(0.45 * self.total_lipids_per_leaflet)
        self.n_psm_outer = int(0.10 * self.total_lipids_per_leaflet)
        self.n_chol_outer = int(0.45 * self.total_lipids_per_leaflet)
        
        # Inner leaflet  
        self.n_pope_inner = int(0.45 * self.total_lipids_per_leaflet)
        self.n_pops_inner = int(0.15 * self.total_lipids_per_leaflet)
        self.n_chol_inner = int(0.40 * self.total_lipids_per_leaflet)
        
        logger.info(f"RBC Membrane composition:")
        logger.info(f"Outer leaflet: {self.n_popc_outer} POPC, {self.n_psm_outer} PSM, {self.n_chol_outer} CHOL")
        logger.info(f"Inner leaflet: {self.n_pope_inner} POPE, {self.n_pops_inner} POPS, {self.n_chol_inner} CHOL")
        
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
        
        # Copy peptide topology with correct name
        top_src = Path(f"data/processed/topologies/{self.peptide_id}/{self.peptide_id}.itp")
        if not top_src.exists():
            # Try old name
            top_src = Path(f"data/processed/topologies/{self.peptide_id}/molecule_0.itp")
        
        if top_src.exists():
            shutil.copy(top_src, self.output_dir / f"{self.peptide_id}.itp")
    
    def create_topology(self):
        """Create the system topology for RBC membrane."""
        logger.info("Creating RBC membrane topology")
        
        # Check if we have sterols file for cholesterol
        has_sterols = (self.output_dir / "martini_v3.0_sterols_v1.0.itp").exists()
        
        topology = f"""
; Include force field parameters
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
"""
        
        if has_sterols:
            topology += '#include "martini_v3.0_sterols_v1.0.itp"\n'
        
        topology += f"""
; Include peptide topology
#include "{self.peptide_id}.itp"

[ system ]
; Name
16x {self.peptide_id} on asymmetric RBC membrane

[ molecules ]
; Compound        #mols
{self.peptide_id}    16

; Outer leaflet (top)
POPC             {self.n_popc_outer}
; PSM would go here if available in Martini 3
POPC             {self.n_psm_outer}  ; Using POPC as PSM substitute
"""
        
        if has_sterols:
            topology += f"CHOL             {self.n_chol_outer}\n"
        else:
            topology += f"; CHOL not available - need martini_v3.0_sterols_v1.0.itp\n"
        
        topology += f"""
; Inner leaflet (bottom)
POPE             {self.n_pope_inner}
POPS             {self.n_pops_inner}
"""
        
        if has_sterols:
            topology += f"CHOL             {self.n_chol_inner}\n"
        else:
            topology += f"; CHOL not available - need martini_v3.0_sterols_v1.0.itp\n"
        
        topology += """
; Solvent and ions (to be added)
W                20000    ; Water beads
NA+              100      ; Sodium ions
CL-              100      ; Chloride ions
"""
        
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("RBC membrane topology created")
        
        # Check lipid availability
        self.check_lipid_availability()
        
        return True
    
    def check_lipid_availability(self):
        """Check which lipids are available in the force field."""
        phospholipids_file = self.output_dir / "martini_v3.0.0_phospholipids_v1.itp"
        
        if phospholipids_file.exists():
            with open(phospholipids_file, 'r') as f:
                content = f.read()
                
            available_lipids = []
            missing_lipids = []
            
            lipids_to_check = ['POPC', 'POPE', 'POPS', 'PSM', 'CHOL']
            
            for lipid in lipids_to_check:
                if f"[ moleculetype ]\n; molname      nrexcl\n  {lipid}" in content or f"{lipid} " in content:
                    available_lipids.append(lipid)
                else:
                    missing_lipids.append(lipid)
            
            logger.info(f"Available lipids: {', '.join(available_lipids)}")
            if missing_lipids:
                logger.warning(f"Missing lipids: {', '.join(missing_lipids)}")
                logger.warning("PSM (Sphingomyelin) might need to be substituted with similar lipid")
    
    def create_building_instructions(self):
        """Create detailed instructions for building the RBC membrane."""
        
        instructions = f"""
RBC MEMBRANE BUILDING INSTRUCTIONS
==================================

Target composition:
- Outer leaflet: ~45% POPC, 10% PSM, 45% Cholesterol  
- Inner leaflet: ~45% POPE, 15% POPS, 40% Cholesterol

Calculated numbers for {self.box_x}x{self.box_y} nm membrane:
- Outer: {self.n_popc_outer} POPC, {self.n_psm_outer} PSM, {self.n_chol_outer} CHOL
- Inner: {self.n_pope_inner} POPE, {self.n_pops_inner} POPS, {self.n_chol_inner} CHOL

METHOD 1: CHARMM-GUI Martini Maker (Recommended)
------------------------------------------------
1. Go to http://www.charmm-gui.org/
2. Select "Martini Maker" → "Membrane Builder"
3. Choose "Asymmetric Bilayer"
4. Set membrane composition:
   - Upper leaflet: POPC (45%), PSM (10%), Cholesterol (45%)
   - Lower leaflet: POPE (45%), POPS (15%), Cholesterol (40%)
5. Upload peptides_16x.pdb as "Protein"
6. Set box size to 12x12x15 nm
7. Add 150 mM NaCl
8. Download and convert to GROMACS format

METHOD 2: insane.py (if available)
----------------------------------
# Build asymmetric membrane
python insane.py \\
    -f peptides_16x.pdb \\
    -o system.gro \\
    -p system.top \\
    -x 12 -y 12 -z 15 \\
    -l POPC:{self.n_popc_outer} -l POPE:{self.n_pope_inner} \\
    -l POPS:{self.n_pops_inner} -l CHOL:{self.n_chol_outer + self.n_chol_inner} \\
    -u POPC:45 -u PSM:10 -u CHOL:45 \\
    -l POPE:45 -l POPS:15 -l CHOL:40 \\
    -sol W -salt 0.15

METHOD 3: Manual Construction with PACKMOL
------------------------------------------
1. Use PACKMOL to arrange lipids in two leaflets
2. Place peptides above membrane
3. Solvate with Martini water beads
4. Add ions to 150 mM NaCl

NOTES:
- PSM (Sphingomyelin) may not be available in Martini 3
  Consider using DPPC or additional POPC as substitute
- Ensure cholesterol flip-flop is restricted during equilibration
- Use semi-isotropic pressure coupling for membrane systems
- Apply position restraints on peptides during initial equilibration

EQUILIBRATION PROTOCOL:
1. Energy minimization (10,000 steps)
2. NVT equilibration (1 ns) - restrain peptides and lipid headgroups
3. NPT equilibration (5 ns) - semi-isotropic pressure coupling
4. NPT equilibration (5 ns) - release peptide restraints
5. Production run (100-500 ns)
"""
        
        with open(self.output_dir / "RBC_MEMBRANE_INSTRUCTIONS.txt", 'w') as f:
            f.write(instructions)
        
        logger.info("Created detailed RBC membrane building instructions")
        return True
    
    def create_mdp_files(self):
        """Create MDP files optimized for asymmetric RBC membrane."""
        
        # Energy minimization - very gentle for complex membrane
        em_mdp = """
; Energy minimization for RBC membrane with peptides
define                  = -DPOSRES  ; Position restrain peptides initially
integrator              = steep
nsteps                  = 20000
emtol                   = 100.0     ; Lower tolerance for complex system
emstep                  = 0.001     ; Smaller step size

; Output
nstlog                  = 1000
nstenergy               = 1000

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
        
        # NVT equilibration - careful temperature coupling
        nvt_mdp = """
; NVT equilibration for asymmetric RBC membrane
define                  = -DPOSRES -DPOSRES_FC=1000  ; Restrain peptides
integrator              = md
dt                      = 0.01      ; Smaller timestep initially
nsteps                  = 100000    ; 1 ns

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 5000

; Temperature coupling - separate groups for asymmetric membrane
tcoupl                  = v-rescale
tc-grps                 = Protein POPC_POPE_PSM POPS CHOL W_ION
tau_t                   = 1.0 1.0 1.0 1.0 1.0
ref_t                   = 310 310 310 310 310

; No pressure coupling in NVT
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
gen_seed                = -1
"""
        
        # NPT equilibration - semi-isotropic for membrane
        npt_mdp = """
; NPT equilibration for asymmetric RBC membrane
define                  = -DPOSRES_FC=500  ; Weaker restraints
integrator              = md
dt                      = 0.02
nsteps                  = 250000    ; 5 ns

; Output
nstlog                  = 1000
nstenergy               = 1000
nstxout-compressed      = 5000

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein POPC_POPE_PSM POPS CHOL W_ION
tau_t                   = 1.0 1.0 1.0 1.0 1.0
ref_t                   = 310 310 310 310 310

; Pressure coupling - semi-isotropic for membrane
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 5.0       ; Slower coupling for equilibration
ref_p                   = 1.0 1.0   ; x-y and z separately
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
        
        # Production run - no restraints
        md_mdp = """
; Production MD for RBC membrane with peptides
integrator              = md
dt                      = 0.025     ; 25 fs for Martini 3
nsteps                  = 20000000  ; 500 ns

; Output - reduced frequency
nstlog                  = 10000
nstenergy               = 1000
nstxout-compressed      = 10000
compressed-x-grps       = System

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein Membrane W_ION
tau_t                   = 1.0 1.0 1.0
ref_t                   = 310 310 310

; Pressure coupling - Parrinello-Rahman for production
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
        
        logger.info("Created MDP files for RBC membrane simulation")
        return True
    
    def create_analysis_scripts(self):
        """Create analysis scripts for membrane-peptide interactions."""
        
        analysis_script = """#!/bin/bash
# Analysis for RBC membrane-peptide system

echo "RBC Membrane-Peptide Analysis"
echo "============================="

# 1. System setup info
echo "Analyzing system composition..."
echo -e "0\\n" | gmx check -f md.xtc 2>&1 | grep -E "Coords|Box"

# 2. Peptide-membrane distance
echo "Calculating peptide-membrane distances..."
echo -e "1\\n2\\n" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi

# 3. Membrane thickness
echo "Analyzing membrane thickness..."
# Would need gmx density or similar tool

# 4. Peptide orientation
echo "Analyzing peptide orientation..."
echo -e "1\\n" | gmx gangle -f md.xtc -s md.tpr -g1 vector -g2 z -oav peptide_angle.xvg

# 5. Lipid order parameters
echo "Note: For detailed lipid analysis, use MDAnalysis or similar tools"

# 6. Create visualization files
echo "Creating visualization files..."
echo -e "1\\n2\\n0\\n" | gmx trjconv -f md.xtc -s md.tpr -o viz_peptide_membrane.xtc -pbc mol -ur compact -dt 1000

# Final structure
echo -e "0\\n" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo ""
echo "Key files for visualization:"
echo "- final_structure.pdb: Final frame in PDB format"
echo "- viz_peptide_membrane.xtc: Trajectory for VMD (every 1 ns)"
echo ""
echo "VMD visualization commands:"
echo "vmd final_structure.pdb viz_peptide_membrane.xtc"
echo ""
echo "In VMD:"
echo "- Peptides: 'resname $PEPTIDE_ID'"
echo "- Membrane: 'resname POPC POPE POPS CHOL'"
echo "- Use 'Drawing Method' -> 'QuickSurf' for membrane"
"""
        
        with open(self.output_dir / "analyze_membrane.sh", 'w') as f:
            f.write(analysis_script.replace('$PEPTIDE_ID', self.peptide_id))
        os.chmod(self.output_dir / "analyze_membrane.sh", 0o755)
        
        # Create feature extraction script
        feature_script = f"""#!/usr/bin/env python3
'''Extract features from RBC membrane simulation for {self.peptide_id}'''

import numpy as np
import subprocess
import json

def extract_features():
    features = {{}}
    
    # Basic trajectory info
    cmd = ['gmx', 'check', '-f', 'md.xtc']
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Parse peptide-membrane distances
    if os.path.exists('peptide_membrane_dist.xvg'):
        distances = []
        with open('peptide_membrane_dist.xvg', 'r') as f:
            for line in f:
                if not line.startswith('#') and not line.startswith('@'):
                    parts = line.split()
                    if len(parts) >= 2:
                        distances.append(float(parts[1]))
        
        if distances:
            features['mean_distance'] = np.mean(distances)
            features['min_distance'] = np.min(distances)
            features['contact_frequency'] = sum(1 for d in distances if d < 0.5) / len(distances)
    
    # Save features
    with open('membrane_features.json', 'w') as f:
        json.dump(features, f, indent=2)
    
    print(f"Extracted features: {{features}}")
    return features

if __name__ == "__main__":
    extract_features()
"""
        
        with open(self.output_dir / "extract_features.py", 'w') as f:
            f.write(feature_script)
        os.chmod(self.output_dir / "extract_features.py", 0o755)
        
        logger.info("Created analysis scripts")
        return True
    
    def run(self):
        """Run the complete RBC membrane building pipeline."""
        logger.info(f"Building RBC membrane system for 16x {self.peptide_id}")
        
        # Create topology
        self.create_topology()
        
        # Create building instructions
        self.create_building_instructions()
        
        # Create MDP files
        self.create_mdp_files()
        
        # Create analysis scripts
        self.create_analysis_scripts()
        
        logger.info(f"RBC membrane setup complete in {self.output_dir}")
        logger.info("Next steps:")
        logger.info("1. Build membrane using CHARMM-GUI following RBC_MEMBRANE_INSTRUCTIONS.txt")
        logger.info("2. Run simulations with the provided MDP files")
        logger.info("3. Analyze with analyze_membrane.sh")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide-id', default='SOLVIA_1')
    parser.add_argument('--output-dir', default='simulations/rbc_membrane')
    args = parser.parse_args()
    
    builder = RBCMembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
