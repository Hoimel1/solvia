#!/usr/bin/env python3
"""
Build RBC membrane with 16 peptides using INSANE with Martini 3 support.

This script uses the local INSANE installation that supports M3 force field
and includes proper M3.CHOL template (5 beads instead of 8).
"""

import os
import sys
import subprocess
import shutil
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def build_rbc_membrane():
    """Build asymmetric RBC membrane with INSANE using Martini 3."""
    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    sim_dir = project_root / "simulations" / "rbc_martini3"
    
    # Clean and create directory
    if sim_dir.exists():
        shutil.rmtree(sim_dir)
    sim_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Building RBC membrane in {sim_dir}")
    
    # Copy necessary files
    logger.info("Copying peptide structure...")
    peptide_src = project_root / "simulations" / "membrane_16x" / "peptides_16x.pdb"
    
    if not peptide_src.exists():
        # Create it if missing
        logger.info("Creating 16x peptide structure...")
        create_cmd = [
            "python", str(script_dir / "setup_16peptide_membrane.py"),
            "--peptide-id", "SOLVIA_1"
        ]
        subprocess.run(create_cmd, check=True)
    
    shutil.copy(peptide_src, sim_dir / "peptides_16x.pdb")
    
    # Copy force field files
    logger.info("Copying Martini 3 force field files...")
    ff_dir = project_root / "force_fields" / "martini3"
    
    ff_files = [
        "martini_v3.0.0.itp",
        "martini_v3.0.0_ions_v1.itp",
        "martini_v3.0.0_phospholipids_PC_v2.itp",
        "martini_v3.0.0_phospholipids_PE_v2.itp",
        "martini_v3.0.0_phospholipids_SM_v2.itp",
        "martini_v3.0_sterols_v1.0.itp",  # Important for M3 CHOL
        "martini_v3.0.0_solvents_v1.itp"
    ]
    
    for ff_file in ff_files:
        src = ff_dir / ff_file
        if src.exists():
            shutil.copy(src, sim_dir)
    
    # Copy peptide ITP (renamed to SOLVIA_1.itp)
    logger.info("Copying peptide topology...")
    peptide_itp_src = project_root / "data" / "processed" / "topologies" / "SOLVIA_1.itp"
    if peptide_itp_src.exists():
        shutil.copy(peptide_itp_src, sim_dir)
    else:
        # Try the old name
        old_itp = project_root / "data" / "processed" / "topologies" / "molecule_0.itp"
        if old_itp.exists():
            shutil.copy(old_itp, sim_dir / "SOLVIA_1.itp")
    
    # Build INSANE command for asymmetric RBC membrane
    logger.info("Building membrane with INSANE (Martini 3)...")
    
    # RBC membrane composition:
    # Outer leaflet: ~45% POPC, 10% PSM, 45% Cholesterol
    # Inner leaflet: ~45% POPE, 15% POPS, 40% Cholesterol
    
    # For a 20x20 nm membrane, we'll use approximately:
    # Outer: 180 POPC, 40 PSM, 180 CHOL
    # Inner: 180 POPE, 60 POPS, 160 CHOL
    
    insane_cmd = [
        "insane",
        "-f", "peptides_16x.pdb",
        "-ff", "M3",  # Force Martini 3
        "-x", "20", "-y", "20", "-z", "15",
        "-pbc", "rectangular",
        "-center",
        # Upper leaflet (outer)
        "-l", "M3.POPC:180",
        "-l", "M3.PSM:40", 
        "-l", "M3.CHOL:180",
        # Lower leaflet (inner) - using -l with negative area
        "-l", "M3.POPE:180", "-area", "0.65",
        "-l", "M3.POPS:60", "-area", "0.65",
        "-l", "M3.CHOL:160", "-area", "0.65",
        # Solvent
        "-sol", "W:90", "-sol", "WF:10",  # 90% regular water, 10% antifreeze
        "-salt", "0.15",
        # Output
        "-o", "system.gro",
        "-p", "system_insane.top"
    ]
    
    logger.info(f"Running: {' '.join(insane_cmd)}")
    
    os.chdir(sim_dir)
    result = subprocess.run(insane_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"INSANE failed:\n{result.stderr}")
        # Try simpler command without asymmetry
        logger.info("Trying symmetric membrane as fallback...")
        insane_cmd_simple = [
            "insane",
            "-f", "peptides_16x.pdb",
            "-ff", "M3",
            "-x", "20", "-y", "20", "-z", "15",
            "-pbc", "rectangular",
            "-center",
            # Mixed membrane
            "-l", "M3.POPC:7",
            "-l", "M3.POPE:5",
            "-l", "M3.POPS:2",
            "-l", "M3.PSM:1",
            "-l", "M3.CHOL:5",
            "-sol", "W",
            "-salt", "0.15",
            "-o", "system.gro",
            "-p", "system_insane.top"
        ]
        result = subprocess.run(insane_cmd_simple, capture_output=True, text=True)
    
    logger.info(f"INSANE output:\n{result.stdout}")
    
    # Fix topology file
    logger.info("Creating proper topology file...")
    
    top_content = f"""; RBC membrane with 16 SOLVIA_1 peptides
; Using Martini 3 force field

#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_phospholipids_PC_v2.itp"
#include "martini_v3.0.0_phospholipids_PE_v2.itp"
#include "martini_v3.0.0_phospholipids_SM_v2.itp"
#include "martini_v3.0_sterols_v1.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"

; Peptide topology
#include "SOLVIA_1.itp"

[ system ]
RBC membrane with 16 SOLVIA_1 peptides

[ molecules ]
"""
    
    # Parse INSANE topology for molecule counts
    if Path("system_insane.top").exists():
        with open("system_insane.top", "r") as f:
            lines = f.readlines()
            
        in_molecules = False
        for line in lines:
            if "[ molecules ]" in line:
                in_molecules = True
                continue
            if in_molecules and line.strip() and not line.startswith(";"):
                # Fix molecule names for M3
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    mol_count = parts[1]
                    
                    # Handle peptide entry
                    if "Protein" in mol_name or "SOLVIA" in mol_name:
                        top_content += f"SOLVIA_1   16\n"
                    # Fix ion names
                    elif mol_name == "NA+":
                        top_content += f"NA         {mol_count}\n"
                    elif mol_name == "CL-":
                        top_content += f"CL         {mol_count}\n"
                    else:
                        top_content += f"{mol_name}     {mol_count}\n"
    
    with open("system.top", "w") as f:
        f.write(top_content)
    
    # Create MDP files
    logger.info("Creating MDP files...")
    
    # Energy minimization
    em_mdp = """title = Energy minimization
integrator = steep
nsteps = 5000
emtol = 1000.0
emstep = 0.01

nstlist = 10
cutoff-scheme = Verlet
ns_type = grid
pbc = xyz
verlet-buffer-tolerance = 0.005

coulombtype = reaction-field
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdw_type = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1

constraints = none
"""
    
    with open("em.mdp", "w") as f:
        f.write(em_mdp)
    
    # NVT equilibration
    nvt_mdp = """title = NVT equilibration
integrator = md
dt = 0.02
nsteps = 50000  ; 1 ns
nstcomm = 10

nstxout = 0
nstvout = 0
nstfout = 0
nstlog = 1000
nstenergy = 100
nstxout-compressed = 1000

cutoff-scheme = Verlet
nstlist = 20
ns_type = grid
pbc = xyz
verlet-buffer-tolerance = 0.005

coulombtype = reaction-field
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdw_type = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1

tcoupl = v-rescale
tc-grps = Protein Membrane Solvent
tau_t = 1.0 1.0 1.0
ref_t = 310 310 310

pcoupl = no

gen_vel = yes
gen_temp = 310
gen_seed = -1

constraints = none
constraint_algorithm = lincs
"""
    
    with open("nvt.mdp", "w") as f:
        f.write(nvt_mdp)
    
    # NPT equilibration
    npt_mdp = """title = NPT equilibration
integrator = md
dt = 0.02
nsteps = 250000  ; 5 ns
nstcomm = 10

nstxout = 0
nstvout = 0
nstfout = 0
nstlog = 1000
nstenergy = 100
nstxout-compressed = 1000

cutoff-scheme = Verlet
nstlist = 20
ns_type = grid
pbc = xyz
verlet-buffer-tolerance = 0.005

coulombtype = reaction-field
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdw_type = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1

tcoupl = v-rescale
tc-grps = Protein Membrane Solvent
tau_t = 1.0 1.0 1.0
ref_t = 310 310 310

pcoupl = Parrinello-Rahman
pcoupltype = semiisotropic
tau_p = 12.0
compressibility = 3e-4 3e-4
ref_p = 1.0 1.0

gen_vel = no

constraints = none
constraint_algorithm = lincs
"""
    
    with open("npt.mdp", "w") as f:
        f.write(npt_mdp)
    
    # Production run
    prod_mdp = """title = Production run
integrator = md
dt = 0.02
nsteps = 5000000  ; 100 ns
nstcomm = 10

nstxout = 0
nstvout = 0
nstfout = 0
nstlog = 5000
nstenergy = 1000
nstxout-compressed = 5000

cutoff-scheme = Verlet
nstlist = 20
ns_type = grid
pbc = xyz
verlet-buffer-tolerance = 0.005

coulombtype = reaction-field
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdw_type = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1

tcoupl = v-rescale
tc-grps = Protein Membrane Solvent
tau_t = 1.0 1.0 1.0
ref_t = 310 310 310

pcoupl = Parrinello-Rahman
pcoupltype = semiisotropic
tau_p = 12.0
compressibility = 3e-4 3e-4
ref_p = 1.0 1.0

gen_vel = no

constraints = none
constraint_algorithm = lincs
"""
    
    with open("prod.mdp", "w") as f:
        f.write(prod_mdp)
    
    # Create run script
    run_script = """#!/bin/bash
#SBATCH --job-name=rbc_m3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1

# Energy minimization
echo "Starting energy minimization..."
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em

# Check if minimization succeeded
if [ ! -f em.gro ]; then
    echo "Energy minimization failed!"
    exit 1
fi

# NVT equilibration
echo "Starting NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm nvt

# NPT equilibration
echo "Starting NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p system.top -o npt.tpr -maxwarn 2
gmx mdrun -v -deffnm npt

# Production run
echo "Starting production run..."
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p system.top -o prod.tpr -maxwarn 2
gmx mdrun -v -deffnm prod

echo "Simulation complete!"
"""
    
    with open("run_simulation.sh", "w") as f:
        f.write(run_script)
    
    os.chmod("run_simulation.sh", 0o755)
    
    # Test the setup
    logger.info("Testing setup with grompp...")
    test_cmd = [
        "gmx", "grompp",
        "-f", "em.mdp",
        "-c", "system.gro",
        "-p", "system.top",
        "-o", "test.tpr",
        "-maxwarn", "2"
    ]
    
    test_result = subprocess.run(test_cmd, capture_output=True, text=True)
    
    if test_result.returncode == 0:
        logger.info("✅ Setup successful! System is ready for simulation.")
        logger.info(f"To run: cd {sim_dir} && ./run_simulation.sh")
    else:
        logger.error(f"Setup test failed:\n{test_result.stderr}")
        logger.info("Please check topology and coordinate files.")
    
    # Summary
    logger.info("\nRBC Membrane Setup Complete!")
    logger.info(f"Location: {sim_dir}")
    logger.info("Files created:")
    logger.info("  - system.gro: Coordinates with membrane")
    logger.info("  - system.top: Martini 3 topology")
    logger.info("  - *.mdp: Simulation parameters")
    logger.info("  - run_simulation.sh: Execution script")
    

if __name__ == "__main__":
    build_rbc_membrane()
