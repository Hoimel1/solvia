#!/usr/bin/env python3
"""
SOLVIA System Equilibration
Performs energy minimization, NVT and NPT equilibration
"""

import os
import sys
import yaml
import argparse
import subprocess
from pathlib import Path

def load_config():
    """Load SOLVIA configuration"""
    config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def create_mdp_files(run_dir, config):
    """Create MDP files for equilibration"""
    mdp_dir = os.path.join(run_dir, "equilibration", "mdp")
    os.makedirs(mdp_dir, exist_ok=True)
    
    # Energy minimization MDP
    em_mdp = f"""; Energy minimization for Martini 3 membrane-peptide system
integrator              = steep
emtol                   = {config['simulation']['em']['emtol']}
emstep                  = 0.01
nsteps                  = {config['simulation']['em']['nsteps']}

; Position restraints on peptides to prevent them from moving below membrane
define                  = -DPOSRES

; Neighbor searching
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = Reaction-Field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; No temperature or pressure coupling
tcoupl                  = no
pcoupl                  = no

; No constraints during EM
constraints             = none
"""
    
    # NVT equilibration MDP
    nvt_mdp = f"""; NVT equilibration for Martini 3 membrane-peptide system
integrator              = md
dt                      = {config['simulation']['timestep']}
nsteps                  = {int(config['simulation']['nvt']['time'] / config['simulation']['timestep'])}
nstcomm                 = 100

; Position restraints on peptides
define                  = -DPOSRES

; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstcalcenergy          = 100
nstenergy              = 500
nstxout-compressed     = 500
compressed-x-precision = 100

; Neighbor searching
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = Reaction-Field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau-t                   = {config['simulation']['nvt']['tau_t']}
ref-t                   = {config['simulation']['temperature']}

; Pressure coupling off
pcoupl                  = no

; Generate velocities
gen-vel                 = yes
gen-temp                = {config['simulation']['temperature']}
gen-seed                = -1

; Constraints
constraints             = none
"""
    
    # NPT equilibration MDP
    npt_mdp = f"""; NPT equilibration for Martini 3 membrane-peptide system
integrator              = md
dt                      = {config['simulation']['timestep']}
nsteps                  = {int(config['simulation']['npt']['time'] / config['simulation']['timestep'])}
nstcomm                 = 100

; Position restraints on peptides (can be removed or weakened)
define                  = -DPOSRES

; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstcalcenergy          = 100
nstenergy              = 500
nstxout-compressed     = 500
compressed-x-precision = 100

; Neighbor searching
nstlist                 = 20
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = Reaction-Field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = System
tau-t                   = {config['simulation']['npt']['tau_t']}
ref-t                   = {config['simulation']['temperature']}

; Pressure coupling
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau-p                   = {config['simulation']['npt']['tau_p']}
ref-p                   = {config['simulation']['pressure']} {config['simulation']['pressure']}
compressibility         = {config['simulation']['npt']['compressibility']} {config['simulation']['npt']['compressibility']}

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
"""
    
    # Write MDP files
    with open(os.path.join(mdp_dir, "em.mdp"), 'w') as f:
        f.write(em_mdp)
    with open(os.path.join(mdp_dir, "nvt.mdp"), 'w') as f:
        f.write(nvt_mdp)
    with open(os.path.join(mdp_dir, "npt.mdp"), 'w') as f:
        f.write(npt_mdp)
    
    return mdp_dir

def run_grompp(mdp_file, coord_file, top_file, output_tpr, 
               ref_file=None, cpt_file=None, maxwarn=2):
    """Run GROMACS grompp"""
    cmd = [
        "gmx", "grompp",
        "-f", mdp_file,
        "-c", coord_file,
        "-p", top_file,
        "-o", output_tpr,
        "-maxwarn", str(maxwarn)
    ]
    
    if ref_file:
        cmd.extend(["-r", ref_file])
    if cpt_file and os.path.exists(cpt_file):
        cmd.extend(["-t", cpt_file])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"GROMPP Error:\n{result.stderr}")
        return False
    return True

def run_mdrun(tpr_file, output_prefix, nt=None):
    """Run GROMACS mdrun"""
    cmd = [
        "gmx", "mdrun",
        "-v",
        "-deffnm", output_prefix
    ]
    
    if nt:
        cmd.extend(["-nt", str(nt)])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"MDRUN Error:\n{result.stderr}")
        return False
    return True

def equilibrate_system(run_dir, occupancy="low", tag: str | None = None):
    """Run complete equilibration workflow"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Input files (support custom tag for replicates, e.g. n1_rep1)
    use_tag = tag if tag else occupancy
    system_gro = os.path.join(run_dir, "system", f"system_{use_tag}.gro")
    system_top = os.path.join(run_dir, "system", f"system_{use_tag}.top")
    
    if not os.path.exists(system_gro):
        print(f"Error: System file not found: {system_gro}")
        print("Run peptide insertion first: python 05_insert_peptides.py")
        sys.exit(1)
    
    # Create MDP files
    mdp_dir = create_mdp_files(run_dir, config)
    
    # Energy minimization
    print("\n=== Energy Minimization ===")
    em_dir = os.path.join(run_dir, "equilibration", "em")
    em_mdp = os.path.join(mdp_dir, "em.mdp")
    em_tpr = os.path.join(em_dir, "em.tpr")
    
    print("Running grompp...")
    if not run_grompp(em_mdp, system_gro, system_top, em_tpr, ref_file=system_gro):
        print("EM grompp failed")
        sys.exit(1)
    
    print("Running energy minimization...")
    if not run_mdrun(em_tpr, os.path.join(em_dir, "em"), 
                     config['performance']['cpu_threads']):
        print("Energy minimization failed")
        sys.exit(1)
    
    # Check EM results
    em_gro = os.path.join(em_dir, "em.gro")
    if os.path.exists(em_gro):
        print("✓ Energy minimization completed successfully")
    else:
        print("✗ Energy minimization output not found")
        sys.exit(1)
    
    # NVT equilibration
    print("\n=== NVT Equilibration ===")
    nvt_dir = os.path.join(run_dir, "equilibration", "nvt")
    nvt_mdp = os.path.join(mdp_dir, "nvt.mdp")
    nvt_tpr = os.path.join(nvt_dir, "nvt.tpr")
    
    print("Running grompp...")
    if not run_grompp(nvt_mdp, em_gro, system_top, nvt_tpr, ref_file=em_gro):
        print("NVT grompp failed")
        sys.exit(1)
    
    print(f"Running NVT equilibration ({config['simulation']['nvt']['time']} ps)...")
    if not run_mdrun(nvt_tpr, os.path.join(nvt_dir, "nvt"),
                     config['performance']['cpu_threads']):
        print("NVT equilibration failed")
        sys.exit(1)
    
    nvt_gro = os.path.join(nvt_dir, "nvt.gro")
    nvt_cpt = os.path.join(nvt_dir, "nvt.cpt")
    if os.path.exists(nvt_gro):
        print("✓ NVT equilibration completed successfully")
    else:
        print("✗ NVT equilibration output not found")
        sys.exit(1)
    
    # NPT equilibration
    print("\n=== NPT Equilibration ===")
    npt_dir = os.path.join(run_dir, "equilibration", "npt")
    npt_mdp = os.path.join(mdp_dir, "npt.mdp")
    npt_tpr = os.path.join(npt_dir, "npt.tpr")
    
    print("Running grompp...")
    if not run_grompp(npt_mdp, nvt_gro, system_top, npt_tpr, 
                      ref_file=nvt_gro, cpt_file=nvt_cpt):
        print("NPT grompp failed")
        sys.exit(1)
    
    print(f"Running NPT equilibration ({config['simulation']['npt']['time']} ps)...")
    if not run_mdrun(npt_tpr, os.path.join(npt_dir, "npt"),
                     config['performance']['cpu_threads']):
        print("NPT equilibration failed")
        sys.exit(1)
    
    npt_gro = os.path.join(npt_dir, "npt.gro")
    if os.path.exists(npt_gro):
        print("✓ NPT equilibration completed successfully")
    else:
        print("✗ NPT equilibration output not found")
        sys.exit(1)
    
    print(f"\n✓ Equilibration complete for {metadata['peptide_id']} ({use_tag})")
    
    # Save equilibration summary
    summary = {
        'peptide_id': metadata['peptide_id'],
        'tag': use_tag,
        'em_completed': os.path.exists(em_gro),
        'nvt_completed': os.path.exists(nvt_gro),
        'npt_completed': os.path.exists(npt_gro),
        'final_structure': 'equilibration/npt/npt.gro',
        'final_checkpoint': 'equilibration/npt/npt.cpt'
    }
    
    with open(os.path.join(run_dir, "equilibration", "summary.yaml"), 'w') as f:
        yaml.dump(summary, f, default_flow_style=False)
    
    return npt_gro

def main():
    parser = argparse.ArgumentParser(
        description="Run equilibration for SOLVIA system"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory"
    )
    parser.add_argument(
        "--occupancy",
        choices=["low", "medium", "high"],
        default="low",
        help="Peptide occupancy level (default: low)"
    )
    parser.add_argument(
        "--tag",
        help="Custom system tag, e.g. n1_rep1 (overrides occupancy-based filenames)"
    )
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Run equilibration
    final_gro = equilibrate_system(args.run_dir, args.occupancy, tag=args.tag)
    
    next_tag = args.tag if args.tag else args.occupancy
    print(f"\nNext step: Run production simulation")
    print(f"Command: python 07_run_production.py {args.run_dir} --tag {next_tag}")

if __name__ == "__main__":
    main()
