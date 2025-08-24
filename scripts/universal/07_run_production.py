#!/usr/bin/env python3
"""
SOLVIA Production Simulation
Runs production MD simulation for hemolytic toxicity analysis
"""

import os
import sys
import yaml
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

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

def create_production_mdp(run_dir, config, simulation_time=None):
    """Create MDP file for production run"""
    if simulation_time is None:
        simulation_time = config['simulation']['production']['time']
    
    # Calculate number of steps
    timestep = config['simulation']['timestep']
    nsteps = int(simulation_time * 1000 / timestep)  # Convert ns to ps
    
    # Output frequency
    output_freq = int(config['simulation']['production']['output_frequency'] / timestep)
    
    production_mdp = f"""; Production MD for Martini 3 membrane-peptide system
; {simulation_time} ns simulation at {config['simulation']['temperature']}K, {config['simulation']['pressure']} bar

integrator              = md
dt                      = {timestep}
nsteps                  = {nsteps}
nstcomm                 = 100

; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = {output_freq}
nstcalcenergy          = 100
nstenergy              = {output_freq}
nstxout-compressed     = {output_freq}
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
tau-t                   = 1.0
ref-t                   = {config['simulation']['temperature']}

; Pressure coupling
pcoupl                  = parrinello-rahman
pcoupltype              = semiisotropic
tau-p                   = 12.0
ref-p                   = {config['simulation']['pressure']} {config['simulation']['pressure']}
compressibility         = 4.5e-5 4.5e-5

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
"""
    
    # Write MDP file
    mdp_file = os.path.join(run_dir, "production", "production.mdp")
    with open(mdp_file, 'w') as f:
        f.write(production_mdp)
    
    return mdp_file, nsteps

def run_production(run_dir, occupancy="low", time_ns=None, gpu=True, tag: str | None = None):
    """Run production simulation"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Check if equilibration was completed
    equil_summary = os.path.join(run_dir, "equilibration", "summary.yaml")
    if not os.path.exists(equil_summary):
        print("Error: Equilibration not completed. Run equilibration first.")
        sys.exit(1)
    
    with open(equil_summary, 'r') as f:
        equil_data = yaml.safe_load(f)
    
    if not equil_data.get('npt_completed'):
        print("Error: NPT equilibration not completed")
        sys.exit(1)
    
    # Input files
    npt_gro = os.path.join(run_dir, equil_data['final_structure'])
    npt_cpt = os.path.join(run_dir, equil_data['final_checkpoint'])
    use_tag = tag if tag else occupancy
    system_top = os.path.join(run_dir, "system", f"system_{use_tag}.top")
    
    # Output directory
    prod_dir = os.path.join(run_dir, "production")
    
    # Create production MDP
    print(f"Setting up production run...")
    mdp_file, nsteps = create_production_mdp(run_dir, config, time_ns)
    
    # Run grompp
    tpr_file = os.path.join(prod_dir, "production.tpr")
    cmd = [
        "gmx", "grompp",
        "-f", mdp_file,
        "-c", npt_gro,
        "-t", npt_cpt,
        "-p", system_top,
        "-o", tpr_file,
        "-maxwarn", "2"
    ]
    
    print("Running grompp...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"GROMPP Error:\n{result.stderr}")
        sys.exit(1)
    
    # Prepare mdrun command
    mdrun_cmd = [
        "gmx", "mdrun",
        "-v",
        "-deffnm", os.path.join(prod_dir, "production")
    ]
    
    # Use GPU if available and requested
    if gpu and config['performance']['gpu']:
        # Note: Martini force field does not use PME for electrostatics
        mdrun_cmd.extend(["-nb", "gpu"])
    else:
        mdrun_cmd.extend(["-nt", str(config['performance']['cpu_threads'])])
    
    # Calculate estimated runtime
    print(f"\n{'='*60}")
    print(f"Starting production simulation:")
    print(f"  Peptide: {metadata['peptide_id']}")
    print(f"  System tag: {use_tag}")
    print(f"  Duration: {time_ns or config['simulation']['production']['time']} ns ({nsteps:,} steps)")
    print(f"  Output every: {config['simulation']['production']['output_frequency']} ps")
    print(f"  GPU acceleration: {'Enabled' if gpu and config['performance']['gpu'] else 'Disabled'}")
    print(f"{'='*60}\n")
    
    # Create simulation status file
    status = {
        'peptide_id': metadata['peptide_id'],
        'occupancy': occupancy,
        'simulation_time_ns': time_ns or config['simulation']['production']['time'],
        'start_time': datetime.now().isoformat(),
        'status': 'running',
        'gpu_enabled': gpu and config['performance']['gpu']
    }
    
    status_file = os.path.join(prod_dir, "simulation_status.yaml")
    with open(status_file, 'w') as f:
        yaml.dump(status, f, default_flow_style=False)
    
    # Run simulation in foreground and stream output to terminal and log
    print("Starting mdrun...")
    print(f"Command: {' '.join(mdrun_cmd)}")
    print("\nOutput is streamed below and also written to production.log")
    log_file = os.path.join(prod_dir, "production.log")
    with open(log_file, 'w') as log:
        process = subprocess.Popen(mdrun_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        try:
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log.write(line)
                log.flush()
        finally:
            process.stdout.close()
            ret = process.wait()
        if ret != 0:
            print("Production simulation failed")
            status['status'] = 'failed'
        else:
            print("\n✓ Production simulation completed successfully")
            status['status'] = 'completed'
            status['end_time'] = datetime.now().isoformat()
    
    # Update status
    with open(status_file, 'w') as f:
        yaml.dump(status, f, default_flow_style=False)
    
    return os.path.join(prod_dir, "production.xtc")

def main():
    parser = argparse.ArgumentParser(
        description="Run production simulation for SOLVIA"
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
        help="Custom system tag to select system_*.top (e.g., n1_rep1)"
    )
    parser.add_argument(
        "--time",
        type=float,
        help="Simulation time in nanoseconds (default: from config)"
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable GPU acceleration"
    )
    
    args = parser.parse_args()
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Run production
    trajectory = run_production(
        args.run_dir, 
        args.occupancy, 
        args.time,
        gpu=not args.no_gpu,
        tag=args.tag
    )
    
    print(f"\nWhen simulation completes:")
    print(f"  1. Check trajectory: gmx check -f {trajectory}")
    print(f"  2. Run analysis: python 08_analyze_trajectory.py {args.run_dir} --occupancy {args.occupancy}")

if __name__ == "__main__":
    main()
