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
import logging
from pathlib import Path

def setup_logging(level=logging.INFO):
    """Setup logging configuration"""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

def load_config(config_file=None):
    """Load SOLVIA configuration with defaults"""
    config_dir = Path(__file__).parent.parent.parent / "config"
    
    if config_file:
        # Use specified config file
        if os.path.isabs(config_file):
            config_path = Path(config_file)
        else:
            config_path = config_dir / config_file
    else:
        # Try available config files in order of preference
        config_files = [
            "config.yaml",
            "pmf_standard_config.yaml"
        ]
        
        config_path = None
        for config_file in config_files:
            potential_path = config_dir / config_file
            if potential_path.exists():
                config_path = potential_path
                break
    
    if not config_path or not config_path.exists():
        logging.warning("No config file found, using defaults")
        return get_default_config()
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    # Add missing default values
    config = add_config_defaults(config)
    logging.info(f"Loaded configuration from: {config_path}")
    return config

def get_default_config():
    """Get default configuration for equilibration"""
    return {
        'simulation': {
            'temperature': 310,  # K
            'pressure': 1.0,     # bar
            'timestep': 0.02,    # ps (20 fs for Martini 3)
            'em': {
                'emtol': 10.0,
                'nsteps': 50000
            },
            'nvt': {
                'time': 2000.0,  # ps (>=2 ns typical)
                'tau_t': 1.0     # ps
            },
            'npt': {
                'time': 20000.0,  # ps (>=20 ns typical)
                'tau_t': 1.0,    # ps
                'tau_p': 12.0,   # ps
                'compressibility': 4.5e-5  # bar^-1 (RBC study reference)
            },
            # Position restraint policy per phase (on/weak/off)
            'posres': {
                'em': 'on',
                'nvt': 'on',
                'npt': 'off'
            },
            # Optional thermostat grouping (defaults to single group 'System')
            'thermostat': {
                'algorithm': 'v-rescale',
                'tc_grps': 'System'
                # When multiple groups are used, set tau_t/ref_t arrays to match group count.
            },
            # Optional multistage NPT (Berendsen -> Parrinello-Rahman)
            'npt_multistage': {
                'enabled': True,
                'initial_fraction': 0.3,
                'berendsen_tau_p': 5.0
            }
        },
        'performance': {
            'cpu_threads': 4
        }
    }

def add_config_defaults(config):
    """Add default values for missing configuration keys"""
    defaults = get_default_config()
    
    # Deep merge defaults
    def merge_dict(base, updates):
        for key, value in updates.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                merge_dict(base[key], value)
            elif key not in base:
                base[key] = value
    
    merge_dict(config, defaults)
    return config

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def create_mdp_files(run_dir, config):
    """Create MDP files for equilibration"""
    mdp_dir = os.path.join(run_dir, "equilibration", "mdp")
    os.makedirs(mdp_dir, exist_ok=True)
    
    def posres_define(phase: str) -> str:
        pol = (
            config.get('simulation', {})
                  .get('posres', {})
                  .get(phase, 'on')
        )
        pol = str(pol).strip().lower()
        if pol in ('off', 'none', 'false', '0', 'no'):
            return '; Position restraints: off\n'
        # weak/strong both rely on POSRES macro here; strength is controlled in topology/itp
        return 'define                  = -DPOSRES\n'

    def thermostat_block(phase: str) -> str:
        sim = config.get('simulation', {})
        tcfg = (sim.get('thermostat') or {})
        alg = tcfg.get('algorithm', 'v-rescale')
        groups = tcfg.get('tc_grps', 'System')
        temp = float(sim.get('temperature', 310))
        if isinstance(groups, (list, tuple)):
            groups = ' '.join(map(str, groups))
        group_list = groups.split()
        ngrp = max(1, len(group_list))
        tau_cfg = tcfg.get('tau_t', None)
        if tau_cfg is None:
            base_tau = float(sim.get(phase, {}).get('tau_t', 1.0))
            tau_str = ' '.join([f"{base_tau}"] * ngrp)
        else:
            if isinstance(tau_cfg, (list, tuple)):
                tau_str = ' '.join(map(str, tau_cfg))
            else:
                tau_str = str(tau_cfg)
        ref_cfg = tcfg.get('ref_t', None)
        if ref_cfg is None:
            ref_str = ' '.join([f"{temp}"] * ngrp)
        else:
            if isinstance(ref_cfg, (list, tuple)):
                ref_str = ' '.join(map(str, ref_cfg))
            else:
                ref_str = str(ref_cfg)
        return f"""tcoupl                  = {alg}
tc-grps                = {groups}
tau-t                  = {tau_str}
ref-t                  = {ref_str}
"""
    
    # Energy minimization MDP
    em_mdp = f"""; Energy minimization for Martini 3 membrane-peptide system
integrator              = steep
emtol                   = {config['simulation']['em']['emtol']}
emstep                  = 0.01
nsteps                  = {config['simulation']['em']['nsteps']}

; Position restraints on peptides to prevent them from moving below membrane
{posres_define('em')}

; Neighbor searching
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW (study reference)
coulombtype             = reaction-field
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
constraint-algorithm    = lincs
"""
    
    # NVT equilibration MDP
    nvt_mdp = f"""; NVT equilibration for Martini 3 membrane-peptide system
integrator              = md
dt                      = {config['simulation']['timestep']}
nsteps                  = {int(config['simulation']['nvt']['time'] / config['simulation']['timestep'])}
nstcomm                 = 100

; Position restraints on peptides
{posres_define('nvt')}

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
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW (study reference)
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
{thermostat_block('nvt')}

; Pressure coupling off
pcoupl                  = no

; Generate velocities
gen-vel                 = yes
gen-temp                = {config['simulation']['temperature']}
gen-seed                = -1

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""
    
    # NPT equilibration MDP  
    npt_mdp = f"""; NPT equilibration for Martini 3 membrane-peptide system
integrator              = md
dt                      = {config['simulation']['timestep']}
nsteps                  = {int(config['simulation']['npt']['time'] / config['simulation']['timestep'])}
nstcomm                 = 100

; Position restraints on peptides (can be removed or weakened)
{posres_define('npt')}

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
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW (study reference)
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
{thermostat_block('npt')}

; Pressure coupling (study reference)
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau-p                   = {config['simulation']['npt']['tau_p']}
ref-p                   = {config['simulation']['pressure']} {config['simulation']['pressure']}
compressibility         = {config['simulation']['npt']['compressibility']} 0

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""
    # Also prepare multistage NPT (Berendsen -> Parrinello-Rahman)
    sim = config.get('simulation', {})
    dt = float(sim.get('timestep', 0.02))
    total_steps = int(sim.get('npt', {}).get('time', 20000.0) / dt)
    ms = sim.get('npt_multistage', {}) or {}
    ms_enabled = bool(ms.get('enabled', True))
    frac = float(ms.get('initial_fraction', 0.3))
    steps_ber = max(0, int(total_steps * frac)) if ms_enabled else 0
    steps_pr = max(1, total_steps - steps_ber)
    tau_p_ber = float(ms.get('berendsen_tau_p', 5.0))
    tau_p_pr = float(sim.get('npt', {}).get('tau_p', 12.0))
    comp = sim.get('npt', {}).get('compressibility', 4.5e-5)
    press = sim.get('pressure', 1.0)

    npt_ber_mdp = f"""; NPT initial (Berendsen)
integrator              = md
dt                      = {dt}
nsteps                  = {steps_ber}
nstcomm                 = 100

; Position restraints (initial stage)
{posres_define('npt_initial') if 'npt_initial' in sim.get('posres', {}) else posres_define('npt')}

; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstcalcenergy           = 100
nstenergy               = 500
nstxout-compressed      = 500
compressed-x-precision  = 100

; Neighbor searching
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
{thermostat_block('npt')}

; Pressure coupling (initial)
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau-p                   = {tau_p_ber}
ref-p                   = {press} {press}
compressibility         = {comp} 0

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""

    npt_pr_mdp = f"""; NPT final (Parrinello-Rahman)
integrator              = md
dt                      = {dt}
nsteps                  = {steps_pr}
nstcomm                 = 100

; Position restraints (final stage)
{posres_define('npt_final') if 'npt_final' in sim.get('posres', {}) else posres_define('npt')}

; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 500
nstcalcenergy           = 100
nstenergy               = 500
nstxout-compressed      = 500
compressed-x-precision  = 100

; Neighbor searching
nstlist                 = 20
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
{thermostat_block('npt')}

; Pressure coupling (final)
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau-p                   = {tau_p_pr}
ref-p                   = {press} {press}
compressibility         = {comp} 0

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""
    
    # Write MDP files
    with open(os.path.join(mdp_dir, "em.mdp"), 'w') as f:
        f.write(em_mdp)
    with open(os.path.join(mdp_dir, "nvt.mdp"), 'w') as f:
        f.write(nvt_mdp)
    with open(os.path.join(mdp_dir, "npt.mdp"), 'w') as f:
        f.write(npt_mdp)
    # Write multistage variants as well
    with open(os.path.join(mdp_dir, "npt_ber.mdp"), 'w') as f:
        f.write(npt_ber_mdp)
    with open(os.path.join(mdp_dir, "npt_pr.mdp"), 'w') as f:
        f.write(npt_pr_mdp)
    
    return mdp_dir

def run_grompp(mdp_file, coord_file, top_file, output_tpr, 
               ref_file=None, cpt_file=None, maxwarn=2):
    """Run GROMACS grompp via Docker"""
    from pathlib import Path
    
    # Get relative paths from project root
    project_root = Path(__file__).parent.parent.parent
    rel_mdp = os.path.relpath(mdp_file, project_root)
    rel_coord = os.path.relpath(coord_file, project_root)
    rel_top = os.path.relpath(top_file, project_root)
    rel_tpr = os.path.relpath(output_tpr, project_root)
    
    # Build Docker command (gromacs container has ENTRYPOINT ["gmx"])
    cmd = [
        "docker", "compose", "run", "--rm", "gromacs",
        "grompp",  # Just the subcommand - Docker already provides "gmx"
        "-f", f"/work/{rel_mdp}",
        "-c", f"/work/{rel_coord}",
        "-p", f"/work/{rel_top}",
        "-o", f"/work/{rel_tpr}",
        "-maxwarn", str(maxwarn)
    ]
    # Add include paths so local and centralized ITPs resolve cleanly
    # Always include force_fields and the run's coarse_grain directory
    # Note: GROMACS 2023.5 does not support -I on grompp. We rely on absolute
    # include paths in generated .top files instead.
    
    if ref_file:
        rel_ref = os.path.relpath(ref_file, project_root)
        cmd.extend(["-r", f"/work/{rel_ref}"])
    if cpt_file and os.path.exists(cpt_file):
        rel_cpt = os.path.relpath(cpt_file, project_root)
        cmd.extend(["-t", f"/work/{rel_cpt}"])
    
    # Run from project root (where docker-compose.yml is)
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(project_root))
    if result.returncode != 0:
        logging.error("GROMPP failed:")
        if result.stdout:
            logging.error(f"STDOUT: {result.stdout}")
        if result.stderr:
            logging.error(f"STDERR: {result.stderr}")
        return False
    return True

def run_mdrun(tpr_file, output_prefix, nt=None, append_restart=False):
    """Run GROMACS mdrun via Docker"""
    from pathlib import Path
    
    # Get relative paths from project root
    project_root = Path(__file__).parent.parent.parent
    rel_tpr = os.path.relpath(tpr_file, project_root)
    rel_prefix = os.path.relpath(output_prefix, project_root)
    
    # Build Docker command (gromacs container has ENTRYPOINT ["gmx"])
    cmd = [
        "docker", "compose", "run", "--rm", "gromacs",
        "mdrun",  # Just the subcommand - Docker already provides "gmx"
        "-v",
        "-deffnm", f"/work/{rel_prefix}"
    ]
    
    if nt:
        cmd.extend(["-nt", str(nt)])
    
    if append_restart:
        cmd.append("-append")
        logging.info("Using -append flag for checkpoint restart")
    
    # Run from project root (where docker-compose.yml is)
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(project_root))
    if result.returncode != 0:
        logging.error(f"MDRUN failed:")
        if result.stdout:
            logging.error(f"STDOUT: {result.stdout}")
        if result.stderr:
            logging.error(f"STDERR: {result.stderr}")
        return False
    return True

def equilibrate_system(run_dir, occupancy="low", tag: str | None = None, config=None):
    """Run complete equilibration workflow"""
    logger = logging.getLogger(__name__)
    
    if config is None:
        config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Input files (support custom tag for replicates, e.g. n1_rep1)
    use_tag = tag if tag else occupancy
    system_dir = os.path.join(run_dir, "system")
    system_gro = os.path.join(system_dir, f"system_{use_tag}.gro")
    system_top = os.path.join(system_dir, f"system_{use_tag}.top")
    
    # If not found, try to find any system file
    if not os.path.exists(system_gro):
        # Look for any system_*.gro file
        import glob
        system_files = glob.glob(os.path.join(system_dir, "system_*.gro"))
        if system_files:
            # Use the first (or most recent) system file found
            system_gro = system_files[0]
            basename = os.path.basename(system_gro)
            use_tag = basename.replace("system_", "").replace(".gro", "")
            system_top = os.path.join(system_dir, f"system_{use_tag}.top")
            logger.info(f"Using system file: {basename}")
        else:
            logger.error(f"No system files found in {system_dir}")
            logger.error("Run peptide insertion first: python 05_insert_peptides.py")
            sys.exit(1)
    
    if not os.path.exists(system_gro):
        logger.error(f"System file not found: {system_gro}")
        sys.exit(1)
    
    # Create MDP files
    mdp_dir = create_mdp_files(run_dir, config)
    
    # Energy minimization
    logger.info("Starting Energy Minimization")
    em_dir = os.path.join(run_dir, "equilibration", "em")
    os.makedirs(em_dir, exist_ok=True)  # Ensure directory exists
    em_mdp = os.path.join(mdp_dir, "em.mdp")
    em_tpr = os.path.join(em_dir, "em.tpr")
    
    logger.info("Running GROMPP for energy minimization")
    if not run_grompp(em_mdp, system_gro, system_top, em_tpr, ref_file=system_gro):
        logger.error("EM grompp failed")
        sys.exit(1)
    
    logger.info("Running energy minimization")
    if not run_mdrun(em_tpr, os.path.join(em_dir, "em"), 
                     config['performance']['cpu_threads']):
        logger.error("Energy minimization failed")
        sys.exit(1)
    
    # Check EM results
    em_gro = os.path.join(em_dir, "em.gro")
    if os.path.exists(em_gro):
        logger.info("✓ Energy minimization completed successfully")
    else:
        logger.error("✗ Energy minimization output not found")
        sys.exit(1)
    
    # NVT equilibration
    logger.info("Starting NVT Equilibration")
    nvt_dir = os.path.join(run_dir, "equilibration", "nvt")
    os.makedirs(nvt_dir, exist_ok=True)  # Ensure directory exists
    nvt_mdp = os.path.join(mdp_dir, "nvt.mdp")
    nvt_tpr = os.path.join(nvt_dir, "nvt.tpr")
    
    logger.info("Running GROMPP for NVT equilibration")
    if not run_grompp(nvt_mdp, em_gro, system_top, nvt_tpr, ref_file=em_gro):
        logger.error("NVT grompp failed")
        sys.exit(1)
    
    logger.info(f"Running NVT equilibration ({config['simulation']['nvt']['time']} ps)")
    if not run_mdrun(nvt_tpr, os.path.join(nvt_dir, "nvt"),
                     config['performance']['cpu_threads']):
        logger.error("NVT equilibration failed")
        sys.exit(1)
    
    nvt_gro = os.path.join(nvt_dir, "nvt.gro")
    nvt_cpt = os.path.join(nvt_dir, "nvt.cpt")
    if os.path.exists(nvt_gro):
        logger.info("✓ NVT equilibration completed successfully")
    else:
        logger.error("✗ NVT equilibration output not found")
        sys.exit(1)
    
    # NPT equilibration (multistage Berendsen -> Parrinello-Rahman)
    logger.info("Starting NPT Equilibration (multistage)")
    # Stage 1: Berendsen
    npt1_dir = os.path.join(run_dir, "equilibration", "npt_ber")
    os.makedirs(npt1_dir, exist_ok=True)
    npt1_mdp = os.path.join(mdp_dir, "npt_ber.mdp")
    npt1_tpr = os.path.join(npt1_dir, "npt_ber.tpr")
    logger.info("Running GROMPP for NPT (Berendsen)")
    if not run_grompp(npt1_mdp, nvt_gro, system_top, npt1_tpr, ref_file=nvt_gro, cpt_file=nvt_cpt):
        logger.error("NPT (Berendsen) grompp failed")
        sys.exit(1)
    logger.info("Running NPT (Berendsen) stage")
    if not run_mdrun(npt1_tpr, os.path.join(npt1_dir, "npt_ber"), config['performance']['cpu_threads']):
        logger.error("NPT (Berendsen) failed")
        sys.exit(1)
    npt1_gro = os.path.join(npt1_dir, "npt_ber.gro")
    npt1_cpt = os.path.join(npt1_dir, "npt_ber.cpt")
    if os.path.exists(npt1_gro):
        logger.info("✓ NPT (Berendsen) completed successfully")
    else:
        logger.error("✗ NPT (Berendsen) output not found")
        sys.exit(1)

    # Stage 2: Parrinello-Rahman
    npt2_dir = os.path.join(run_dir, "equilibration", "npt_pr")
    os.makedirs(npt2_dir, exist_ok=True)
    npt2_mdp = os.path.join(mdp_dir, "npt_pr.mdp")
    npt2_tpr = os.path.join(npt2_dir, "npt_pr.tpr")
    logger.info("Running GROMPP for NPT (Parrinello-Rahman)")
    if not run_grompp(npt2_mdp, npt1_gro, system_top, npt2_tpr, ref_file=npt1_gro, cpt_file=npt1_cpt):
        logger.error("NPT (Parrinello-Rahman) grompp failed")
        sys.exit(1)
    logger.info("Running NPT (Parrinello-Rahman) stage")
    if not run_mdrun(npt2_tpr, os.path.join(npt2_dir, "npt_pr"), config['performance']['cpu_threads']):
        logger.error("NPT (Parrinello-Rahman) failed")
        sys.exit(1)
    npt2_gro = os.path.join(npt2_dir, "npt_pr.gro")
    if os.path.exists(npt2_gro):
        logger.info("✓ NPT (Parrinello-Rahman) completed successfully")
    else:
        logger.error("✗ NPT (Parrinello-Rahman) output not found")
        sys.exit(1)
    
    logger.info(f"\u2713 Equilibration complete for {metadata['peptide_id']} ({use_tag})")
    
    # Save equilibration summary
    summary = {
        'peptide_id': metadata['peptide_id'],
        'tag': use_tag,
        'em_completed': os.path.exists(em_gro),
        'nvt_completed': os.path.exists(nvt_gro),
        'npt_completed': os.path.exists(npt2_gro),
        'final_structure': 'equilibration/npt_pr/npt_pr.gro',
        'final_checkpoint': 'equilibration/npt_pr/npt_pr.cpt'
    }
    
    with open(os.path.join(run_dir, "equilibration", "summary.yaml"), 'w') as f:
        yaml.dump(summary, f, default_flow_style=False)
    
    return npt2_gro

def main():
    parser = argparse.ArgumentParser(
        description="Run equilibration for SOLVIA system"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory"
    )
    parser.add_argument(
        "--config",
        help="Configuration file to use (default: search for config.yaml or pmf_standard_config.yaml)"
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
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(log_level)
    
    # Check if run directory exists
    if not os.path.exists(args.run_dir):
        logger.error(f"Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Load configuration
    config = load_config(args.config)
    
    # Run equilibration
    final_gro = equilibrate_system(args.run_dir, args.occupancy, tag=args.tag, config=config)
    
    next_tag = args.tag if args.tag else args.occupancy
    logger.info(f"Next step: Run production simulation")
    logger.info(f"Command: python 07_run_production.py {args.run_dir} --tag {next_tag}")

if __name__ == "__main__":
    main()
