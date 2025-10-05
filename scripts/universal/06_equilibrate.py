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
import math
from copy import deepcopy
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
            "pmf_standard_config.yaml",
            "config.yaml"
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

    # Add missing default values (non-mutating)
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

def deep_merge(defaults: dict, overrides: dict | None) -> dict:
    """Return a deep-merged copy of defaults overlaid by overrides."""
    result = deepcopy(defaults)
    def rec(dst, src):
        for k, v in (src or {}).items():
            if isinstance(v, dict) and isinstance(dst.get(k), dict):
                rec(dst[k], v)
            else:
                dst[k] = v
    rec(result, overrides or {})
    return result


def add_config_defaults(config: dict | None) -> dict:
    """Return merged config (defaults + overrides) without mutating input."""
    defaults = get_default_config()
    return deep_merge(defaults, config or {})


def validate_and_normalize_config(cfg: dict) -> dict:
    """Validate basic simulation parameters and normalize types.

    Raises ValueError on invalid inputs. Returns cfg with a validation flag.
    """
    sim = cfg.get('simulation', {}) or {}
    # required positive floats
    for key in ('temperature', 'pressure', 'timestep'):
        try:
            v = float(sim.get(key, 0))
        except Exception:
            raise ValueError(f"simulation.{key} must be a float")
        if not (v > 0):
            raise ValueError(f"simulation.{key} must be > 0")
    # times (ps)
    for phase in ('nvt', 'npt'):
        try:
            t = float(sim.get(phase, {}).get('time', 0.0))
        except Exception:
            raise ValueError(f"simulation.{phase}.time must be a float")
        if not (t > 0):
            raise ValueError(f"simulation.{phase}.time must be > 0 ps")
    # tau_p
    try:
        tau_p = float(sim.get('npt', {}).get('tau_p', 0.0))
    except Exception:
        raise ValueError("simulation.npt.tau_p must be a float")
    if tau_p <= 0:
        raise ValueError("simulation.npt.tau_p must be > 0")
    # compressibility
    try:
        float(sim.get('npt', {}).get('compressibility', 4.5e-5))
    except Exception:
        raise ValueError("simulation.npt.compressibility must be float")

    perf = cfg.setdefault('performance', {})
    try:
        cpu_threads = int(perf.get('cpu_threads', 4))
    except Exception as exc:
        raise ValueError("performance.cpu_threads must be an integer") from exc
    perf['cpu_threads'] = max(1, cpu_threads)

    cfg['__validated__'] = True
    return cfg


def _nsteps(time_ps: float, dt_ps: float) -> int:
    """Robust integer step count avoiding FP off-by-one via floor with epsilon."""
    return int(math.floor((float(time_ps) / float(dt_ps)) + 1e-9))


def _docker_available(project_root: Path) -> bool:
    """Return True if docker compose is available and compose file exists."""
    try:
        r = subprocess.run(["docker", "compose", "version"], capture_output=True, text=True)
        return r.returncode == 0 and (project_root / "docker-compose.yml").exists()
    except Exception:
        return False


def _has_lincs_issue(diagnostics: dict | None) -> bool:
    """Detect classic LINCS/constraint failures from mdrun diagnostics."""
    if not diagnostics:
        return False
    text = diagnostics.get('combined') or ''
    lowered = text.lower()
    lincs_markers = (
        'lincs warning',
        'constraint error',
        'constraint failure',
        'relative constraint deviation',
        'bonds that rotated more than',
    )
    return any(marker in lowered for marker in lincs_markers)


class EquilibrationError(RuntimeError):
    pass

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def create_mdp_files(run_dir, config):
    """Create MDP files for equilibration"""
    mdp_dir = os.path.join(run_dir, "equilibration", "mdp")
    os.makedirs(mdp_dir, exist_ok=True)

    # Normalized locals for templating (avoid FP artifacts and injection)
    sim = config.get('simulation', {})
    dt = float(sim.get('timestep', 0.02))
    temp = float(sim.get('temperature', 310.0))
    press = float(sim.get('pressure', 1.0))
    nvt_steps = _nsteps(float(sim.get('nvt', {}).get('time', 2000.0)), dt)
    npt_steps = _nsteps(float(sim.get('npt', {}).get('time', 20000.0)), dt)
    tau_p_glob = float(sim.get('npt', {}).get('tau_p', 12.0))
    comp_glob = float(sim.get('npt', {}).get('compressibility', 4.5e-5))

    def posres_define(phase: str) -> str:
        pol = (
            config.get('simulation', {})
                  .get('posres', {})
                  .get(phase, 'on')
        )
        pol = str(pol).strip().lower()
        if pol in ('off', 'none', 'false', '0', 'no'):
            return '; Position restraints: off\n'
        if pol in ('weak', 'soft'):
            return 'define                  = -DPOSRES_WEAK\n'
        # default strong
        return 'define                  = -DPOSRES\n'

    def thermostat_block(phase: str) -> str:
        sim_local = config.get('simulation', {})
        tcfg = (sim_local.get('thermostat') or {})
        alg = str(tcfg.get('algorithm', 'v-rescale')).strip().lower()
        allowed_algos = {'v-rescale', 'nose-hoover'}
        if alg not in allowed_algos:
            logging.warning(f"Unsupported thermostat algorithm '{alg}', falling back to v-rescale")
            alg = 'v-rescale'
        groups = tcfg.get('tc_grps', 'System')
        if isinstance(groups, (list, tuple)):
            groups = ' '.join(map(str, groups))
        group_list = groups.split()
        ngrp = max(1, len(group_list))
        tau_cfg = tcfg.get('tau_t', None)
        if tau_cfg is None:
            base_tau = float(sim_local.get(phase, {}).get('tau_t', 1.0))
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
; epsilon_r=15, rcoulomb=1.1, rvdw=1.1 are Martini 3 recommendations.
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

; Electrostatics and VdW (Martini 3 standard)
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
; epsilon_r=15, rcoulomb=1.1, rvdw=1.1 are Martini 3 recommendations.
integrator              = md
dt                      = {dt}
nsteps                  = {nvt_steps}
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

; Electrostatics and VdW (Martini 3 standard)
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
gen-temp                = {temp}
gen-seed                = -1

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""

    # NPT equilibration MDP
    npt_mdp = f"""; NPT equilibration for Martini 3 membrane-peptide system
; epsilon_r=15, rcoulomb=1.1, rvdw=1.1 are Martini 3 recommendations.
integrator              = md
dt                      = {dt}
nsteps                  = {npt_steps}
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

; Electrostatics and VdW (Martini 3 standard)
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
vdw-type                = Cut-off
vdw-modifier            = Potential-shift-verlet
rvdw                    = 1.1
cutoff-scheme           = Verlet

; Temperature coupling
{thermostat_block('npt')}

; Pressure coupling (semiisotropic)
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau-p                   = {tau_p_glob}
ref-p                   = {press} {press}
compressibility         = {comp_glob} 0

; No velocity generation
gen-vel                 = no

; Constraints
constraints             = none
constraint-algorithm    = lincs
"""
    # Also prepare multistage NPT (Berendsen -> Parrinello-Rahman)
    sim = config.get('simulation', {})
    dt = float(sim.get('timestep', 0.02))
    total_steps = _nsteps(float(sim.get('npt', {}).get('time', 20000.0)), dt)
    ms = sim.get('npt_multistage', {}) or {}
    ms_enabled = bool(ms.get('enabled', True))
    frac = float(ms.get('initial_fraction', 0.3))
    steps_ber = max(0, int(total_steps * frac)) if ms_enabled else 0
    steps_pr = max(1, total_steps - steps_ber)
    tau_p_ber = float(ms.get('berendsen_tau_p', 5.0))
    tau_p_pr = float(sim.get('npt', {}).get('tau_p', 12.0))
    comp = sim.get('npt', {}).get('compressibility', 4.5e-5)
    press = sim.get('pressure', 1.0)
    logging.info(f"NPT multistage: total={total_steps} steps; berendsen={steps_ber}; pr={steps_pr} (frac={frac:.2f})")

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
    """Run GROMACS grompp via Docker (preferred) or host fallback."""
    # Get relative paths from project root
    project_root = Path(__file__).parent.parent.parent
    rel_mdp = os.path.relpath(mdp_file, project_root)
    rel_coord = os.path.relpath(coord_file, project_root)
    rel_top = os.path.relpath(top_file, project_root)
    rel_tpr = os.path.relpath(output_tpr, project_root)
    # Container working dir: place emergency files (e.g., step*.pdb) in target folder
    out_dir = os.path.dirname(output_tpr)
    rel_out_dir = os.path.relpath(out_dir, project_root)

    use_docker = _docker_available(project_root) and (os.environ.get("USE_DOCKER", "1") != "0")

    if use_docker:
        # gromacs container has ENTRYPOINT ["gmx"]. We call subcommand only.
        cmd = [
            "docker", "compose", "run", "--rm", "--workdir", f"/work/{rel_out_dir}", "gromacs",
            "grompp",
            "-f", f"/work/{rel_mdp}",
            "-c", f"/work/{rel_coord}",
            "-p", f"/work/{rel_top}",
            "-o", f"/work/{rel_tpr}",
            "-maxwarn", str(maxwarn)
        ]
        if ref_file:
            rel_ref = os.path.relpath(ref_file, project_root)
            cmd.extend(["-r", f"/work/{rel_ref}"])
        if cpt_file and os.path.exists(cpt_file):
            rel_cpt = os.path.relpath(cpt_file, project_root)
            cmd.extend(["-t", f"/work/{rel_cpt}"])
    else:
        # Host fallback: requires gmx in PATH
        cmd = [
            "gmx", "grompp",
            "-f", os.path.abspath(mdp_file),
            "-c", os.path.abspath(coord_file),
            "-p", os.path.abspath(top_file),
            "-o", os.path.abspath(output_tpr),
            "-maxwarn", str(maxwarn)
        ]
        if ref_file:
            cmd.extend(["-r", os.path.abspath(ref_file)])
        if cpt_file and os.path.exists(cpt_file):
            cmd.extend(["-t", os.path.abspath(cpt_file)])

    # Run; set cwd to out_dir so any emergency files (step*.pdb) land there
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=(str(project_root) if use_docker else out_dir or str(project_root)))
    if result.returncode != 0:
        logging.error("GROMPP failed:")
        if result.stdout:
            logging.error(f"STDOUT: {result.stdout}")
        if result.stderr:
            logging.error(f"STDERR: {result.stderr}")
        return False
    return True

def run_mdrun(tpr_file, output_prefix, nt=None, append_restart=False):
    """Run GROMACS mdrun via Docker (preferred) or host fallback.

    Returns (success: bool, diagnostics: dict).
    """
    project_root = Path(__file__).parent.parent.parent
    rel_prefix = os.path.relpath(output_prefix, project_root)
    out_dir = os.path.dirname(output_prefix)
    rel_out_dir = os.path.relpath(out_dir, project_root)

    use_docker = _docker_available(project_root) and (os.environ.get("USE_DOCKER", "1") != "0")
    base_cmd = [
        "docker", "compose", "run", "--rm", "--workdir", f"/work/{rel_out_dir}", "gromacs"
    ] if use_docker else ["gmx"]
    deffnm_arg = f"/work/{rel_prefix}" if use_docker else os.path.abspath(output_prefix)
    workdir = str(project_root) if use_docker else (out_dir or str(project_root))

    attempted_single_thread = False
    effective_nt = nt
    append_logged = False

    while True:
        cmd = base_cmd + ["mdrun", "-v", "-deffnm", deffnm_arg]
        if effective_nt:
            cmd.extend(["-nt", str(effective_nt)])
        if append_restart:
            cmd.append("-append")
            if not append_logged:
                logging.info("Using -append flag for checkpoint restart")
                append_logged = True

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
        stdout_text = result.stdout or ""
        stderr_text = result.stderr or ""
        combined = (stdout_text + "\n" + stderr_text).strip()
        diagnostics = {
            'stdout': stdout_text,
            'stderr': stderr_text,
            'combined': combined,
            'cmd': cmd,
            'nt': effective_nt,
        }

        if result.returncode == 0:
            return True, diagnostics

        lowered = combined.lower()
        domain_error = "domain decomposition" in lowered or "minimum cell size" in lowered
        if domain_error and not attempted_single_thread and effective_nt and int(effective_nt) > 1:
            logging.warning(
                "MDRUN domain decomposition failed with nt=%s; retrying with single thread.",
                effective_nt,
            )
            attempted_single_thread = True
            effective_nt = 1
            continue

        logging.error("MDRUN failed:")
        if stdout_text:
            logging.error(f"STDOUT: {stdout_text}")
        if stderr_text:
            logging.error(f"STDERR: {stderr_text}")
        return False, diagnostics

def equilibrate_system(run_dir, occupancy="low", tag: str | None = None, config=None):
    """Run complete equilibration workflow with robust error handling."""
    logger = logging.getLogger(__name__)

    if config is None:
        config = load_config()
    # ensure defaults and validate
    config = validate_and_normalize_config(add_config_defaults(config))
    metadata = load_run_metadata(run_dir)

    # Input files (support custom tag for replicates, e.g. n1_rep1)
    use_tag = tag if tag else occupancy
    system_dir = os.path.join(run_dir, "system")
    system_gro = os.path.join(system_dir, f"system_{use_tag}.gro")
    system_top = os.path.join(system_dir, f"system_{use_tag}.top")

    # Track progress for summary writing
    summary = {
        'peptide_id': metadata.get('peptide_id'),
        'tag': use_tag,
        'em_completed': False,
        'nvt_completed': False,
        'npt_completed': False,
        'final_structure': None,
        'final_checkpoint': None,
    }

    try:
        # If not found, try to find newest system_*.gro file
        if not os.path.exists(system_gro):
            import glob
            system_files = glob.glob(os.path.join(system_dir, "system_*.gro"))
            if system_files:
                system_files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
                system_gro = system_files[0]
                basename = os.path.basename(system_gro)
                use_tag = basename.replace("system_", "").replace(".gro", "")
                system_top = os.path.join(system_dir, f"system_{use_tag}.top")
                summary['tag'] = use_tag
                logger.info(f"Using system file: {basename}")
            else:
                raise EquilibrationError(f"No system files found in {system_dir}. Run peptide insertion first.")

        if not os.path.exists(system_gro):
            raise EquilibrationError(f"System file not found: {system_gro}")

        # Track current timestep for adaptive retries
        current_dt = float(config.get('simulation', {}).get('timestep', 0.02))
        initial_dt = current_dt
        config.setdefault('simulation', {})['timestep'] = current_dt
        max_stage_attempts = 3

        # Create MDP files
        mdp_dir = create_mdp_files(run_dir, config)

        def reduce_timestep(stage_label: str) -> bool:
            """Halve the timestep (down to 0.005 ps) and regenerate MDPs."""
            nonlocal current_dt, mdp_dir
            new_dt = max(current_dt / 2.0, 0.005)
            if new_dt >= current_dt - 1e-12:
                logger.error(
                    "%s failed and timestep cannot be reduced further (%.4f ps).",
                    stage_label,
                    current_dt,
                )
                return False
            current_dt = new_dt
            config['simulation']['timestep'] = current_dt
            logger.warning(
                "%s suffered LINCS issues; reducing timestep to %.4f ps and regenerating MDPs.",
                stage_label,
                current_dt,
            )
            mdp_dir = create_mdp_files(run_dir, config)
            return True

        # Energy minimization
        logger.info("Starting Energy Minimization")
        em_dir = os.path.join(run_dir, "equilibration", "em")
        os.makedirs(em_dir, exist_ok=True)
        em_mdp = os.path.join(mdp_dir, "em.mdp")
        em_tpr = os.path.join(em_dir, "em.tpr")

        logger.info("Running GROMPP for energy minimization")
        if not run_grompp(em_mdp, system_gro, system_top, em_tpr, ref_file=system_gro):
            raise EquilibrationError("EM grompp failed")

        logger.info("Running energy minimization")
        success, _ = run_mdrun(em_tpr, os.path.join(em_dir, "em"), config['performance']['cpu_threads'])
        if not success:
            raise EquilibrationError("Energy minimization failed")

        # Check EM results
        em_gro = os.path.join(em_dir, "em.gro")
        if os.path.exists(em_gro):
            logger.info("✓ Energy minimization completed successfully")
            summary['em_completed'] = True
        else:
            raise EquilibrationError("Energy minimization output not found")

        # NVT equilibration (allow timestep fallback on LINCS issues)
        logger.info("Starting NVT Equilibration")
        nvt_dir = os.path.join(run_dir, "equilibration", "nvt")
        os.makedirs(nvt_dir, exist_ok=True)
        nvt_time_ps = float(config.get('simulation', {}).get('nvt', {}).get('time', 2000.0))
        nvt_attempt = 0

        while True:
            nvt_attempt += 1
            nvt_mdp = os.path.join(mdp_dir, "nvt.mdp")
            nvt_tpr = os.path.join(nvt_dir, "nvt.tpr")

            logger.info("Running GROMPP for NVT equilibration (attempt %d)", nvt_attempt)
            if not run_grompp(nvt_mdp, em_gro, system_top, nvt_tpr, ref_file=em_gro):
                raise EquilibrationError("NVT grompp failed")

            logger.info(
                "Running NVT equilibration (%.1f ps @ %.4f ps timestep)",
                nvt_time_ps,
                current_dt,
            )
            nvt_success, nvt_diag = run_mdrun(
                nvt_tpr,
                os.path.join(nvt_dir, "nvt"),
                config['performance']['cpu_threads'],
            )

            if nvt_success:
                break

            if (
                _has_lincs_issue(nvt_diag)
                and nvt_attempt < max_stage_attempts
                and reduce_timestep("NVT equilibration")
            ):
                continue

            raise EquilibrationError("NVT equilibration failed")

        nvt_gro = os.path.join(nvt_dir, "nvt.gro")
        nvt_cpt = os.path.join(nvt_dir, "nvt.cpt")
        if os.path.exists(nvt_gro):
            logger.info("✓ NVT equilibration completed successfully")
            summary['nvt_completed'] = True
        else:
            raise EquilibrationError("NVT equilibration output not found")

        # NPT equilibration (multistage Berendsen -> Parrinello-Rahman)
        logger.info("Starting NPT Equilibration (multistage)")
        # Stage 1: Berendsen
        npt1_dir = os.path.join(run_dir, "equilibration", "npt_ber")
        os.makedirs(npt1_dir, exist_ok=True)
        npt1_attempt = 0

        while True:
            npt1_attempt += 1
            npt1_mdp = os.path.join(mdp_dir, "npt_ber.mdp")
            npt1_tpr = os.path.join(npt1_dir, "npt_ber.tpr")
            logger.info("Running GROMPP for NPT (Berendsen) (attempt %d)", npt1_attempt)
            if not run_grompp(npt1_mdp, nvt_gro, system_top, npt1_tpr, ref_file=nvt_gro, cpt_file=nvt_cpt):
                raise EquilibrationError("NPT (Berendsen) grompp failed")
            logger.info("Running NPT (Berendsen) stage (dt=%.4f ps)", current_dt)
            npt1_success, npt1_diag = run_mdrun(
                npt1_tpr,
                os.path.join(npt1_dir, "npt_ber"),
                config['performance']['cpu_threads'],
            )

            if npt1_success:
                break

            if (
                _has_lincs_issue(npt1_diag)
                and npt1_attempt < max_stage_attempts
                and reduce_timestep("NPT (Berendsen)")
            ):
                continue

            raise EquilibrationError("NPT (Berendsen) failed")

        npt1_gro = os.path.join(npt1_dir, "npt_ber.gro")
        npt1_cpt = os.path.join(npt1_dir, "npt_ber.cpt")
        if os.path.exists(npt1_gro):
            logger.info("✓ NPT (Berendsen) completed successfully")
        else:
            raise EquilibrationError("NPT (Berendsen) output not found")

        # Stage 2: Parrinello-Rahman
        npt2_dir = os.path.join(run_dir, "equilibration", "npt_pr")
        os.makedirs(npt2_dir, exist_ok=True)
        npt2_attempt = 0

        while True:
            npt2_attempt += 1
            npt2_mdp = os.path.join(mdp_dir, "npt_pr.mdp")
            npt2_tpr = os.path.join(npt2_dir, "npt_pr.tpr")
            logger.info("Running GROMPP for NPT (Parrinello-Rahman) (attempt %d)", npt2_attempt)
            if not run_grompp(npt2_mdp, npt1_gro, system_top, npt2_tpr, ref_file=npt1_gro, cpt_file=npt1_cpt):
                raise EquilibrationError("NPT (Parrinello-Rahman) grompp failed")
            logger.info("Running NPT (Parrinello-Rahman) stage (dt=%.4f ps)", current_dt)
            npt2_success, npt2_diag = run_mdrun(
                npt2_tpr,
                os.path.join(npt2_dir, "npt_pr"),
                config['performance']['cpu_threads'],
            )

            if npt2_success:
                break

            if (
                _has_lincs_issue(npt2_diag)
                and npt2_attempt < max_stage_attempts
                and reduce_timestep("NPT (Parrinello-Rahman)")
            ):
                continue

            raise EquilibrationError("NPT (Parrinello-Rahman) failed")

        npt2_gro = os.path.join(npt2_dir, "npt_pr.gro")
        npt2_cpt = os.path.join(npt2_dir, "npt_pr.cpt")
        if os.path.exists(npt2_gro):
            logger.info("✓ NPT (Parrinello-Rahman) completed successfully")
            summary['npt_completed'] = True
            summary['final_structure'] = 'equilibration/npt_pr/npt_pr.gro'
            summary['final_checkpoint'] = 'equilibration/npt_pr/npt_pr.cpt'
        else:
            raise EquilibrationError("NPT (Parrinello-Rahman) output not found")

        if abs(current_dt - initial_dt) > 1e-6:
            logger.info(
                "Adaptive timestep: reduced from %.4f ps to %.4f ps during equilibration.",
                initial_dt,
                current_dt,
            )

        logger.info(f"\u2713 Equilibration complete for {metadata['peptide_id']} ({use_tag})")
        return npt2_gro

    finally:
        # Always persist a summary of what completed
        try:
            equil_dir = os.path.join(run_dir, "equilibration")
            os.makedirs(equil_dir, exist_ok=True)
            with open(os.path.join(equil_dir, "summary.yaml"), 'w') as f:
                yaml.dump(summary, f, default_flow_style=False)
        except Exception as e:
            logger.warning(f"Could not write equilibration summary: {e}")

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

    # Run equilibration with error handling
    try:
        final_gro = equilibrate_system(args.run_dir, args.occupancy, tag=args.tag, config=config)
    except EquilibrationError as e:
        logger.error(str(e))
        sys.exit(1)
    
    next_tag = args.tag if args.tag else args.occupancy
    logger.info(f"Next step: Run production simulation")
    logger.info(f"Command: python 07_run_production.py {args.run_dir} --tag {next_tag}")

if __name__ == "__main__":
    main()
