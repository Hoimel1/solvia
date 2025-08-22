"""
MD simulation workflow using GROMACS
Based on core scientific concept (Section 4.1) and best practices
"""

# MDP templates directory
MDP_DIR = Path("config/mdp_templates")

rule prepare_mdp_files:
    output:
        em_mdp = MDP_DIR / "em.mdp",
        nvt_mdp = MDP_DIR / "nvt.mdp",
        npt_mdp = MDP_DIR / "npt.mdp",
        prod_mdp = MDP_DIR / "prod.mdp"
    params:
        temp = config["simulation"]["temperature"],
        timestep = config["simulation"]["timestep"],
        prod_time = config["simulation"]["production_time"],
        output_freq = config["simulation"]["output_frequency"]
    run:
        os.makedirs(MDP_DIR, exist_ok=True)
        
        # Energy minimization
        with open(output.em_mdp, 'w') as f:
            f.write("""
; Energy minimization
integrator              = steep
emtol                   = 500.0
emstep                  = 0.01
nsteps                  = 50000

; Neighbor searching
nstlist                 = 20
cutoff-scheme           = Verlet
ns_type                 = grid
pbc                     = xyz
verlet-buffer-tolerance = 0.005

; Electrostatics and VdW
coulombtype             = reaction-field
rcoulomb                = 1.1
epsilon_r               = 15
epsilon_rf              = 0
vdw_type                = cutoff
vdw-modifier            = potential-shift-verlet
rvdw                    = 1.1
""")
        
        # NVT equilibration
        with open(output.nvt_mdp, 'w') as f:
            f.write(f"""
; NVT equilibration
define                  = -DPOSRES -DPOSRES_FULL
integrator              = md
dt                      = {params.timestep}
nsteps                  = 2500000  ; 50 ns
tinit                   = 0
init-step               = 0

; Output
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 50000
nstxout-compressed      = 50000
compressed-x-grps       = System
nstenergy               = 50000

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein Membrane Solvent
tau_t                   = 1.0 1.0 1.0
ref_t                   = {params.temp} {params.temp} {params.temp}

; Pressure coupling
pcoupl                  = no

; Constraints and parameters as before
""")
        
        # NPT equilibration
        with open(output.npt_mdp, 'w') as f:
            f.write(f"""
; NPT equilibration
define                  = -DPOSRES -DPOSRES_WEAK
integrator              = md
dt                      = {params.timestep}
nsteps                  = 5000000  ; 100 ns

; Output control
nstxout-compressed      = 50000
nstenergy               = 50000

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein Membrane Solvent
tau_t                   = 1.0 1.0 1.0
ref_t                   = {params.temp} {params.temp} {params.temp}

; Pressure coupling
pcoupl                  = berendsen
pcoupltype              = semiisotropic
tau_p                   = 12.0
ref_p                   = 1.0 1.0
compressibility         = 3e-4 3e-4
""")
        
        # Production
        with open(output.prod_mdp, 'w') as f:
            f.write(f"""
; Production run
integrator              = md
dt                      = {params.timestep}
nsteps                  = {int(params.prod_time * 1000 / params.timestep)}
tinit                   = 0

; Output
nstxout-compressed      = {int(params.output_freq / params.timestep)}
nstenergy               = {int(params.output_freq / params.timestep)}

; Temperature coupling
tcoupl                  = v-rescale
tc-grps                 = Protein Membrane Solvent
tau_t                   = 1.0 1.0 1.0
ref_t                   = {params.temp} {params.temp} {params.temp}

; Pressure coupling
pcoupl                  = parrinello-rahman
pcoupltype              = semiisotropic
tau_p                   = 12.0
ref_p                   = 1.0 1.0
compressibility         = 3e-4 3e-4
""")

rule run_md_simulation:
    input:
        system_gro = DATA_DIR / "intermediate" / "systems" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "system.gro",
        system_top = DATA_DIR / "intermediate" / "systems" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "system.top",
        em_mdp = MDP_DIR / "em.mdp",
        nvt_mdp = MDP_DIR / "nvt.mdp",
        npt_mdp = MDP_DIR / "npt.mdp",
        prod_mdp = MDP_DIR / "prod.mdp"
    output:
        traj = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "prod.xtc",
        edr = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "prod.edr",
        tpr = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "prod.tpr"
    params:
        prefix = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}",
        posres_script = "src/simulation/generate_posres.py"
    container:
        config["containers"]["gromacs"]
    threads: 8
    resources:
        gpu = 1,
        mem_gb = 16,
        time_hours = 24
    log:
        "logs/md_simulation/{amp_id}_{occupancy}_rep{replicate}.log"
    shell:
        """
        cd {params.prefix}
        
        # Generate position restraints
        python {params.posres_script} {input.system_gro}
        
        # Energy minimization
        gmx grompp -f {input.em_mdp} -c {input.system_gro} -p {input.system_top} \
                   -o em.tpr -maxwarn 2 >> {log} 2>&1
        gmx mdrun -deffnm em -v -nb gpu >> {log} 2>&1
        
        # NVT equilibration
        gmx grompp -f {input.nvt_mdp} -c em.gro -r em.gro -p {input.system_top} \
                   -o nvt.tpr -maxwarn 2 >> {log} 2>&1
        gmx mdrun -deffnm nvt -v -nb gpu -bonded gpu >> {log} 2>&1
        
        # NPT equilibration
        gmx grompp -f {input.npt_mdp} -c nvt.gro -r nvt.gro -t nvt.cpt -p {input.system_top} \
                   -o npt.tpr -maxwarn 2 >> {log} 2>&1
        gmx mdrun -deffnm npt -v -nb gpu -bonded gpu >> {log} 2>&1
        
        # Production
        gmx grompp -f {input.prod_mdp} -c npt.gro -t npt.cpt -p {input.system_top} \
                   -o {output.tpr} -maxwarn 2 >> {log} 2>&1
        gmx mdrun -deffnm prod -v -nb gpu -bonded gpu >> {log} 2>&1
        
        # Move outputs
        mv prod.xtc {output.traj}
        mv prod.edr {output.edr}
        """

