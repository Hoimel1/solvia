# SOLVIA Snakemake Pipeline
# Automated workflow for hemolytic toxicity prediction of antimicrobial peptides

import os
import yaml
from pathlib import Path
import glob

# Load configuration
configfile: "config/solvia_config.yaml"

# Define input FASTA files
FASTA_DIR = "data/raw/fasta_split"
ALL_PEPTIDES = sorted([f.replace('.fasta', '') for f in os.listdir(FASTA_DIR) if f.endswith('.fasta')])
# For testing, start with SOLVIA_1
PEPTIDES = ["SOLVIA_1"] if "SOLVIA_1" in ALL_PEPTIDES else ALL_PEPTIDES[:1]

# PMF validation study peptides
PMF_TOXIC = ["SOLVIA_1", "SOLVIA_8", "SOLVIA_14", "SOLVIA_215", "SOLVIA_164", "SOLVIA_126", 
             "SOLVIA_68", "SOLVIA_32", "SOLVIA_482", "SOLVIA_490", "SOLVIA_515", "SOLVIA_524",
             "SOLVIA_527", "SOLVIA_617", "SOLVIA_624", "SOLVIA_850", "SOLVIA_858", "SOLVIA_941",
             "SOLVIA_974", "SOLVIA_1023", "SOLVIA_1045"]
PMF_NONTOXIC = ["SOLVIA_12", "SOLVIA_398", "SOLVIA_1051", "SOLVIA_1125", "SOLVIA_1219", "SOLVIA_1315",
                "SOLVIA_1343", "SOLVIA_1363", "SOLVIA_1564", "SOLVIA_1587", "SOLVIA_1663", "SOLVIA_1680",
                "SOLVIA_1684", "SOLVIA_1743", "SOLVIA_1844", "SOLVIA_1941", "SOLVIA_1952", "SOLVIA_1962",
                "SOLVIA_2012", "SOLVIA_1115", "SOLVIA_794"]
PMF_PEPTIDES = PMF_TOXIC + PMF_NONTOXIC

# Use PMF peptides if running PMF study
if "pmf_study" in config.get("targets", []):
    PEPTIDES = PMF_PEPTIDES

# Create log directory
os.makedirs("logs", exist_ok=True)

# Default parameters
OCCUPANCIES = ["low"]  # For testing, just use low occupancy
PRODUCTION_TIME = config["simulation"]["production"]["time"]

# Determine run number
# Can be overridden with: snakemake --config run=2
if "run" in config:
    RUNS = [str(config["run"])]
else:
    # Auto-detect next run number
    def get_next_run_number(peptide):
        existing_runs = glob.glob(f"simulations/{peptide.lower()}_run_*")
        if not existing_runs:
            return 1
        run_numbers = []
        for run_dir in existing_runs:
            try:
                run_num = int(run_dir.split('_run_')[-1])
                run_numbers.append(run_num)
            except ValueError:
                continue
        return max(run_numbers) + 1 if run_numbers else 1
    
    # For now, use the same run number for all peptides
    RUNS = [str(get_next_run_number(PEPTIDES[0]))]

print(f"Processing run {RUNS[0]} for peptides: {', '.join(PEPTIDES)}")

# Rule 1: Setup run directory structure
rule setup_run:
    input:
        fasta = lambda wildcards: f"{FASTA_DIR}/{wildcards.peptide}.fasta"
    output:
        metadata = "simulations/{peptide}_run_{run}/metadata.yaml"
    params:
        script = "scripts/universal/01_setup_run.py"
    log:
        "logs/{peptide}_run_{run}_setup.log"
    wildcard_constraints:
        run = r'\d+'
    run:
        import subprocess
        
        # Run setup script
        result = subprocess.run(
            [sys.executable, params.script, input.fasta],
            capture_output=True,
            text=True
        )
        
        # Write log
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Setup failed for {wildcards.peptide}")
        
        # Handle directory naming
        peptide_lower = wildcards.peptide.lower()
        expected_dir = f"simulations/{peptide_lower}_run_1"
        target_dir = f"simulations/{peptide_lower}_run_{wildcards.run}"
        
        # Rename if needed
        if os.path.exists(expected_dir) and not os.path.exists(target_dir):
            os.rename(expected_dir, target_dir)
        
        # Create symlink for uppercase
        upper_dir = f"simulations/{wildcards.peptide}_run_{wildcards.run}"
        if not os.path.exists(upper_dir) and target_dir != upper_dir:
            os.symlink(os.path.basename(target_dir), upper_dir)

# Rule 2: Run ColabFold
rule colabfold:
    input:
        metadata = "simulations/{peptide}_run_{run}/metadata.yaml"
    output:
        selection = "simulations/{peptide}_run_{run}/colabfold/model_selection.yaml"
    params:
        script = "scripts/universal/02_run_colabfold.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}"
    threads: 4
    log:
        "logs/{peptide}_run_{run}_colabfold.log"
    run:
        import subprocess
        
        print(f"Running ColabFold for {wildcards.peptide}... This may take 10-30 minutes.")
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"ColabFold failed for {wildcards.peptide}")
        
        # Create symlink for uppercase if needed
        lower_output = f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}/colabfold/model_selection.yaml"
        if os.path.exists(lower_output) and not os.path.exists(output.selection):
            os.makedirs(os.path.dirname(output.selection), exist_ok=True)
            os.symlink(
                f"../../{wildcards.peptide.lower()}_run_{wildcards.run}/colabfold/model_selection.yaml",
                output.selection
            )

# Rule 3: Coarse-graining with Martinize2
rule coarse_grain:
    input:
        selection = "simulations/{peptide}_run_{run}/colabfold/model_selection.yaml"
    output:
        cg_pdb = "simulations/{peptide}_run_{run}/cg_pdb/{peptide}_cg.pdb",
        cg_itp = "simulations/{peptide}_run_{run}/cg_pdb/{peptide}.itp"
    params:
        script = "scripts/universal/03_coarse_grain.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}"
    log:
        "logs/{peptide}_run_{run}_coarse_grain.log"
    run:
        import subprocess
        
        print(f"Running Martinize2 for {wildcards.peptide}...")
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Coarse-graining failed for {wildcards.peptide}")
        
        # Create symlinks for outputs
        real_cg_pdb = f"{params.run_dir}/cg_pdb/{wildcards.peptide}_cg.pdb"
        real_cg_itp = f"{params.run_dir}/cg_pdb/{wildcards.peptide}.itp"
        
        for real_file, output_file in [(real_cg_pdb, output.cg_pdb), (real_cg_itp, output.cg_itp)]:
            if os.path.exists(real_file) and not os.path.exists(output_file):
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                rel_path = os.path.relpath(real_file, os.path.dirname(output_file))
                os.symlink(rel_path, output_file)

# Rule 4: Build membrane template
rule build_membrane:
    input:
        cg_pdb = "simulations/{peptide}_run_{run}/cg_pdb/{peptide}_cg.pdb"
    output:
        membrane_gro = "simulations/{peptide}_run_{run}/membrane_template/membrane_template.gro",
        membrane_top = "simulations/{peptide}_run_{run}/membrane_template/membrane_template.top"
    params:
        script = "scripts/universal/04_build_membrane.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}"
    log:
        "logs/{peptide}_run_{run}_membrane.log"
    run:
        import subprocess
        import shutil
        
        print(f"Building membrane for {wildcards.peptide}...")
        
        # Create necessary directories
        os.makedirs(os.path.join(params.run_dir, "membrane_template"), exist_ok=True)
        os.makedirs(os.path.join(params.run_dir, "logs"), exist_ok=True)
        
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Membrane building failed for {wildcards.peptide}")

# Rule 5: Insert peptides into membrane
rule insert_peptides:
    input:
        cg_pdb = "simulations/{peptide}_run_{run}/cg_pdb/{peptide}_cg.pdb",
        membrane_gro = "simulations/{peptide}_run_{run}/membrane_template/membrane_template.gro"
    output:
        system_gro = "simulations/{peptide}_run_{run}/system/system_{occupancy}.gro",
        system_top = "simulations/{peptide}_run_{run}/system/system_{occupancy}.top"
    params:
        script = "scripts/universal/05_insert_peptides.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}"
    log:
        "logs/{peptide}_run_{run}_insert_{occupancy}.log"
    run:
        import subprocess
        
        print(f"Inserting {wildcards.peptide} peptides ({wildcards.occupancy} occupancy)...")
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir, "--occupancy", wildcards.occupancy],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Peptide insertion failed for {wildcards.peptide}")
        
        # Handle symlinks for output files
        for suffix in ['.gro', '.top']:
            real_file = f"{params.run_dir}/system/system_{wildcards.occupancy}{suffix}"
            output_file = output.system_gro if suffix == '.gro' else output.system_top
            
            if os.path.exists(real_file) and not os.path.exists(output_file):
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                rel_path = os.path.relpath(real_file, os.path.dirname(output_file))
                os.symlink(rel_path, output_file)

# Rule 6: Equilibration (EM, NVT, NPT)
rule equilibrate:
    input:
        system_gro = "simulations/{peptide}_run_{run}/system/system_{occupancy}.gro",
        system_top = "simulations/{peptide}_run_{run}/system/system_{occupancy}.top"
    output:
        npt_gro = "simulations/{peptide}_run_{run}/equilibration/npt/npt_{occupancy}.gro",
        summary = "simulations/{peptide}_run_{run}/equilibration/summary_{occupancy}.yaml"
    params:
        script = "scripts/universal/06_equilibrate.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}"
    log:
        "logs/{peptide}_run_{run}_equilibrate_{occupancy}.log"
    threads: 4
    run:
        import subprocess
        
        print(f"Starting equilibration for {wildcards.peptide} ({wildcards.occupancy} occupancy)...")
        print("  1. Energy minimization (with position restraints)")
        print("  2. NVT equilibration (50 ps)")
        print("  3. NPT equilibration (100 ps)")
        
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir, "--occupancy", wildcards.occupancy],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Equilibration failed for {wildcards.peptide}")
        
        # Create copies for occupancy-specific output names
        import shutil
        npt_src = os.path.join(params.run_dir, "equilibration/npt/npt.gro")
        summary_src = os.path.join(params.run_dir, "equilibration/summary.yaml")
        
        if os.path.exists(npt_src) and not os.path.exists(output.npt_gro):
            os.makedirs(os.path.dirname(output.npt_gro), exist_ok=True)
            shutil.copy(npt_src, output.npt_gro)
        if os.path.exists(summary_src) and not os.path.exists(output.summary):
            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            shutil.copy(summary_src, output.summary)

# Rule 7: Production simulation
rule production:
    input:
        npt_gro = "simulations/{peptide}_run_{run}/equilibration/npt/npt_{occupancy}.gro",
        system_top = "simulations/{peptide}_run_{run}/system/system_{occupancy}.top"
    output:
        tpr = "simulations/{peptide}_run_{run}/production/{occupancy}/production.tpr",
        xtc = "simulations/{peptide}_run_{run}/production/{occupancy}/production.xtc",
        edr = "simulations/{peptide}_run_{run}/production/{occupancy}/production.edr",
        gro = "simulations/{peptide}_run_{run}/production/{occupancy}/production.gro",
        cpt = "simulations/{peptide}_run_{run}/production/{occupancy}/production.cpt"
    params:
        script = "scripts/universal/07_run_production.py",
        run_dir = lambda wildcards: f"simulations/{wildcards.peptide.lower()}_run_{wildcards.run}",
        time = PRODUCTION_TIME
    log:
        "logs/{peptide}_run_{run}_production_{occupancy}.log"
    threads: 4
    run:
        import subprocess
        
        print(f"Starting {params.time} ns production simulation for {wildcards.peptide} ({wildcards.occupancy} occupancy)...")
        print("This will take several hours to days depending on system size and hardware.")
        
        result = subprocess.run(
            [sys.executable, params.script, params.run_dir, 
             "--occupancy", wildcards.occupancy,
             "--time", str(params.time)],
            capture_output=True,
            text=True
        )
        
        with open(log[0], 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n")
                f.write(result.stderr)
        
        if result.returncode != 0:
            raise RuntimeError(f"Production simulation failed for {wildcards.peptide}")

# Rule 8: Setup PMF umbrella windows (for PMF validation study only)
rule setup_pmf_windows:
    input:
        equilibrated = "simulations/{peptide}_run_{run}/equilibration/npt/npt_{occupancy}.gro",
        peptide_gro = "simulations/{peptide}_run_{run}/cg_pdb/{peptide}.gro",
        topology = "simulations/{peptide}_run_{run}/system/system_{occupancy}.top"
    output:
        flag = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/setup_complete.flag"
    params:
        script = "scripts/pmf/01_setup_umbrella_windows.py"
    wildcard_constraints:
        run = r'\d+',
        replicate = r'[12]'
    log:
        "logs/{peptide}_run_{run}_pmf_setup_rep{replicate}_{occupancy}.log"
    shell:
        """
        python {params.script} simulations/{wildcards.peptide}_run_{wildcards.run} \
            --replicate {wildcards.replicate} > {log} 2>&1 && \
        touch {output.flag}
        """

# Rule 9: Run PMF umbrella sampling
rule run_pmf_windows:
    input:
        setup_flag = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/setup_complete.flag"
    output:
        complete_flag = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/windows_complete.flag"
    params:
        script = "scripts/pmf/02_run_umbrella_windows.py",
        parallel = 4
    wildcard_constraints:
        run = r'\d+',
        replicate = r'[12]'
    log:
        "logs/{peptide}_run_{run}_pmf_windows_rep{replicate}.log"
    threads: 4
    shell:
        """
        python {params.script} simulations/{wildcards.peptide}_run_{wildcards.run} \
            --replicate {wildcards.replicate} --parallel {params.parallel} > {log} 2>&1 && \
        touch {output.complete_flag}
        """

# Rule 10: Analyze PMF results
rule analyze_pmf:
    input:
        windows_flag = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/windows_complete.flag"
    output:
        results = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/pmf_results.yaml",
        plot = "simulations/{peptide}_run_{run}/pmf/rep{replicate}/pmf_plot.png"
    params:
        script = "scripts/pmf/03_analyze_pmf.py"
    wildcard_constraints:
        run = r'\d+',
        replicate = r'[12]'
    log:
        "logs/{peptide}_run_{run}_pmf_analysis_rep{replicate}.log"
    shell:
        """
        python {params.script} simulations/{wildcards.peptide}_run_{wildcards.run} \
            --replicate {wildcards.replicate} > {log} 2>&1
        """

# Test rules with explicit run numbers
rule test_setup:
    input:
        expand("simulations/{peptide}_run_{run}/metadata.yaml", 
               peptide=PEPTIDES, run=RUNS)

rule test_colabfold:
    input:
        expand("simulations/{peptide}_run_{run}/colabfold/model_selection.yaml", 
               peptide=PEPTIDES, run=RUNS)

rule test_coarse_grain:
    input:
        expand("simulations/{peptide}_run_{run}/cg_pdb/{peptide}_cg.pdb", 
               peptide=PEPTIDES, run=RUNS)

rule test_membrane:
    input:
        expand("simulations/{peptide}_run_{run}/membrane_template/membrane_template.gro", 
               peptide=PEPTIDES, run=RUNS)

rule test_insert:
    input:
        expand("simulations/{peptide}_run_{run}/system/system_{occupancy}.gro", 
               peptide=PEPTIDES, run=RUNS, occupancy=OCCUPANCIES)

rule test_equilibrate:
    input:
        expand("simulations/{peptide}_run_{run}/equilibration/npt/npt_{occupancy}.gro",
               peptide=PEPTIDES, run=RUNS, occupancy=OCCUPANCIES)

rule test_production:
    input:
        expand("simulations/{peptide}_run_{run}/production/{occupancy}/production.xtc",
               peptide=PEPTIDES, run=RUNS, occupancy=OCCUPANCIES)

rule test_pmf:
    input:
        expand("simulations/{peptide}_run_{run}/pmf/rep{replicate}/pmf_results.yaml",
               peptide=PEPTIDES, run=RUNS, replicate=[1, 2], occupancy=OCCUPANCIES[0])

# Quick test (up to equilibration)
rule test_quick:
    input:
        rules.test_equilibrate.input

# Full pipeline (including production)
rule all:
    input:
        rules.test_production.input

# PMF validation study pipeline
rule pmf_study:
    input:
        expand("simulations/{peptide}_run_{run}/pmf/rep{replicate}/pmf_results.yaml",
               peptide=PMF_PEPTIDES if 'PMF_PEPTIDES' in globals() else [],
               run=RUNS, replicate=[1, 2], occupancy=OCCUPANCIES[0])

# Clean rule for specific run
rule clean_run:
    params:
        run = lambda wildcards: RUNS[0]
    shell:
        """
        echo "Cleaning run {params.run} temporary files..."
        find simulations -name "*_run_{params.run}" -type d -exec find {{}} -name "\\#*\\#" -delete \\;
        find simulations -name "*_run_{params.run}" -type d -exec find {{}} -name "*.log.[0-9]*" -delete \\;
        echo "Done."
        """