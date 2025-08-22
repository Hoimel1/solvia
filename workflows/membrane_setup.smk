"""
Membrane setup workflow using INSANE
Based on biophysical principles (Section 4.1)
"""

rule build_membrane:
    output:
        gro = DATA_DIR / "intermediate" / "membranes" / "membrane_base.gro",
        top = DATA_DIR / "intermediate" / "membranes" / "membrane_base.top"
    params:
        box = config["membrane"]["box_size"],
        outer = config["membrane"]["composition"]["outer_leaflet"],
        inner = config["membrane"]["composition"]["inner_leaflet"],
        salt = config["membrane"]["salt_concentration"]
    container:
        config["containers"]["insane"]
    log:
        "logs/membrane_setup/base_membrane.log"
    shell:
        """
        # Build asymmetric RBC-like membrane
        # Outer leaflet: POPC:PSM:CHOL
        # Inner leaflet: POPE:POPS:CHOL
        
        python /opt/insane/insane.py \
            -o {output.gro} \
            -p {output.top} \
            -box {params.box[0]} {params.box[1]} {params.box[2]} \
            -sol W \
            -salt {params.salt} \
            -u "POPC:{params.outer[POPC]} PSM:{params.outer[PSM]} CHOL:{params.outer[CHOL]}" \
            -l "POPE:{params.inner[POPE]} POPS:{params.inner[POPS]} CHOL:{params.inner[CHOL]}" \
            -asym 0.2 \
            > {log} 2>&1
        """

rule insert_peptides:
    input:
        membrane_gro = DATA_DIR / "intermediate" / "membranes" / "membrane_base.gro",
        membrane_top = DATA_DIR / "intermediate" / "membranes" / "membrane_base.top",
        peptide_pdb = DATA_DIR / "intermediate" / "coarse_grained" / "{amp_id}" / "cg.pdb",
        peptide_itp = DATA_DIR / "intermediate" / "coarse_grained" / "{amp_id}" / "molecule.itp"
    output:
        system_gro = DATA_DIR / "intermediate" / "systems" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "system.gro",
        system_top = DATA_DIR / "intermediate" / "systems" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "system.top"
    params:
        n_peptides = lambda wildcards: config["simulation"]["occupancies"][wildcards.occupancy],
        seed = lambda wildcards: int(wildcards.replicate) * 42
    log:
        "logs/membrane_setup/{amp_id}_{occupancy}_rep{replicate}.log"
    script:
        "../src/simulation/insert_peptides.py"

