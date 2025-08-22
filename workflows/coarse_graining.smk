"""
Coarse-graining workflow using martinize2
Based on core scientific concept (Section 4.1)
"""

rule coarse_grain:
    input:
        pdb = DATA_DIR / "intermediate" / "structures" / "{amp_id}" / "ranked_0.pdb",
        plddt = DATA_DIR / "intermediate" / "structures" / "{amp_id}" / "plddt_scores.json"
    output:
        cg_pdb = DATA_DIR / "intermediate" / "coarse_grained" / "{amp_id}" / "cg.pdb",
        topology = DATA_DIR / "intermediate" / "coarse_grained" / "{amp_id}" / "molecule.itp"
    params:
        forcefield = config["coarse_graining"]["forcefield"],
        plddt_threshold = config["structure_prediction"]["plddt_threshold"]
    container:
        config["containers"]["martinize"]
    log:
        "logs/coarse_graining/{amp_id}.log"
    shell:
        """
        # Check if IDP flag is needed based on pLDDT
        mean_plddt=$(python -c "import json; print(json.load(open('{input.plddt}'))['mean_plddt'])")
        
        if (( $(echo "$mean_plddt < {params.plddt_threshold}" | bc -l) )); then
            IDP_FLAG="--martini3-idp"
        else
            IDP_FLAG=""
        fi
        
        martinize2 \
            -f {input.pdb} \
            -o {output.topology} \
            -x {output.cg_pdb} \
            -ff {params.forcefield} \
            $IDP_FLAG \
            -dssp \
            > {log} 2>&1
        """

