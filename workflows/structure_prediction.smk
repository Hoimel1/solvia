"""
Structure prediction workflow using ColabFold
Based on implementation guide (Section 10)
"""

rule predict_structure:
    input:
        fasta = DATA_DIR / "input" / "sequences" / "{amp_id}.fasta"
    output:
        pdb = DATA_DIR / "intermediate" / "structures" / "{amp_id}" / "ranked_0.pdb",
        plddt = DATA_DIR / "intermediate" / "structures" / "{amp_id}" / "plddt_scores.json"
    params:
        output_dir = DATA_DIR / "intermediate" / "structures" / "{amp_id}",
        num_seeds = config["structure_prediction"]["num_seeds"],
        model_preset = config["structure_prediction"]["model_preset"]
    container:
        config["containers"]["colabfold"]
    threads: 4
    resources:
        gpu = 1,
        mem_gb = 16
    log:
        "logs/structure_prediction/{amp_id}.log"
    shell:
        """
        colabfold_batch \
            {input.fasta} \
            {params.output_dir} \
            --num-seeds {params.num_seeds} \
            --model-preset {params.model_preset} \
            --use-gpu \
            > {log} 2>&1
            
        # Extract pLDDT scores
        python -c "
import json
import numpy as np
from Bio import PDB

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure('pdb', '{output.pdb}')
plddt_scores = []

for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.name == 'CA':
                    plddt_scores.append(atom.bfactor)

mean_plddt = np.mean(plddt_scores)
with open('{output.plddt}', 'w') as f:
    json.dump({{'mean_plddt': mean_plddt, 'scores': plddt_scores}}, f)
"
        """

