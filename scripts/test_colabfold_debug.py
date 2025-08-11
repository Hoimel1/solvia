#!/usr/bin/env python3
"""Debug test for ColabFold subprocess call"""

import subprocess
import os
from pathlib import Path

# Set up paths
fasta_path = Path("data/raw/fasta_split/SOLVIA_1.fasta")
output_path = Path("test_debug_output")
colabfold_path = Path("localcolabfold/localcolabfold/colabfold-conda")

# Set up environment
env = os.environ.copy()
env['PATH'] = f"{colabfold_path}/bin:{env['PATH']}"

# ColabFold command
cmd = [
    'colabfold_batch',
    '--num-models', '1',
    '--num-recycle', '3',
    str(fasta_path),
    str(output_path)
]

print(f"Running command: {' '.join(cmd)}")
print(f"PATH includes: {colabfold_path}/bin")

try:
    result = subprocess.run(
        cmd, 
        capture_output=True, 
        text=True, 
        env=env,
        timeout=300
    )
    
    print(f"Return code: {result.returncode}")
    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")
        
except Exception as e:
    print(f"Exception type: {type(e)}")
    print(f"Exception: {e}")
    import traceback
    traceback.print_exc()
