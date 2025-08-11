#!/usr/bin/env python3
"""
Rename molecule_0.itp files to peptide-specific names (e.g., SOLVIA_1.itp)
"""

import os
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def rename_itp_files():
    """Rename all molecule_0.itp files to peptide_id.itp"""
    
    topology_dir = Path("data/processed/topologies")
    
    for peptide_dir in topology_dir.iterdir():
        if peptide_dir.is_dir():
            peptide_id = peptide_dir.name
            old_itp = peptide_dir / "molecule_0.itp"
            new_itp = peptide_dir / f"{peptide_id}.itp"
            
            if old_itp.exists():
                # Read the file and update the molecule name
                with open(old_itp, 'r') as f:
                    content = f.read()
                
                # Replace molecule_0 with peptide_id in the content
                content = content.replace('molecule_0', peptide_id)
                
                # Write to new file
                with open(new_itp, 'w') as f:
                    f.write(content)
                
                # Update topol.top to reference the new file
                topol_file = peptide_dir / "topol.top"
                if topol_file.exists():
                    with open(topol_file, 'r') as f:
                        topol_content = f.read()
                    
                    topol_content = topol_content.replace(
                        '#include "molecule_0.itp"', 
                        f'#include "{peptide_id}.itp"'
                    )
                    topol_content = topol_content.replace('molecule_0', peptide_id)
                    
                    with open(topol_file, 'w') as f:
                        f.write(topol_content)
                
                # Remove old file
                old_itp.unlink()
                
                logger.info(f"Renamed {old_itp} to {new_itp}")

if __name__ == "__main__":
    rename_itp_files()
    logger.info("ITP renaming complete!")
