#!/bin/bash
# Basic analysis

echo "Extracting trajectory..."
# Extract protein and membrane for visualization
echo -e "1\n13\n" | gmx trjconv -f md.xtc -s md.tpr -o md_protein_membrane.xtc -pbc mol -ur compact

echo "Calculating RMSD..."
echo -e "4\n4\n" | gmx rms -f md.xtc -s md.tpr -o rmsd.xvg

echo "Calculating RMSF..."
echo -e "4\n" | gmx rmsf -f md.xtc -s md.tpr -o rmsf.xvg -res

echo "Creating PDB for VMD..."
echo -e "0\n" | gmx trjconv -f md.gro -s md.tpr -o final_frame.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo "Load final_frame.pdb and md_protein_membrane.xtc in VMD for visualization"
