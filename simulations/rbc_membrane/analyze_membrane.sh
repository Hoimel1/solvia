#!/bin/bash
# Analysis for RBC membrane-peptide system

echo "RBC Membrane-Peptide Analysis"
echo "============================="

# 1. System setup info
echo "Analyzing system composition..."
echo -e "0\n" | gmx check -f md.xtc 2>&1 | grep -E "Coords|Box"

# 2. Peptide-membrane distance
echo "Calculating peptide-membrane distances..."
echo -e "1\n2\n" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi

# 3. Membrane thickness
echo "Analyzing membrane thickness..."
# Would need gmx density or similar tool

# 4. Peptide orientation
echo "Analyzing peptide orientation..."
echo -e "1\n" | gmx gangle -f md.xtc -s md.tpr -g1 vector -g2 z -oav peptide_angle.xvg

# 5. Lipid order parameters
echo "Note: For detailed lipid analysis, use MDAnalysis or similar tools"

# 6. Create visualization files
echo "Creating visualization files..."
echo -e "1\n2\n0\n" | gmx trjconv -f md.xtc -s md.tpr -o viz_peptide_membrane.xtc -pbc mol -ur compact -dt 1000

# Final structure
echo -e "0\n" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo ""
echo "Key files for visualization:"
echo "- final_structure.pdb: Final frame in PDB format"
echo "- viz_peptide_membrane.xtc: Trajectory for VMD (every 1 ns)"
echo ""
echo "VMD visualization commands:"
echo "vmd final_structure.pdb viz_peptide_membrane.xtc"
echo ""
echo "In VMD:"
echo "- Peptides: 'resname SOLVIA_1'"
echo "- Membrane: 'resname POPC POPE POPS CHOL'"
echo "- Use 'Drawing Method' -> 'QuickSurf' for membrane"
