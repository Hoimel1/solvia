#!/bin/bash
# Test Martini 3 system

echo "Testing Martini 3 peptide system..."

# 1. Test without water first
echo "Step 1: Testing peptide topology..."
gmx grompp -f em_martini3.mdp -c peptides.gro -p test_nowat.top -o test.tpr -maxwarn 5
if [ $? -eq 0 ]; then
    echo "Topology test passed!"
else
    echo "Topology error - check SOLVIA_1_16x.itp"
    exit 1
fi

echo ""
echo "Peptide system is ready for membrane insertion."
echo "See MEMBRANE_INSTRUCTIONS.txt for next steps."
