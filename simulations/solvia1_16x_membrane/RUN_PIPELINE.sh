#!/bin/bash
# Master script to run complete simulation

echo "SOLVIA Membrane Simulation Pipeline"
echo "==================================="
echo "Peptide: SOLVIA_1
echo "16 copies on POPC membrane"
echo ""
echo "NOTE: This requires a properly built membrane system!"
echo "See MEMBRANE_BUILDING_INSTRUCTIONS.txt"
echo ""
echo "Once you have a complete system.gro file with membrane+peptides+water+ions:"
echo ""
echo "./01_minimize.sh"
echo "./02_equilibrate.sh"
echo "./03_production.sh"
echo "./04_analyze.sh"
