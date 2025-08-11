#!/bin/bash
# Run membrane simulation with 16 peptides

echo "Starting membrane simulation with 16x SOLVIA_1"
date

# Note: Full pipeline requires:
# 1. INSANE or packmol to build proper membrane
# 2. Combining peptides with membrane
# 3. Solvation with Martini water
# 4. Ion addition

# For now, we have a simplified setup
echo "This is a simplified setup. For production:"
echo "1. Use INSANE to build POPC bilayer"
echo "2. Position peptides using gmx insert-molecules"
echo "3. Solvate with Martini W beads"
echo "4. Add ions with gmx genion"

# Example commands (to be run after proper setup):
# gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr
# gmx mdrun -v -deffnm em

# gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro
# gmx mdrun -v -deffnm nvt

# gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt
# gmx mdrun -v -deffnm npt

# gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -r npt.gro -t npt.cpt
# gmx mdrun -v -deffnm md

echo "Setup complete. Ready for manual membrane building."
