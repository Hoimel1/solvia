#!/bin/bash
# Production MD
echo "Running production MD..."
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -r npt.gro -t npt.cpt -maxwarn 1
gmx mdrun -v -deffnm md -nt 8
