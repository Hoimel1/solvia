#!/bin/bash
# Energy minimization
echo "Running energy minimization..."
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em -nt 8
