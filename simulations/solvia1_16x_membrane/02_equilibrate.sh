#!/bin/bash
# Equilibration protocol for membrane system

echo "Running NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -maxwarn 1
gmx mdrun -v -deffnm nvt -nt 8

echo "Running NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 1
gmx mdrun -v -deffnm npt -nt 8
