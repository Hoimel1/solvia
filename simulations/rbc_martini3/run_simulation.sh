#!/bin/bash
#SBATCH --job-name=rbc_m3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1

# Energy minimization
echo "Starting energy minimization..."
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em

# Check if minimization succeeded
if [ ! -f em.gro ]; then
    echo "Energy minimization failed!"
    exit 1
fi

# NVT equilibration
echo "Starting NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm nvt

# NPT equilibration
echo "Starting NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p system.top -o npt.tpr -maxwarn 2
gmx mdrun -v -deffnm npt

# Production run
echo "Starting production run..."
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p system.top -o prod.tpr -maxwarn 2
gmx mdrun -v -deffnm prod

echo "Simulation complete!"
