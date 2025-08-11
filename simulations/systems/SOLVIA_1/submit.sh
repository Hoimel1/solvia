#!/bin/bash
#SBATCH --job-name=SOLVIA_1_md
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --output=SOLVIA_1_%j.out
#SBATCH --error=SOLVIA_1_%j.err

# Load GROMACS module (adjust as needed)
module load gromacs/2021.4

# Set number of OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Starting simulation for SOLVIA_1"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on: $HOSTNAME"
echo "Start time: $(date)"

# Energy minimization
gmx grompp -f em.mdp -c SOLVIA_1_box.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v

# NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v

# NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v

# Production MD
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md -v

echo "Simulation completed"
echo "End time: $(date)"
