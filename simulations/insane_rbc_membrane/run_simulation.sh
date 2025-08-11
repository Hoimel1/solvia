#!/bin/bash
# Run RBC membrane simulation for SOLVIA_1

echo "Starting RBC membrane simulation pipeline..."

# 1. Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8

# 2. NVT equilibration  
echo "Step 2: NVT equilibration (1 ns)"
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -maxwarn 10
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8

# 3. NPT equilibration
echo "Step 3: NPT equilibration (5 ns)"  
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -maxwarn 10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8

# 4. Production run
echo "Step 4: Production run (500 ns)"
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt -maxwarn 10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 8

echo "Simulation complete!"

# 5. Basic analysis
echo "Step 5: Analysis"
echo -e "1\n1\n" | gmx mindist -f md.xtc -s md.tpr -o peptide_membrane_dist.xvg -pi
echo -e "0\n" | gmx trjconv -f md.xtc -s md.tpr -o viz_full.xtc -pbc mol -ur compact
echo -e "0\n" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

echo "Analysis complete!"
echo "Visualize with: vmd final_structure.pdb viz_full.xtc"
