#!/bin/bash
# Martini 3 demo simulation

echo "Running Martini 3 demo for SOLVIA_1..."

# 1. Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 5
gmx mdrun -deffnm em -v

# 2. NVT equilibration
echo "Step 2: NVT equilibration (100 ps)"
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5
gmx mdrun -deffnm nvt -v

# 3. Production
echo "Step 3: Production run (1 ns)"
gmx grompp -f md.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -deffnm md -v

# 4. Analysis
echo "Step 4: Creating visualization files"
echo "0" | gmx trjconv -f md.xtc -s md.tpr -o traj_centered.xtc -pbc mol -center
echo "0" | gmx trjconv -f md.gro -s md.tpr -o final.pdb -pbc mol -center

echo ""
echo "Demo complete!"
echo "View trajectory: vmd final.pdb traj_centered.xtc"
echo ""
echo "This was a simple demo. For membrane simulations:"
echo "1. Use CHARMM-GUI Martini Maker"
echo "2. Select Martini 3 force field"
echo "3. Build your RBC membrane composition"
