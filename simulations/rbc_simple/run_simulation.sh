#!/bin/bash
# Run RBC-like membrane simulation for SOLVIA_1

echo "Starting RBC membrane simulation pipeline..."
echo "System: 16x SOLVIA_1 on RBC-like membrane"
echo "Note: Using DPPC as substitute for PSM"

# 1. Energy minimization
echo "========================================="
echo "Step 1: Energy minimization"
echo "========================================="
gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 10
if [ $? -ne 0 ]; then
    echo "Error in grompp for minimization"
    exit 1
fi

gmx mdrun -deffnm em -ntmpi 1 -ntomp 8
if [ $? -ne 0 ]; then
    echo "Error in energy minimization"
    exit 1
fi

# Check if minimization converged
tail -n 20 em.log | grep -E "Steepest|converged"

# 2. Generate index file for temperature coupling groups
echo "Creating index file..."
cat << EOF | gmx make_ndx -f em.gro -o index.ndx
q
EOF

# 3. NVT equilibration (1 ns)
echo "========================================="
echo "Step 2: NVT equilibration (1 ns)"
echo "========================================="
gmx grompp -f nvt.mdp -c em.gro -p system.top -o nvt.tpr -r em.gro -n index.ndx -maxwarn 10
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8

# 4. NPT equilibration (5 ns)
echo "========================================="
echo "Step 3: NPT equilibration (5 ns)"  
echo "========================================="
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr -r nvt.gro -t nvt.cpt -n index.ndx -maxwarn 10
gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8

# 5. Short production run (10 ns for testing)
echo "========================================="
echo "Step 4: Short production run (10 ns)"
echo "========================================="
# Create modified mdp for short run
sed 's/nsteps.*=.*/nsteps = 400000  ; 10 ns/' md.mdp > md_short.mdp

gmx grompp -f md_short.mdp -c npt.gro -p system.top -o md.tpr -t npt.cpt -n index.ndx -maxwarn 10
gmx mdrun -deffnm md -ntmpi 1 -ntomp 8

echo "========================================="
echo "Simulation complete!"
echo "========================================="

# 6. Basic analysis
echo "Step 5: Basic analysis"

# Create visualization trajectory (every 100 ps)
echo "0" | gmx trjconv -f md.xtc -s md.tpr -o viz.xtc -pbc mol -ur compact -dt 100

# Final structure in PDB format
echo "0" | gmx trjconv -f md.gro -s md.tpr -o final_structure.pdb -pbc mol -ur compact

echo ""
echo "Analysis files created:"
echo "- viz.xtc: Trajectory for visualization (every 100 ps)"
echo "- final_structure.pdb: Final structure"
echo ""
echo "To visualize in VMD:"
echo "vmd final_structure.pdb viz.xtc"
echo ""
echo "To run feature extraction:"
echo "python ../../scripts/extract_features.py"
