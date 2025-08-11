#!/usr/bin/env python3
"""
Build RBC membrane using insane with standard lipids.
Since PSM might not be available, we'll use DPPC as a substitute.
"""

import os
import subprocess
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SimpleRBCMembraneBuilder:
    def __init__(self, peptide_id, output_dir):
        self.peptide_id = peptide_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use the CG PDB array
        self.peptide_array = Path(f"simulations/membrane_16x/{peptide_id}_16x.pdb")
        
        # Check if we already have the 16x array
        if not self.peptide_array.exists():
            logger.info(f"Creating peptide array first...")
            subprocess.run([
                "python", "scripts/setup_16peptide_membrane.py",
                "--peptide-id", peptide_id,
                "--num-peptides", "16"
            ])
            
        # Copy necessary files
        self.prepare_files()
        
    def prepare_files(self):
        """Copy all necessary files to working directory."""
        # Copy peptide array  
        shutil.copy(self.peptide_array, self.output_dir / "peptides_16x.pdb")
        
        # Copy force field files
        ff_dir = Path("force_fields/martini3")
        for ff_file in ff_dir.glob("*.itp"):
            shutil.copy(ff_file, self.output_dir / ff_file.name)
        
        # Copy peptide topology with correct name
        top_src = Path(f"data/processed/topologies/{self.peptide_id}/{self.peptide_id}.itp")
        if top_src.exists():
            shutil.copy(top_src, self.output_dir / f"{self.peptide_id}.itp")
    
    def build_rbc_membrane(self):
        """Build the RBC membrane using insane with standard lipids."""
        logger.info("Building RBC membrane with insane (standard lipids)")
        
        # Calculate lipid numbers for 12x12 nm membrane
        # Outer leaflet: 45% POPC, 10% DPPC (substitute for PSM), 45% CHOL
        n_popc = 97
        n_dppc = 22  # Substitute for PSM
        n_chol_outer = 97
        
        # Inner leaflet: 45% POPE, 15% POPS, 40% CHOL  
        n_pope = 97
        n_pops = 32
        n_chol_inner = 86
        
        # Total cholesterol
        n_chol_total = n_chol_outer + n_chol_inner
        
        # First, try a simple symmetric membrane
        cmd = [
            "insane",
            "-f", "peptides_16x.pdb",           # Input peptides
            "-o", "system.gro",                 # Output coordinates
            "-p", "system_insane.top",          # Output topology
            "-x", "12",                         # Box X dimension
            "-y", "12",                         # Box Y dimension  
            "-z", "15",                         # Box Z dimension
            "-center",                          # Center the system
            "-orient",                          # Orient proteins in membrane
            "-sol", "W",                        # Use Martini water
            "-salt", "0.15",                    # 150 mM salt
            # Lipid composition - simplified
            "-l", f"POPC:{n_popc + 50}",       # More POPC
            "-l", f"POPE:{n_pope}",
            "-l", f"POPS:{n_pops}",
            "-l", f"DPPC:{n_dppc}",             # DPPC instead of PSM
            "-l", f"CHOL:{n_chol_total}",
            "-pbc", "rectangular",              # PBC type
            "-dm", "3"                          # Position proteins 3nm from membrane center
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
        # Run insane
        result = subprocess.run(
            cmd,
            cwd=self.output_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            logger.info("Insane completed successfully!")
            logger.info("Output summary:")
            # Extract key information from output
            output_lines = result.stdout.split('\n') + result.stderr.split('\n')
            for line in output_lines:
                if any(keyword in line for keyword in ['Adding', 'beads', 'lipids', 'Protein', 'oriented']):
                    logger.info(f"  {line.strip()}")
                    
            # Save the output for debugging
            with open(self.output_dir / "insane_output.log", 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
        else:
            logger.error("Insane failed")
            logger.error(f"STDERR: {result.stderr}")
            logger.error(f"STDOUT: {result.stdout}")
            
            # Save error log
            with open(self.output_dir / "insane_error.log", 'w') as f:
                f.write("Command: " + ' '.join(cmd) + "\n\n")
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n") 
                f.write(result.stderr)
            
            return False
            
        # Fix the topology file
        self.fix_topology()
        
        return True
    
    def fix_topology(self):
        """Fix the topology file generated by insane."""
        logger.info("Fixing topology file for Martini 3")
        
        # Read insane-generated topology
        insane_top = self.output_dir / "system_insane.top"
        if not insane_top.exists():
            logger.error("Insane topology not found!")
            return
            
        with open(insane_top, 'r') as f:
            insane_content = f.read()
        
        # Create proper Martini 3 topology
        topology = f"""; RBC-like membrane system with 16x {self.peptide_id}
; Generated with insane
; Note: Using DPPC as substitute for PSM

; Include force field parameters
#include "martini_v3.0.0.itp"

; Include lipid topologies
#include "martini_v3.0.0_phospholipids_v1.itp"
#include "martini_v3.0_sterols_v1.0.itp"

; Include solvent and ion topologies
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"

; Include peptide topology
#include "{self.peptide_id}.itp"

[ system ]
; Name
RBC-like membrane with 16x {self.peptide_id}

"""
        
        # Extract molecule counts from insane topology
        in_molecules = False
        molecule_lines = []
        
        for line in insane_content.split('\n'):
            if "[ molecules ]" in line:
                in_molecules = True
                continue
            elif in_molecules and line.strip() and not line.startswith(';'):
                # Clean up molecule names
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    mol_count = parts[1]
                    
                    # Fix peptide name
                    if mol_name in ['Protein', 'protein', 'molecule_0']:
                        mol_name = self.peptide_id
                    
                    molecule_lines.append(f"{mol_name:<16} {mol_count}")
        
        # Add molecules section
        topology += "[ molecules ]\n"
        topology += "; Compound        #mols\n"
        topology += f"; Target composition:\n"
        topology += f"; Outer: ~45% POPC, 10% DPPC (PSM substitute), 45% CHOL\n"
        topology += f"; Inner: ~45% POPE, 15% POPS, 40% CHOL\n\n"
        
        for line in molecule_lines:
            topology += line + "\n"
        
        # Write fixed topology
        with open(self.output_dir / "system.top", 'w') as f:
            f.write(topology)
        
        logger.info("Topology fixed and saved as system.top")
        
        # Show summary
        logger.info("System composition:")
        for line in molecule_lines:
            logger.info(f"  {line}")
    
    def create_mdp_files(self):
        """Copy MDP files from RBC membrane directory."""
        src_dir = Path("simulations/rbc_membrane")
        
        for mdp in ['em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']:
            src = src_dir / mdp
            if src.exists():
                shutil.copy(src, self.output_dir / mdp)
        
        logger.info("MDP files copied")
    
    def create_run_script(self):
        """Create a script to run the complete simulation."""
        
        run_script = f"""#!/bin/bash
# Run RBC-like membrane simulation for {self.peptide_id}

echo "Starting RBC membrane simulation pipeline..."
echo "System: 16x {self.peptide_id} on RBC-like membrane"
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
"""
        
        script_path = self.output_dir / "run_simulation.sh"
        with open(script_path, 'w') as f:
            f.write(run_script)
        os.chmod(script_path, 0o755)
        
        logger.info(f"Created run script: {script_path}")
    
    def run(self):
        """Run the complete pipeline."""
        logger.info(f"Building RBC-like membrane for {self.peptide_id} using insane")
        logger.info("Note: Using DPPC as substitute for PSM")
        
        # Build membrane
        if not self.build_rbc_membrane():
            logger.error("Failed to build membrane")
            return False
            
        # Copy MDP files
        self.create_mdp_files()
        
        # Create run script
        self.create_run_script()
        
        logger.info(f"\n{'='*60}")
        logger.info(f"RBC-like membrane setup complete in {self.output_dir}")
        logger.info(f"{'='*60}")
        logger.info("\nNext steps:")
        logger.info(f"1. cd {self.output_dir}")
        logger.info("2. Check system.gro in VMD")
        logger.info("3. ./run_simulation.sh")
        logger.info(f"{'='*60}\n")
        
        return True

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build RBC-like membrane with insane")
    parser.add_argument('--peptide-id', default='SOLVIA_1', help='Peptide ID')
    parser.add_argument('--output-dir', default='simulations/rbc_simple', help='Output directory')
    args = parser.parse_args()
    
    builder = SimpleRBCMembraneBuilder(args.peptide_id, args.output_dir)
    builder.run()

if __name__ == "__main__":
    main()
