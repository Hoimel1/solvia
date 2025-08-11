# Force Field Files Required

This directory should contain the Martini force field files needed for the MD simulations.

## Required Files for Martini 3

Please upload the following files to `force_fields/martini3/`:

1. **martini_v3.0.0.itp** - Main force field definitions
2. **martini_v3.0.0_solvents_v1.itp** - Water and ion topologies
3. **martini_v3.0.0_phospholipids_v1.itp** - Lipid topologies (for membrane)
4. **martini_v3.0.0_ions_v1.itp** - Ion parameters (if separate)
5. **martini.itp** - Link file (if used by martinize2)

## Optional Files

- Additional lipid types if using mixed membranes
- Modified parameters if using custom force fields
- Martini 2 files (in `martini2/` directory) if needed for comparison

## Download Sources

Official Martini force field files can be obtained from:
- http://cgmartini.nl/
- The Martini force field repository

## Usage

Once these files are in place:
1. The solvation step will work properly
2. Membrane building can be completed
3. Ion addition will function correctly

The simulation setup scripts will automatically look for these files and include them in the topology.
