#!/usr/bin/env python3
"""
SOLVIA RBC Membrane Template Builder
Creates evidence-based asymmetric RBC membrane using INSANE
"""

import os
import sys
import yaml
import argparse
import shutil
import subprocess
from pathlib import Path

def load_config():
    """Load SOLVIA configuration"""
    config_path = Path(__file__).parent.parent.parent / "config" / "solvia_config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def build_membrane_template(config):
    """Build RBC membrane template using INSANE"""
    # Get membrane parameters from config
    membrane_config = config['membrane']
    
    # Build INSANE command: always use container PATH (ignore host venv)
    insane_cmd = shutil.which("insane") or "insane"
    
    # Format box size as single string
    box_str = f"{membrane_config['box_size'][0]},{membrane_config['box_size'][1]},{membrane_config['box_size'][2]}"
    
    cmd = [
        insane_cmd,
        "-o", "membrane_template.gro",
        "-p", "membrane_template.top",
        "-box", box_str,
        "-sol", membrane_config['solvent'],
        "-salt", str(membrane_config['salt_concentration'])
    ]
    
    # Add upper leaflet lipids
    upper_str = []
    for lipid, fraction in membrane_config['upper_leaflet'].items():
        cmd.extend(["-u", f"{lipid}:{fraction}"])
        upper_str.append(f"{lipid}:{fraction}")
    
    # Add lower leaflet lipids
    lower_str = []
    for lipid, fraction in membrane_config['lower_leaflet'].items():
        cmd.extend(["-l", f"{lipid}:{fraction}"])
        lower_str.append(f"{lipid}:{fraction}")
    
    return cmd, " ".join(upper_str), " ".join(lower_str)

def run_insane(run_dir):
    """Run INSANE to create membrane template"""
    config = load_config()
    original_dir = os.getcwd()
    
    # Output directory
    if run_dir:
        run_dir_abs = os.path.abspath(run_dir)
        output_dir = os.path.join(run_dir_abs, "membrane_template")
    else:
        # Create global template
        templates_base = os.path.abspath(config['directories']['templates'])
        output_dir = os.path.join(templates_base, "membrane")
        os.makedirs(output_dir, exist_ok=True)
    # Ensure output directory exists for run-specific builds
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if already exists
    gro_file = os.path.join(output_dir, "membrane_template.gro")
    top_file = os.path.join(output_dir, "membrane_template.top")
    
    if os.path.exists(gro_file) and os.path.exists(top_file):
        print("Membrane template already exists.")
        return gro_file, top_file
    
    # Build command
    cmd, upper_str, lower_str = build_membrane_template(config)
    
    # Log file
    if run_dir:
        log_file = os.path.join(run_dir_abs, "logs", "insane_membrane.log")
    else:
        log_file = os.path.join(output_dir, "insane_membrane.log")
    # Ensure logs directory exists for run-specific builds
    if run_dir:
        os.makedirs(os.path.join(run_dir_abs, "logs"), exist_ok=True)
    
    print("Building RBC membrane template with INSANE...")
    print(f"Upper leaflet: {upper_str}")
    print(f"Lower leaflet: {lower_str}")
    print(f"Box size: {' x '.join(map(str, config['membrane']['box_size']))} nm")
    print(f"Salt: {config['membrane']['salt_concentration']} M")
    
    # Change to output directory
    os.chdir(output_dir)
    
    # Run INSANE
    with open(log_file, 'w') as log:
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True
            )
            log.write(result.stdout)
            print("✓ INSANE completed successfully")
        except subprocess.CalledProcessError as e:
            log.write(e.stdout if e.stdout else "")
            print(f"✗ INSANE failed. Check log: {log_file}")
            os.chdir(original_dir)
            sys.exit(1)
    
    os.chdir(original_dir)
    
    # Parse membrane composition from log
    parse_membrane_composition(log_file, output_dir)
    
    # Create README for membrane
    create_membrane_readme(output_dir, config)
    
    return gro_file, top_file

def parse_membrane_composition(log_file, output_dir):
    """Parse actual membrane composition from topology file"""
    top_file = os.path.join(output_dir, "membrane_template.top")
    if not os.path.exists(top_file):
        print("Warning: Topology file not found")
        return
    
    composition = {}
    
    # Read topology file
    with open(top_file, 'r') as f:
        lines = f.readlines()
    
    # Find molecules section
    in_molecules = False
    for line in lines:
        line = line.strip()
        if '[ molecules ]' in line:
            in_molecules = True
            continue
        if in_molecules and line and not line.startswith(';'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    molecule = parts[0]
                    count = int(parts[1])
                    # Handle duplicate CHOL entries
                    if molecule in composition:
                        composition[molecule] += count
                    else:
                        composition[molecule] = count
                except:
                    pass
    
    # Save composition
    with open(os.path.join(output_dir, "composition.yaml"), 'w') as f:
        yaml.dump(composition, f, default_flow_style=False)
    
    # Print summary
    if composition:
        print("\nMembrane composition:")
        total_lipids = sum(v for k, v in composition.items() if k not in ['W', 'NA', 'CL'])
        for mol, count in sorted(composition.items()):
            if mol not in ['W', 'NA', 'CL']:
                print(f"  {mol}: {count} ({count/total_lipids*100:.1f}%)")
            else:
                print(f"  {mol}: {count}")

def create_membrane_readme(output_dir, config):
    """Create README for membrane template"""
    readme_content = f"""# RBC Membrane Template

## Description
Evidence-based asymmetric red blood cell (RBC) membrane model for SOLVIA hemolytic toxicity predictions.

## Composition
### Upper Leaflet (Outer)
- POPC: {config['membrane']['upper_leaflet']['POPC']*100:.0f}%
- PSM: {config['membrane']['upper_leaflet']['PSM']*100:.0f}%
- CHOL: {config['membrane']['upper_leaflet']['CHOL']*100:.0f}%

### Lower Leaflet (Inner)
- POPE: {config['membrane']['lower_leaflet']['POPE']*100:.0f}%
- POPS: {config['membrane']['lower_leaflet']['POPS']*100:.0f}%
- CHOL: {config['membrane']['lower_leaflet']['CHOL']*100:.0f}%

## Parameters
- Box size: {' x '.join(map(str, config['membrane']['box_size']))} nm
- Membrane thickness: ~{config['membrane']['thickness']} nm
- Salt concentration: {config['membrane']['salt_concentration']} M (NaCl)
- Total cholesterol: ~40-45 mol%

## Evidence Base
Based on literature values for human RBC membranes:
- High cholesterol content for membrane stability
- Asymmetric distribution with PS in inner leaflet
- Physiological ionic strength

## Files
- `membrane_template.gro`: Structure file
- `membrane_template.top`: Topology file
- `composition.yaml`: Actual lipid counts
- `insane_membrane.log`: Build log

Generated by SOLVIA membrane builder
"""
    
    with open(os.path.join(output_dir, "README.md"), 'w') as f:
        f.write(readme_content)

def main():
    parser = argparse.ArgumentParser(
        description="Build RBC membrane template for SOLVIA"
    )
    parser.add_argument(
        "run_dir",
        nargs='?',
        help="Run directory (optional - if not provided, creates global template)"
    )
    parser.add_argument(
        "--global-template",
        action="store_true",
        help="Create global membrane template in data/templates/membrane"
    )
    
    args = parser.parse_args()
    
    # Determine where to create template
    if args.global_template or not args.run_dir:
        print("Creating global membrane template...")
        run_dir = None
    else:
        if not os.path.exists(args.run_dir):
            print(f"Error: Run directory not found: {args.run_dir}")
            sys.exit(1)
        run_dir = args.run_dir
    
    # Build membrane
    gro_file, top_file = run_insane(run_dir)
    
    print(f"\n✓ Membrane template created:")
    print(f"  Structure: {gro_file}")
    print(f"  Topology: {top_file}")
    
    if run_dir:
        print(f"\nNext step: Insert peptides into membrane")
        print(f"Command: python 05_insert_peptides.py {run_dir}")

if __name__ == "__main__":
    main()
