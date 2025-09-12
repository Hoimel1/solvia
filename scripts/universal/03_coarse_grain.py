#!/usr/bin/env python3
"""
SOLVIA Peptide Coarse-Graining
Converts atomistic structures to Martini 3 coarse-grained representation
Handles all ITP processing and symlink creation automatically
"""

import os
import sys
import yaml
import argparse
import subprocess
import shutil
import re
from pathlib import Path
import numpy as np

def load_config(config_file=None):
    """Load SOLVIA configuration"""
    config_dir = Path(__file__).parent.parent.parent / "config"
    
    if config_file:
        # Use specified config file
        if os.path.isabs(config_file):
            config_path = Path(config_file)
        else:
            config_path = config_dir / config_file
    else:
        # Try available config files in order of preference
        config_files = [
            "config.yaml",
            "pmf_standard_config.yaml"
        ]
        
        config_path = None
        for config_file in config_files:
            potential_path = config_dir / config_file
            if potential_path.exists():
                config_path = potential_path
                break
    
    if not config_path or not config_path.exists():
        print("Warning: No config file found, using defaults")
        return {
            'coarse_graining': {
                'c_terminal': 'NH2-ter',
                'force_field': 'martini3001',
                'elastic_network': True
            }
        }
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    # Add default coarse-graining values if missing
    if 'coarse_graining' not in config:
        config['coarse_graining'] = {}
    
    cg_defaults = {
        'c_terminal': 'NH2-ter',
        'force_field': 'martini3001',
        'force_field_auto': False,  # avoid mixing M2/M3
        'elastic_network': False,   # safer default off for peptides
        'elastic_k': 300,           # if enabled, softer defaults
        'elastic_rc': 0.6,
        'use_dssp': True,           # request DSSP if image supports it
        'docker_service': 'martinize',
        'posres_k': 300.0           # CG position restraints (kJ/mol/nm^2)
    }
    
    for key, default_value in cg_defaults.items():
        if key not in config['coarse_graining']:
            config['coarse_graining'][key] = default_value
            
    print(f"Loaded configuration from: {config_path}")
    return config

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def get_best_structure(run_dir):
    """Get best structure from ColabFold or other source"""
    # Check for model selection file
    selection_file = os.path.join(run_dir, "colabfold", "model_selection.yaml")
    if not os.path.exists(selection_file):
        print("Error: No model selection found. Run ColabFold first.")
        sys.exit(1)
    
    with open(selection_file, 'r') as f:
        selection = yaml.safe_load(f)
    
    # Use best_model (filename) or best_pdb (full path)
    if 'best_model' in selection and selection['best_model']:
        pdb_path = os.path.join(run_dir, "colabfold", selection['best_model'])
    elif 'best_pdb' in selection:
        # best_pdb might be a full path or relative path
        if os.path.isabs(selection['best_pdb']):
            pdb_path = selection['best_pdb']
        else:
            pdb_path = selection['best_pdb']
    else:
        print("Error: No best model found in selection file")
        sys.exit(1)
    
    # Ensure best_model.pdb is an actual file inside the colabfold dir (not a host-absolute symlink)
    best_model_link = os.path.join(run_dir, "colabfold", "best_model.pdb")
    try:
        if os.path.islink(best_model_link) or os.path.exists(best_model_link):
            os.remove(best_model_link)
        shutil.copy2(pdb_path, best_model_link)
    except Exception as e:
        print(f"Warning: Could not stage best_model.pdb: {e}")
    
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found: {pdb_path}")
        sys.exit(1)
    
    return pdb_path, selection.get('best_plddt', 80)

def extract_peptide_id(run_dir):
    """Extract peptide ID from run directory name"""
    run_dir_name = os.path.basename(run_dir)  # e.g., "solvia_1_run_1"
    # Extract peptide_id: solvia_1_run_1 -> SOLVIA_1
    parts = run_dir_name.split('_')
    if len(parts) >= 2:
        peptide_id = '_'.join(parts[:2]).upper()  # "SOLVIA_1"
    else:
        peptide_id = "PEPTIDE"
    return peptide_id

def run_martinize_docker(run_dir, pdb_path, config):
    """Run Martinize via Docker and handle all ITP processing"""
    print(f"\n{'='*50}")
    print("Running Martinize via Docker...")
    print(f"{'='*50}")
    
    # Extract configuration
    cg_config = config.get('coarse_graining', {})
    c_terminal = cg_config.get('c_terminal', 'NH2-ter')
    # Force-field selection (adaptive for short/aromatic peptides)
    def _read_fasta(fa_path: Path):
        try:
            with open(fa_path, 'r') as f:
                lines = [l.strip() for l in f if l.strip()]
            return ''.join(l for l in lines if not l.startswith('>'))
        except Exception:
            return ''

    def _select_optimal_force_field(seq: str) -> tuple[str, str]:
        # Keep the CG ecosystem consistent: prefer Martini 3
        return 'martini3001', 'martini3 (consistent ecosystem)'

    force_field = cg_config.get('force_field', 'martini3001')
    if cg_config.get('force_field_auto', False):
        seq = _read_fasta(Path(run_dir) / 'input' / 'peptide.fasta')
        ff_sel, rationale = _select_optimal_force_field(seq)
        if ff_sel != force_field:
            print(f"  Force-field auto-select: {force_field} -> {ff_sel} ({rationale})")
            force_field = ff_sel
    
    # Get peptide ID
    peptide_id = extract_peptide_id(run_dir)
    
    # Get relative path from project root
    project_root = Path(__file__).parent.parent.parent
    rel_run_dir = os.path.relpath(run_dir, project_root)
    
    # Create output directory
    output_dir = os.path.join(run_dir, "coarse_grain")
    os.makedirs(output_dir, exist_ok=True)
    # Pre-clean stale outputs that could block martinize writes (e.g., symlinks)
    try:
        out = Path(output_dir)
        for name in (f"{peptide_id}.itp", f"{peptide_id}_0.itp"):
            p = out / name
            if p.exists() or p.is_symlink():
                p.unlink()
        # Remove editor/Gromacs backups
        for p in out.iterdir():
            n = p.name
            if n.startswith('#') and n.endswith('#'):
                try: p.unlink()
                except Exception: pass
    except Exception:
        pass
    
    # Build Docker command
    service = str(cg_config.get('docker_service', 'martinize'))
    cmd = [
        'docker', 'compose', 'run', '--rm',
        '-w', f'/work/{rel_run_dir}/coarse_grain',  # ensure all emitted files land in coarse_grain
        service,
        '-f', f'/work/{rel_run_dir}/colabfold/best_model.pdb',
        '-ff', force_field,
        '-x', f'/work/{rel_run_dir}/coarse_grain/{peptide_id}_cg.pdb',
        '-o', f'/work/{rel_run_dir}/coarse_grain/{peptide_id}.itp',
        '-name', peptide_id,
        '-cter', c_terminal,  # C-terminal amidation
        '-maxwarn', '1',
        '-cys', 'auto'         # Auto-detect disulfide bonds
    ]
    if bool(cg_config.get('use_dssp', True)):
        cmd.append('-dssp')
    # Encourage a self-contained ITP with backbone position restraints header
    cmd.extend(['-p', 'backbone', '-pf', '1000'])
    
    # Add elastic network only if enabled in config
    if cg_config.get('elastic_network', False):
        ef = str(int(cg_config.get('elastic_k', 300)))
        eu = str(float(cg_config.get('elastic_rc', 0.6)))
        cmd.extend(['-elastic', '-ef', ef, '-eu', eu])
        print("  Elastic network: enabled")
    else:
        print("  Elastic network: disabled")
    
    print(f"Peptide ID: {peptide_id}")
    print(f"Command: {' '.join(cmd)}")
    
    # Run Martinize via Docker
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(project_root))
    
    if result.returncode != 0:
        print("\nError running Martinize:")
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        sys.exit(1)
    
    print("\nMartinize output:")
    print(result.stdout)
    
    # Handle the ITP files created by martinize
    handle_itp_files(output_dir, peptide_id)
    
    # Create symlinks
    create_symlinks(output_dir, peptide_id)
    
    # Create position restraint file
    create_position_restraints(output_dir, peptide_id, kpos=float(cg_config.get('posres_k', 300.0)))
    
    return peptide_id

def _pdb_read_coords(pdb_file):
    """Read PDB ATOM/HETATM coords and atom names/res ids.

    Returns: list of (resid, atomname, (x,y,z)) in nm.
    """
    recs = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                try:
                    resid = int(line[22:26])
                    atom = line[12:16].strip()
                    x = float(line[30:38]) / 10.0
                    y = float(line[38:46]) / 10.0
                    z = float(line[46:54]) / 10.0
                    recs.append((resid, atom, np.array([x,y,z], dtype=float)))
                except Exception:
                    continue
    except Exception:
        pass
    return recs

def _rg(coords):
    if coords is None or len(coords) == 0:
        return None
    arr = np.vstack(coords)
    com = arr.mean(axis=0)
    return float(np.sqrt(((arr - com)**2).sum(axis=1).mean()))

def _clash_fraction(coords, cutoff_nm=0.38):
    n = len(coords)
    if n < 2:
        return 0.0
    arr = np.vstack(coords)
    clashes = 0
    total = 0
    for i in range(n-1):
        for j in range(i+1, n):
            total += 1
            if np.linalg.norm(arr[i]-arr[j]) < float(cutoff_nm):
                clashes += 1
    return float(clashes/total) if total else 0.0

def _backbone_distances_cg(cg_recs):
    bb = [(res, xyz) for (res, atom, xyz) in cg_recs if atom.upper() == 'BB']
    bb = sorted(bb, key=lambda t: t[0])
    dists = []
    for i in range(len(bb)-1):
        d = float(np.linalg.norm(bb[i+1][1]-bb[i][1]))
        dists.append(d)
    return dists

def _merge_cg_config(cfg: dict) -> dict:
    """Return a unified coarse_graining config from either top-level or pmf.coarse_graining."""
    if cfg is None:
        return {}
    if 'coarse_graining' in cfg:
        base = dict(cfg.get('coarse_graining') or {})
    else:
        base = {}
    pmf_cg = ((cfg.get('pmf') or {}).get('coarse_graining') or {})
    # pmf.coarse_graining keys override base
    base.update(pmf_cg)
    return base

def validate_cg_structure(run_dir, peptide_id, pdb_path, config=None):
    """Validate CG structure geometry and mapping fidelity; write YAML report.

    Checks (configurable via coarse_graining.validation):
      - CG clash fraction (pair distance < clash_cutoff_nm; default 0.38 nm)
      - CG backbone BB distances within [bb_dist_min_nm, bb_dist_max_nm] (default [0.25, 0.50] nm)
      - Optional: Rg deviation AA vs CG (< rg_dev_pct_max), disabled by default without AA→CG mapping
    """
    import yaml as _yaml
    out_dir = Path(run_dir) / 'coarse_grain'
    cg_pdb = out_dir / f"{peptide_id}_cg.pdb"
    report = {'passed': True, 'metrics': {}, 'thresholds': {}}
    # thresholds (configurable)
    cg_cfg = _merge_cg_config(config or {})
    vcfg = (cg_cfg.get('validation') or {})
    thr = {
        'clash_frac_max': float(vcfg.get('clash_frac_max', 0.02)),
        'clash_cutoff_nm': float(vcfg.get('clash_cutoff_nm', 0.38)),
        'bb_dist_min_nm': float(vcfg.get('bb_dist_min_nm', 0.25)),
        'bb_dist_max_nm': float(vcfg.get('bb_dist_max_nm', 0.50)),
        'rg_dev_pct_max': float(vcfg.get('rg_dev_pct_max', 20.0)),
        'rg_check_enabled': bool(vcfg.get('rg_check_enabled', False)),
    }
    report['thresholds'] = thr
    # Read coords
    cg_recs = _pdb_read_coords(str(cg_pdb))
    aa_recs = _pdb_read_coords(str(pdb_path))
    cg_coords = [xyz for (_r, _a, xyz) in cg_recs]
    aa_coords = [xyz for (_r, _a, xyz) in aa_recs]
    # Clash fraction
    cf = _clash_fraction(cg_coords, cutoff_nm=thr['clash_cutoff_nm'])
    report['metrics']['clash_fraction'] = float(cf)
    if cf > thr['clash_frac_max']:
        report['passed'] = False
    # BB distances
    bb_d = _backbone_distances_cg(cg_recs)
    if bb_d:
        d_min = float(np.min(bb_d)); d_max = float(np.max(bb_d)); d_mean = float(np.mean(bb_d))
    else:
        d_min = d_max = d_mean = None
    report['metrics']['bb_dist_min'] = d_min
    report['metrics']['bb_dist_max'] = d_max
    report['metrics']['bb_dist_mean'] = d_mean
    if d_min is not None and (d_min < thr['bb_dist_min_nm'] or d_max > thr['bb_dist_max_nm']):
        report['passed'] = False
    # Rg mapping fidelity (disabled by default; AA→CG mapping recommended for strict checks)
    rg_cg = _rg(cg_coords); rg_aa = _rg(aa_coords)
    report['metrics']['rg_cg_nm'] = rg_cg
    report['metrics']['rg_aa_nm'] = rg_aa
    if thr['rg_check_enabled'] and (rg_cg is not None) and (rg_aa is not None) and rg_aa > 1e-12:
        dev = abs(rg_cg - rg_aa) / rg_aa * 100.0
        report['metrics']['rg_deviation_pct'] = float(dev)
        if dev > thr['rg_dev_pct_max']:
            report['passed'] = False
    else:
        report['metrics']['rg_deviation_pct'] = None
    # Write report
    try:
        with open(out_dir / 'cg_validation.yaml', 'w') as f:
            _yaml.dump(report, f, default_flow_style=False)
    except Exception:
        pass
    # Console summary
    print("\nCG validation:")
    print(f"  Clash fraction: {cf:.4f} (≤ {thr['clash_frac_max']}, cutoff {thr['clash_cutoff_nm']:.2f} nm)")
    if d_min is not None:
        print(f"  BB distance min/mean/max: {d_min:.3f}/{d_mean:.3f}/{d_max:.3f} nm (in [{thr['bb_dist_min_nm']:.2f}, {thr['bb_dist_max_nm']:.2f}])")
    if report['metrics']['rg_deviation_pct'] is not None:
        print(f"  Rg deviation: {report['metrics']['rg_deviation_pct']:.1f}% (≤ {thr['rg_dev_pct_max']}%)")
    print(f"  Passed: {report['passed']}")
    if not report['passed']:
        print("  Recommendation: Inspect CG mapping; consider disabling elastic network or adjusting martinize parameters for peptides.")
    return report

def handle_itp_files(output_dir, peptide_id):
    """Process martinize outputs to a single, clean ITP and centralize it.

    - Prefer a self-contained file without includes.
    - If martinize produced a wrapper {peptide_id}.itp that includes {peptide_id}_0.itp,
      replace {peptide_id}.itp with the contents of {peptide_id}_0.itp.
    - Strip any include lines to global force fields (e.g., "martini.itp").
    - Normalize [ moleculetype ] name to {peptide_id}.
    - Copy final ITP to project force_fields/{peptide_id}.itp.
    """
    print("\nProcessing ITP files (flatten + centralize)...")
    out = Path(output_dir)
    project_root = Path(__file__).parent.parent.parent
    ff_dir = project_root / 'force_fields'
    ff_dir.mkdir(exist_ok=True)

    mol_itp = out / f"{peptide_id}.itp"
    mol0_itp = out / f"{peptide_id}_0.itp"
    # If older martinize dumped _0.itp into project root (container CWD), relocate it
    stray_root = project_root / f"{peptide_id}_0.itp"
    if (not mol0_itp.exists()) and stray_root.exists():
        try:
            stray_root.replace(mol0_itp)
            print(f"  ✓ Moved stray {stray_root.name} into {out.name}/")
        except Exception as e:
            print(f"  Warning: could not move stray {stray_root.name}: {e}")

    # If martinize emitted only _0.itp, adopt it as the main ITP
    if mol0_itp.exists() and (not mol_itp.exists() or mol0_itp.stat().st_size > mol_itp.stat().st_size):
        try:
            content = mol0_itp.read_text()
            mol_itp.write_text(content)
            print(f"  ✓ Adopted {mol0_itp.name} as {mol_itp.name}")
        except Exception as e:
            print(f"  Warning: could not adopt {mol0_itp.name}: {e}")

    # If wrapper detected (includes or system/molecules), try to replace with a real ITP
    if mol_itp.exists():
        text = mol_itp.read_text()
        wrapper_like = ('#include' in text and f'{peptide_id}_0.itp' in text) or ('\n[ system ]' in text)
        if wrapper_like and mol0_itp.exists():
            try:
                text = mol0_itp.read_text()
                mol_itp.write_text(text)
                print(f"  ✓ Flattened wrapper using {mol0_itp.name}")
            except Exception as e:
                print(f"  Warning: could not flatten wrapper: {e}")

    # Report found ITPs
    for p in sorted(out.glob(f"{peptide_id}*.itp")):
        try:
            size = p.stat().st_size
        except Exception:
            size = 0
        print(f"  found: {p.name} ({size} B)")

    # Sanitize: strip includes to force field aliases and normalize moleculetype name
    if mol_itp.exists():
        try:
            lines = mol_itp.read_text().splitlines()
            out_lines = []
            in_mt = False
            mt_fixed = False
            for line in lines:
                s = line.strip()
                # Drop includes to global FF or local wrappers
                if s.startswith('#include') and ('martini' in s or f'{peptide_id}_0.itp' in s):
                    continue
                if s.lower().startswith('[ moleculetype'):
                    in_mt = True
                    out_lines.append(line)
                    continue
                if in_mt:
                    if not s or s.startswith(';'):
                        out_lines.append(line)
                        continue
                    cols = s.split()
                    if cols:
                        cols[0] = peptide_id
                        out_lines.append(' '.join(cols))
                        in_mt = False
                        mt_fixed = True
                        continue
                # Filter out accidental top-level sections if present in wrapper
                if s.startswith('[') and any(h in s.lower() for h in ('[ system', '[ molecules')):
                    # Skip wrapper-only sections
                    continue
                out_lines.append(line)
            mol_itp.write_text('\n'.join(out_lines) + '\n')
            if mt_fixed:
                print(f"  ✓ Normalized moleculetype name to {peptide_id}")
        except Exception as e:
            print(f"  Warning: could not sanitize ITP: {e}")

    # Copy final ITP to force_fields and create convenience symlink back
    if mol_itp.exists():
        try:
            target = ff_dir / f"{peptide_id}.itp"
            # Write/overwrite centralized copy
            target.write_text(mol_itp.read_text())
            print(f"  ✓ Copied to force_fields: {target}")
        except Exception as e:
            print(f"  Warning: could not copy ITP to force_fields: {e}")

    # Ensure only a single authoritative ITP exists to avoid double counting
    # 1) Remove local {peptide}_0.itp if present
    try:
        if mol0_itp.exists():
            mol0_itp.unlink()
            print(f"  ✓ Removed duplicate {mol0_itp.name}")
    except Exception as e:
        print(f"  Warning: could not remove {mol0_itp.name}: {e}")

    # 2) Keep local {peptide}.itp as a regular file for now (martinize may rerun)
    #    We will replace it with a symlink in a finalization step after validation.

    # 3) Remove editor/GROMACS-style backup files (#...#) in coarse_grain
    try:
        for p in out.iterdir():
            n = p.name
            if n.startswith('#') and n.endswith('#'):
                try:
                    p.unlink()
                except Exception:
                    pass
    except Exception:
        pass

def create_symlinks(output_dir, peptide_id):
    """Create symlinks for backward compatibility"""
    print("\nCreating symlinks...")
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    # Prefer centralized ITP in force_fields for all runs
    project_root = Path(__file__).parent.parent.parent
    ff_itp = project_root / 'force_fields' / f'{peptide_id}.itp'

    # Create symlinks to the actual files
    symlink_map = {
        'peptide_cg.pdb': f'{peptide_id}_cg.pdb',
        'peptide_cg.itp': str(ff_itp) if ff_itp.exists() else f'{peptide_id}.itp'
    }
    
    for link, target in symlink_map.items():
        if os.path.exists(link) or os.path.islink(link):
            os.remove(link)
        if os.path.exists(target):
            os.symlink(target, link)
            print(f"  ✓ {link} -> {target}")
    
    # Return to original directory
    os.chdir(original_dir)

def finalize_itp(output_dir, peptide_id):
    """After all martinize runs/validation, make coarse_grain/{peptide}.itp a symlink to force_fields.

    Also ensure duplicates and backups are removed.
    """
    out = Path(output_dir)
    project_root = Path(__file__).parent.parent.parent
    ff_itp = project_root / 'force_fields' / f'{peptide_id}.itp'
    local_itp = out / f'{peptide_id}.itp'
    # Clean backups
    try:
        for p in out.iterdir():
            n = p.name
            if n.startswith('#') and n.endswith('#'):
                try: p.unlink()
                except Exception: pass
    except Exception:
        pass
    # Replace local ITP with symlink if central exists
    try:
        if ff_itp.exists():
            if local_itp.exists() or local_itp.is_symlink():
                try: local_itp.unlink()
                except Exception: pass
            # Create a relative symlink so it resolves both on host and in container (/work mount)
            rel = os.path.relpath(ff_itp, out)
            local_itp.symlink_to(rel)
            print(f"  ✓ Finalized: {local_itp.name} -> {rel}")
    except Exception as e:
        print(f"  Warning: finalize link failed: {e}")

def create_position_restraints(output_dir, peptide_id, kpos=300.0):
    """Create position restraint file (CG)."""
    print("\nGenerating position restraints...")
    
    posre_file = os.path.join(output_dir, "peptide_posre.itp")
    itp_file = os.path.join(output_dir, f"{peptide_id}.itp")
    
    if not os.path.exists(itp_file):
        print(f"  Warning: ITP file not found: {itp_file}")
        return
    
    atom_ids = []
    in_atoms = False
    try:
        with open(itp_file, 'r') as itp:
            for line in itp:
                s = line.strip()
                if not s or s.startswith(';'):
                    continue
                if s.startswith('['):
                    in_atoms = s.lower().startswith('[ atoms')
                    continue
                if in_atoms:
                    cols = s.split()
                    # Typical columns: id type resnr resid atomname cgnr charge
                    if len(cols) >= 6:
                        atomname = cols[4]
                        if atomname == 'BB':
                            atom_ids.append(int(cols[0]))
    except Exception:
        pass
    with open(posre_file, 'w') as f:
        f.write(f"; Position restraint file for {peptide_id}\n")
        f.write("[ position_restraints ]\n")
        f.write("; atom  type      fx      fy      fz\n")
        for aid in atom_ids:
            f.write(f"{aid:6d}     1  {float(kpos):.1f}  {float(kpos):.1f}  {float(kpos):.1f}\n")
    
    print(f"  ✓ Created {posre_file}")

def print_summary(run_dir, peptide_id, config, metadata=None):
    """Print summary of coarse-graining results"""
    rel_run_dir = os.path.relpath(run_dir, Path(__file__).parent.parent.parent)
    c_terminal = config.get('coarse_graining', {}).get('c_terminal', 'NH2-ter')
    elastic_network = config.get('coarse_graining', {}).get('elastic_network', False)
    
    print(f"\n{'='*50}")
    print("✅ Coarse-graining complete!")
    print(f"{'='*50}")
    
    if metadata:
        print(f"Run information:")
        print(f"  Peptide ID: {metadata.get('peptide_id', peptide_id)}")
        print(f"  Run: {metadata.get('run_name', 'Unknown')}")
        print(f"  Created: {metadata.get('created', 'Unknown')}")
        print()
    
    print(f"Output files:")
    print(f"  PDB: {rel_run_dir}/coarse_grain/{peptide_id}_cg.pdb")
    print(f"  ITP: {rel_run_dir}/coarse_grain/{peptide_id}.itp (single file, no double counting)")
    print(f"  Position restraints: {rel_run_dir}/coarse_grain/peptide_posre.itp")
    print(f"\nSymlinks created:")
    print(f"  peptide_cg.pdb -> {peptide_id}_cg.pdb")
    print(f"  peptide_cg.itp -> {peptide_id}.itp")
    # GRO not generated by default; convert explicitly if needed in downstream steps
    print(f"\nMartinize settings:")
    print(f"  Molecule name: {peptide_id}")
    print(f"  C-terminal: {c_terminal} (amidated)")
    # N-terminal: rely on martinize2 defaults unless configured explicitly
    print(f"  Elastic network: {'enabled' if elastic_network else 'disabled'}")
    print(f"{'='*50}")

def main():
    parser = argparse.ArgumentParser(
        description="Coarse-grain peptide structure using Martinize via Docker"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory containing ColabFold output"
    )
    parser.add_argument(
        "--config",
        help="Configuration file to use (default: search for config.yaml or pmf_standard_config.yaml)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run even if output exists"
    )
    
    args = parser.parse_args()
    
    # Convert to absolute path
    run_dir = os.path.abspath(args.run_dir)
    
    # Check run directory
    if not os.path.exists(run_dir):
        print(f"Error: Run directory not found: {run_dir}")
        sys.exit(1)
    
    # Check if already done (unless forced)
    peptide_id = extract_peptide_id(run_dir)
    output_dir = os.path.join(run_dir, "coarse_grain")
    output_pdb = os.path.join(output_dir, f"{peptide_id}_cg.pdb")
    output_itp = os.path.join(output_dir, f"{peptide_id}.itp")
    
    if not args.force and os.path.exists(output_pdb) and os.path.exists(output_itp):
        print("Coarse-graining already completed. Use --force to re-run.")
        try:
            metadata = load_run_metadata(run_dir)
        except:
            metadata = None
        print_summary(run_dir, peptide_id, load_config(args.config), metadata)
        return
    
    # Load configuration
    config = load_config(args.config)
    try:
        metadata = load_run_metadata(run_dir)
    except Exception:
        metadata = None
    
    # Get best structure
    pdb_path, plddt = get_best_structure(run_dir)
    print(f"Using best structure: {os.path.basename(pdb_path)} (pLDDT: {plddt:.1f})")
    
    # Run martinize via Docker
    peptide_id = run_martinize_docker(run_dir, pdb_path, config)
    
    # Validation of CG structure with optional auto-retry
    try:
        report = validate_cg_structure(run_dir, peptide_id, pdb_path, config)
        cg_cfg = _merge_cg_config(config)
        auto = ((cg_cfg.get('validation') or {}).get('auto_retry') or {})
        if (not report.get('passed')) and bool(auto.get('enabled', True)):
            steps = auto.get('steps', ['disable_en'])
            retries = int(auto.get('max_retries', 1))
            attempted = 0
            for step in steps:
                if attempted >= retries:
                    break
                if step == 'disable_en':
                    print("Auto-retry: disabling elastic network and re-running martinize...")
                    # Clone config and turn off EN
                    cfg2 = dict(config)
                    # normalize nesting
                    if 'coarse_graining' in cfg2:
                        cfg2['coarse_graining'] = dict(cfg2['coarse_graining'])
                        cfg2['coarse_graining']['elastic_network'] = False
                    else:
                        cfg2.setdefault('coarse_graining', {})['elastic_network'] = False
                    run_martinize_docker(run_dir, pdb_path, cfg2)
                    report = validate_cg_structure(run_dir, peptide_id, pdb_path, cfg2)
                    attempted += 1
                    if report.get('passed'):
                        break
    except Exception as e:
        print(f"Warning: CG validation failed: {e}")

    # Finalize ITP centralization/symlinks now that martinize runs are done
    try:
        finalize_itp(os.path.join(run_dir, 'coarse_grain'), peptide_id)
    except Exception as e:
        print(f"Warning: finalize step failed: {e}")
    
    # Print summary
    print_summary(run_dir, peptide_id, config, metadata)

if __name__ == "__main__":
    main()
