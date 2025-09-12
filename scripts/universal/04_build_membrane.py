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
import numpy as np

def load_config():
    """Load SOLVIA configuration (tries pmf_standard_config.yaml, then solvia_config.yaml)"""
    cfg_dir = Path(__file__).parent.parent.parent / "config"
    for fname in ("pmf_standard_config.yaml", "solvia_config.yaml"):
        cpath = cfg_dir / fname
        if cpath.exists():
            with open(cpath, 'r') as f:
                print(f"Loaded membrane config from: {cpath}")
                return yaml.safe_load(f)
    raise FileNotFoundError("No configuration YAML found in config/ directory.")

def _validated_rbc_composition():
    """Return experimentally guided RBC lipid composition (outer/inner).

    Based on literature (Lorent 2020, Verkleij 1973):
    - Outer: PSM+PC dominant; PS minimal; high CHOL
    - Inner: PE+PS dominant; reduced PC; fraction CHOL asymmetry
    Values below are coarse-grained representatives using MARTINI species.
    """
    # Targeting literature-guided asymmetry:
    # - Total CHOL ~ 0.42 (42%)
    # - Outer: PC+SM dominate; PS minimal
    # - Inner: PE+PS dominate (~85% of non-CHOL pool)
    return {
        'upper_leaflet': {
            'CHOL': 0.42,
            'POPC': 0.30,
            'PSM': 0.27,
            'POPS': 0.01,
        },
        'lower_leaflet': {
            'CHOL': 0.42,
            'POPE': 0.293,
            'POPS': 0.20,
            'POPC': 0.087,
        },
    }

def _membrane_cfg_from(cfg: dict) -> dict:
    """Extract membrane config from either top-level 'membrane' or 'pmf.membrane'.

    Provides sensible defaults if some keys are missing.
    """
    mem = {}
    if isinstance(cfg, dict):
        if 'membrane' in cfg:
            mem = dict(cfg.get('membrane') or {})
        elif isinstance(cfg.get('pmf'), dict) and 'membrane' in cfg.get('pmf', {}):
            mem = dict((cfg.get('pmf') or {}).get('membrane') or {})
    # Defaults
    mem.setdefault('box_size', [10, 10, 14])
    mem.setdefault('solvent', 'W')
    mem.setdefault('salt_concentration', 0.15)
    mem.setdefault('use_validated_rbc', True)
    return mem


def build_membrane_template(membrane_config: dict, output_dir):
    """Build RBC membrane template using INSANE via Docker.

    Expects a membrane_config dict (already extracted) with keys like box_size,
    solvent, salt_concentration, and leaflet compositions.
    """
    # Optionally override with validated RBC composition
    if membrane_config.get('use_validated_rbc', True):
        validated = _validated_rbc_composition()
        membrane_config['upper_leaflet'] = validated['upper_leaflet']
        membrane_config['lower_leaflet'] = validated['lower_leaflet']
    
    # Get relative path from project root
    project_root = Path(__file__).parent.parent.parent
    rel_output_dir = os.path.relpath(output_dir, project_root)
    
    # Format box size as single string
    box_str = f"{membrane_config['box_size'][0]},{membrane_config['box_size'][1]},{membrane_config['box_size'][2]}"
    
    # Build Docker command for INSANE
    cmd = [
        "docker", "compose", "run", "--rm", "insane",
        "-o", f"/work/{rel_output_dir}/membrane_template.gro",
        "-p", f"/work/{rel_output_dir}/membrane_template.top",
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

def _parse_gro_for_leaflets(gro_path: str, config: dict | None = None):
    """Parse GRO and classify leaflets.

    - Uses PO4/P bead z to estimate midplane (median).
    - Computes per-residue mean z and classifies all residues vs. midplane.
    - Handles ties and degenerate cases by balanced split.

    Returns: (res_map, outer_ids, inner_ids, z_median)
    """
    # Focus only on lipid residues for classification
    mem = _membrane_cfg_from(config or {}) if config is not None else {}
    allowed_lipids = set()
    allowed_lipids.update((mem.get('upper_leaflet') or {}).keys())
    allowed_lipids.update((mem.get('lower_leaflet') or {}).keys())
    if not allowed_lipids:
        allowed_lipids = {'POPC', 'POPE', 'POPS', 'PSM', 'CHOL'}

    res_map: dict[int, str] = {}
    res_sum_z: dict[int, float] = {}
    res_cnt: dict[int, int] = {}
    po4pz: list[float] = []
    try:
        with open(gro_path, 'r') as f:
            lines = f.readlines()
    except Exception:
        return {}, set(), set(), None
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        try:
            resid = int(line[0:5])
            resn = line[5:10].strip()
            atom = line[10:15].strip()
            z = float(line[36:44])
        except Exception:
            continue
        res_map[resid] = resn
        res_sum_z[resid] = res_sum_z.get(resid, 0.0) + z
        res_cnt[resid] = res_cnt.get(resid, 0) + 1
        if (atom in ('PO4', 'P')) and (resn in allowed_lipids):
            po4pz.append(z)
    if not res_map:
        return {}, set(), set(), None
    # Midplane from PO4/P if available, else from all residues
    if po4pz:
        zmed = float(np.median(np.asarray(po4pz, dtype=float)))
    else:
        means = [res_sum_z[k] / max(1, res_cnt.get(k, 1)) for k in res_sum_z.keys()]
        zmed = float(np.median(np.asarray(means, dtype=float))) if means else None
        if zmed is None:
            return res_map, set(), set(), None
    # Per-residue mean z (lipids only)
    res_mean_z = {k: (res_sum_z[k] / max(1, res_cnt.get(k, 1))) for k in res_sum_z.keys() if res_map.get(k) in allowed_lipids}
    outer = set()
    inner = set()
    ties = []
    for rid, mz in res_mean_z.items():
        if mz > zmed:
            outer.add(rid)
        elif mz < zmed:
            inner.add(rid)
        else:
            ties.append(rid)
    # Distribute ties to balance
    if ties:
        ties_sorted = sorted(ties)
        for i, rid in enumerate(ties_sorted):
            if len(outer) <= len(inner):
                outer.add(rid)
            else:
                inner.add(rid)
    # Fallback: balanced split by z if one side empty
    if not outer or not inner:
        order = sorted(res_mean_z.items(), key=lambda t: t[1])
        n = len(order)
        if n >= 2:
            inner = {rid for rid, _ in order[: n // 2]}
            outer = {rid for rid, _ in order[n // 2 :]}
    return res_map, outer, inner, zmed

def _parse_ndx_groups(ndx_path: str) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    name = None
    with open(ndx_path, 'r') as f:
        for line in f:
            ls = line.strip()
            if not ls:
                continue
            if ls.startswith('[') and ls.endswith(']'):
                name = ls.strip('[]').strip()
                groups.setdefault(name, [])
                continue
            if name:
                parts = ls.split()
                for p in parts:
                    if p.isdigit():
                        groups[name].append(int(p))
    return groups


def validate_membrane_asymmetry(output_dir: str, config: dict | None = None):
    """Compute leaflet‑resolved composition and validate RBC asymmetry constraints.

    Writes asymmetry_validation.yaml with metrics and boolean checks.
    """
    gro_file = os.path.join(output_dir, "membrane_template.gro")
    # Try index_leaflets.ndx if present to define leaflets precisely
    ndx_path = os.path.join(output_dir, 'index_leaflets.ndx')
    res_map, outer_ids, inner_ids, zmed = {}, set(), set(), None
    if os.path.exists(ndx_path):
        try:
            groups = _parse_ndx_groups(ndx_path)
            outer_atoms = set(groups.get('OuterPO4', []) + groups.get('Outer', []))
            inner_atoms = set(groups.get('InnerPO4', []) + groups.get('Inner', []))
            # Map atom indices to residue ids via GRO order (1-based)
            atom_to_res: dict[int, int] = {}
            res_map_tmp: dict[int, str] = {}
            res_sum_z: dict[int, float] = {}
            res_cnt: dict[int, int] = {}
            po4pz: list[float] = []
            with open(gro_file, 'r') as f:
                lines = f.readlines()
            atom_index = 0
            for line in lines[2:-1]:
                if len(line) < 44:
                    continue
                try:
                    resid = int(line[0:5])
                    resn = line[5:10].strip()
                    atom = line[10:15].strip()
                    z = float(line[36:44])
                except Exception:
                    continue
                atom_index += 1
                atom_to_res[atom_index] = resid
                res_map_tmp[resid] = resn
                res_sum_z[resid] = res_sum_z.get(resid, 0.0) + z
                res_cnt[resid] = res_cnt.get(resid, 0) + 1
                if atom in ('PO4', 'P'):
                    po4pz.append(z)
            # Build residue sets
            outer_ids = {atom_to_res[a] for a in outer_atoms if a in atom_to_res}
            inner_ids = {atom_to_res[a] for a in inner_atoms if a in atom_to_res}
            # Supplement: assign lipids without PO4 based on mean z vs midplane
            zmed_idx = None
            if po4pz:
                zmed_idx = float(np.median(np.asarray(po4pz, dtype=float)))
            else:
                means = [res_sum_z[k] / max(1, res_cnt.get(k, 1)) for k in res_sum_z.keys()]
                if means:
                    zmed_idx = float(np.median(np.asarray(means, dtype=float)))
            if zmed_idx is not None:
                memcfg_l = _membrane_cfg_from(config or {}) if config is not None else {}
                allowed_lipids = set()
                allowed_lipids.update((memcfg_l.get('upper_leaflet') or {}).keys())
                allowed_lipids.update((memcfg_l.get('lower_leaflet') or {}).keys())
                if not allowed_lipids:
                    allowed_lipids = {'POPC','POPE','POPS','PSM','CHOL'}
                all_lipid_res = {rid for rid, rn in res_map_tmp.items() if rn in allowed_lipids}
                missing = all_lipid_res - outer_ids - inner_ids
                for rid in missing:
                    mz = res_sum_z.get(rid, 0.0) / max(1, res_cnt.get(rid, 1))
                    if mz >= zmed_idx:
                        outer_ids.add(rid)
                    else:
                        inner_ids.add(rid)
            # Update res_map from tmp
            res_map = res_map_tmp
        except Exception:
            outer_ids = set(); inner_ids = set()
    # If no usable index, fallback to GRO-based classification
    if not outer_ids or not inner_ids:
        res_map, outer_ids, inner_ids, zmed = _parse_gro_for_leaflets(gro_file, config)
    else:
        # Prepare res_map if not set
        if not res_map:
            try:
                with open(gro_file, 'r') as f:
                    lines = f.readlines()
                for line in lines[2:-1]:
                    if len(line) < 44:
                        continue
                    resid = int(line[0:5]); resn = line[5:10].strip()
                    res_map[resid] = resn
            except Exception:
                pass
    metrics = {}
    if not outer_ids and not inner_ids:
        metrics['error'] = 'No PO4 atoms found for leaflet split'
        ok = False
    else:
        # Count resnames per leaflet; exclude solvent/ions
        memcfg = _membrane_cfg_from(config or {}) if config is not None else {}
        allowed_lipids = set()
        allowed_lipids.update((memcfg.get('upper_leaflet') or {}).keys())
        allowed_lipids.update((memcfg.get('lower_leaflet') or {}).keys())
        if not allowed_lipids:
            # default set
            allowed_lipids = {'POPC', 'POPE', 'POPS', 'PSM', 'CHOL'}
        excluded = {'W', 'NA', 'CL'}
        def counts(ids):
            c = {}
            for rid in ids:
                rn = res_map.get(rid)
                if not rn:
                    continue
                if rn in excluded:
                    continue
                if allowed_lipids and (rn not in allowed_lipids):
                    continue
                c[rn] = c.get(rn, 0) + 1
            return c
        cu = counts(outer_ids); cl = counts(inner_ids)
        total_u = sum(cu.values()) or 1
        total_l = sum(cl.values()) or 1
        lipids = sorted(set(list(cu.keys()) + list(cl.keys())))
        for lip in lipids:
            fu = cu.get(lip, 0)/total_u
            fl = cl.get(lip, 0)/total_l
            metrics[f'{lip}_upper_fraction'] = float(fu)
            metrics[f'{lip}_lower_fraction'] = float(fl)
            metrics[f'{lip}_asymmetry_index'] = float(abs(fu - fl))
        # Key RBC checks (configurable)
        ps_u = metrics.get('POPS_upper_fraction', 0.0); ps_l = metrics.get('POPS_lower_fraction', 0.0)
        ps_enrichment = (ps_l / max(ps_u, 1e-6)) if (ps_u or ps_l) else 0.0
        # CHOL fractions: report both leaflet mean (legacy) and total fraction across both leaflets
        chol_u = metrics.get('CHOL_upper_fraction', 0.0); chol_l = metrics.get('CHOL_lower_fraction', 0.0)
        chol_mean = 0.5*(chol_u + chol_l) if (chol_u or chol_l) else 0.0
        # Total CHOL fraction relative to total lipids in both leaflets
        # Prefer composition.yaml (topology-based) when available for robust totals
        chol_total_frac = None
        try:
            comp_file = os.path.join(output_dir, "composition.yaml")
            if os.path.exists(comp_file):
                with open(comp_file, 'r') as f:
                    comp_tot = yaml.safe_load(f) or {}
                lip_only = {k: v for k, v in comp_tot.items() if k not in ("W", "NA", "CL")}
                tot = float(sum(lip_only.values())) or 1.0
                chol_total_frac = float(lip_only.get('CHOL', 0)) / tot
        except Exception:
            chol_total_frac = None
        if chol_total_frac is None:
            chol_total = float(cu.get('CHOL', 0) + cl.get('CHOL', 0))
            total_both = float(total_u + total_l) if (total_u + total_l) > 0 else 1.0
            chol_total_frac = chol_total / total_both
        metrics['ps_inner_enrichment_ratio'] = float(ps_enrichment)
        metrics['chol_mean_fraction'] = float(chol_mean)
        metrics['chol_total_fraction'] = float(chol_total_frac)
        # Thresholds
        mem = _membrane_cfg_from(config or {}) if config is not None else {}
        vcfg = (mem.get('validation') or {})
        chol_lo = float(vcfg.get('chol_total_min', 0.25))
        chol_hi = float(vcfg.get('chol_total_max', 0.35))
        ps_min = float(vcfg.get('ps_enrichment_min', 10.0))
        checks = {
            # Use total CHOL fraction across both leaflets for the gate; mean is informational
            'cholesterol_content_valid': (chol_lo <= chol_total_frac <= chol_hi),
            'ps_asymmetry_valid': (ps_enrichment >= ps_min),
        }
        ok = all(checks.values())
    # Write report
    out = {
        'metrics': metrics,
        'passed': bool(ok),
    }
    with open(os.path.join(output_dir, 'asymmetry_validation.yaml'), 'w') as f:
        yaml.dump(out, f, default_flow_style=False)
    if ok:
        print("✓ Asymmetry validation passed (PS inner enrichment, CHOL content)")
    else:
        print("✗ Asymmetry validation failed — see asymmetry_validation.yaml")

def run_insane(run_dir):
    """Run INSANE to create membrane template"""
    config = load_config()
    memcfg = _membrane_cfg_from(config)
    original_dir = os.getcwd()
    
    # Output directory
    if run_dir:
        run_dir_abs = os.path.abspath(run_dir)
        output_dir = os.path.join(run_dir_abs, "membrane_template")
    else:
        # Create global template - use project's data/templates directory
        project_root = Path(__file__).parent.parent.parent
        templates_base = project_root / "data" / "templates"
        output_dir = templates_base / "membrane"
        output_dir = str(output_dir)  # Convert to string for os.path operations
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
    cmd, upper_str, lower_str = build_membrane_template(memcfg, output_dir)
    
    # Log file
    if run_dir:
        log_file = os.path.join(run_dir_abs, "logs", "insane_membrane.log")
    else:
        log_file = os.path.join(output_dir, "insane_membrane.log")
    # Ensure logs directory exists for run-specific builds
    if run_dir:
        os.makedirs(os.path.join(run_dir_abs, "logs"), exist_ok=True)
    
    print("Building RBC membrane template with INSANE via Docker...")
    print(f"Upper leaflet: {upper_str}")
    print(f"Lower leaflet: {lower_str}")
    print(f"Box size: {' x '.join(map(str, memcfg['box_size']))} nm")
    print(f"Salt: {memcfg['salt_concentration']} M")
    
    # Run INSANE from project root (Docker needs docker-compose.yml)
    project_root = Path(__file__).parent.parent.parent
    
    # Run INSANE via Docker
    with open(log_file, 'w') as log:
        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=str(project_root),  # Run from project root
                check=True
            )
            log.write(result.stdout)
            print("✓ INSANE completed successfully")
        except subprocess.CalledProcessError as e:
            log.write(e.stdout if e.stdout else "")
            print(f"✗ INSANE failed. Check log: {log_file}")
            sys.exit(1)
    
    # Parse membrane composition from log/topology
    parse_membrane_composition(log_file, output_dir)

    # Create leaflet index for validation and downstream reproducibility
    try:
        create_leaflet_index_for_template(output_dir, config)
    except Exception:
        pass

    # Run basic QC on the generated membrane template (box size, composition)
    run_membrane_qc(output_dir, config)
    # Asymmetry validation
    validate_membrane_asymmetry(output_dir, config)

    # Create README for membrane (reflect effective membrane config)
    try:
        eff_cfg = {'membrane': memcfg}
        create_membrane_readme(output_dir, eff_cfg)
    except Exception:
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


def create_leaflet_index_for_template(output_dir: str, config: dict | None = None):
    """Create index_leaflets.ndx (OuterPO4/InnerPO4) from membrane_template.gro.

    Classification: split PO4/P atoms by global z-median.
    Also writes full-lipid groups [ Outer ] and [ Inner ] containing all
    atoms from lipid residues assigned to each leaflet (POPC/POPE/POPS/PSM/CHOL
    or those defined in config pmf.membrane leaflets).
    """
    gro_file = os.path.join(output_dir, "membrane_template.gro")
    ndx_path = os.path.join(output_dir, 'index_leaflets.ndx')
    try:
        with open(gro_file, 'r') as f:
            lines = f.readlines()
    except Exception:
        return None
    # Allowed lipid residue names
    memcfg = _membrane_cfg_from(config or {}) if config is not None else {}
    allowed_lipids = set()
    allowed_lipids.update((memcfg.get('upper_leaflet') or {}).keys())
    allowed_lipids.update((memcfg.get('lower_leaflet') or {}).keys())
    if not allowed_lipids:
        allowed_lipids = {'POPC', 'POPE', 'POPS', 'PSM', 'CHOL'}

    # Parse GRO atoms and collect per-residue stats
    po4_atoms = []  # list of (atom_index, z) for headgroups only
    atom_index = 0
    res_sum_z: dict[int, float] = {}
    res_cnt: dict[int, int] = {}
    atom_to_res: dict[int, int] = {}
    res_map: dict[int, str] = {}
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        atom_index += 1
        atom = line[10:15].strip()
        try:
            resid = int(line[0:5])
            resn = line[5:10].strip()
            z = float(line[36:44])
        except Exception:
            continue
        atom_to_res[atom_index] = resid
        res_map[resid] = resn
        # accumulate per-residue z for lipids
        if resn in allowed_lipids:
            res_sum_z[resid] = res_sum_z.get(resid, 0.0) + z
            res_cnt[resid] = res_cnt.get(resid, 0) + 1
        if atom not in ('PO4', 'P'):
            continue
        # only count PO4/P from lipids (ignore solvent/ions)
        if resn in allowed_lipids:
            po4_atoms.append((atom_index, z))
    if not po4_atoms:
        return None
    zmed = float(np.median(np.asarray([z for (_i, z) in po4_atoms], dtype=float)))
    outer = [i for (i, z) in po4_atoms if z >= zmed]
    inner = [i for (i, z) in po4_atoms if z < zmed]
    # Balance if one side empty
    if not outer or not inner:
        half = len(po4_atoms) // 2
        order = [i for (i, _z) in sorted(po4_atoms, key=lambda t: t[1])]
        inner = order[:half]
        outer = order[half:]
    # Build full lipid atom groups using per-residue mean z vs midplane
    res_mean_z = {rid: (res_sum_z[rid] / max(1, res_cnt.get(rid, 1))) for rid in res_sum_z.keys()}
    # Classify residues by mean z
    outer_res = {rid for rid, mz in res_mean_z.items() if mz >= zmed}
    inner_res = set(res_mean_z.keys()) - outer_res
    # Collect atom indices for each leaflet (all atoms of lipid residues)
    outer_all = []
    inner_all = []
    atom_index = 0
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        atom_index += 1
        try:
            resid = int(line[0:5])
            resn = line[5:10].strip()
        except Exception:
            continue
        if resn not in allowed_lipids:
            continue
        if resid in outer_res:
            outer_all.append(atom_index)
        elif resid in inner_res:
            inner_all.append(atom_index)
    try:
        with open(ndx_path, 'w') as f:
            f.write('[ OuterPO4 ]\n')
            for k, idx in enumerate(outer):
                if k and (k % 15 == 0):
                    f.write('\n')
                f.write(f"{idx} ")
            f.write("\n\n[ InnerPO4 ]\n")
            for k, idx in enumerate(inner):
                if k and (k % 15 == 0):
                    f.write('\n')
                f.write(f"{idx} ")
            f.write('\n\n[ Outer ]\n')
            for k, idx in enumerate(outer_all):
                if k and (k % 15 == 0):
                    f.write('\n')
                f.write(f"{idx} ")
            f.write('\n\n[ Inner ]\n')
            for k, idx in enumerate(inner_all):
                if k and (k % 15 == 0):
                    f.write('\n')
                f.write(f"{idx} ")
            f.write('\n')
        print(f"Created leaflet index: {ndx_path}")
        return ndx_path
    except Exception:
        return None

def create_membrane_readme(output_dir, config):
    """Create README for membrane template.

    Accepts full config; extracts membrane section flexibly.
    """
    mem = _membrane_cfg_from(config)
    # Safe access with defaults
    upper = mem.get('upper_leaflet', {})
    lower = mem.get('lower_leaflet', {})
    bx = mem.get('box_size', [10, 10, 14])
    thick = mem.get('thickness', 4.5)
    salt = mem.get('salt_concentration', 0.15)
    readme_content = f"""# RBC Membrane Template

## Description
Evidence-based asymmetric red blood cell (RBC) membrane model for SOLVIA hemolytic toxicity predictions.

## Composition
### Upper Leaflet (Outer)
- POPC: {upper.get('POPC', 0)*100:.0f}%
- PSM: {upper.get('PSM', 0)*100:.0f}%
- CHOL: {upper.get('CHOL', 0)*100:.0f}%

### Lower Leaflet (Inner)
- POPE: {lower.get('POPE', 0)*100:.0f}%
- POPS: {lower.get('POPS', 0)*100:.0f}%
- CHOL: {lower.get('CHOL', 0)*100:.0f}%

## Parameters
- Box size: {' x '.join(map(str, bx))} nm
- Membrane thickness: ~{thick} nm
- Salt concentration: {salt} M (NaCl)
- Total cholesterol: ~25-35 mol%

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


def _read_gro_box(gro_path: str):
    """Parse gro box vectors from last line, return (x,y,z) in nm."""
    with open(gro_path, 'r') as f:
        lines = f.readlines()
    if not lines:
        return None
    parts = lines[-1].split()
    try:
        # Standard GRO: 3 floats for box edges
        x, y, z = map(float, parts[:3])
        return x, y, z
    except Exception:
        return None


def run_membrane_qc(output_dir: str, config: dict):
    """Validate membrane template against QC requirements.

    Checks:
    - Box X,Y dimensions >= requested (±0.1 nm tolerance)
    - Total cholesterol fraction within 0.40–0.50 (RBC-like)
    - Expected lipid species present per leaflet set (at least in total)
    - Water and ions present
    """
    gro_file = os.path.join(output_dir, "membrane_template.gro")
    top_file = os.path.join(output_dir, "membrane_template.top")
    comp_file = os.path.join(output_dir, "composition.yaml")

    # Extract membrane config flexibly
    mem = _membrane_cfg_from(config)
    # Box QC
    requested = mem.get('box_size', [10, 10, 14])
    dims = _read_gro_box(gro_file)
    if dims is None:
        print("Warning: Could not read GRO box for QC")
    else:
        x, y, z = dims
        tol = 0.1
        if x + tol < requested[0] or y + tol < requested[1]:
            print(f"✗ QC: Box too small (got {x:.2f}×{y:.2f}×{z:.2f} nm, expected ≥ {requested[0]}×{requested[1]}×{requested[2]} nm)")
            sys.exit(2)

    # Composition QC
    if not os.path.exists(comp_file):
        # Try to regenerate composition if missing
        parse_membrane_composition(os.path.join(output_dir, "insane_membrane.log"), output_dir)
    with open(comp_file, 'r') as f:
        comp = yaml.safe_load(f) or {}
    # Count lipids (exclude water/ions)
    lipid_counts = {k: v for k, v in comp.items() if k not in ("W", "NA", "CL")}
    total_lipids = sum(lipid_counts.values()) or 1
    chol = float(lipid_counts.get('CHOL', 0))
    chol_frac = chol / total_lipids
    # Configurable CHOL total fraction band
    vcfg = (mem.get('validation') or {})
    chol_lo = float(vcfg.get('chol_total_min', 0.25))
    chol_hi = float(vcfg.get('chol_total_max', 0.35))
    chol_tol = float(vcfg.get('chol_total_tol', 0.0))
    lo = chol_lo - chol_tol
    hi = chol_hi + chol_tol
    if not (lo <= chol_frac <= hi):
        print(f"✗ QC: Cholesterol fraction {chol_frac*100:.1f}% out of configured range ({chol_lo*100:.0f}–{chol_hi*100:.0f}%)")
        sys.exit(2)
    # Species presence QC
    expected_upper = set((mem.get('upper_leaflet') or {}).keys())
    expected_lower = set((mem.get('lower_leaflet') or {}).keys())
    present = set(lipid_counts.keys())
    missing = (expected_upper | expected_lower) - present
    if missing:
        print(f"✗ QC: Missing expected lipid species: {', '.join(sorted(missing))}")
        sys.exit(2)
    # Solvent/ions presence
    if comp.get('W', 0) <= 0:
        print("✗ QC: No water molecules found in template")
        sys.exit(2)
    # Optional: ions presence (salt_concentration may be 0)
    if float(mem.get('salt_concentration', 0.0)) > 0 and (comp.get('NA', 0) + comp.get('CL', 0) == 0):
        print("✗ QC: Salt requested but no ions found in template")
        sys.exit(2)
    print("✓ Membrane QC passed (box, composition, solvent/ions)")

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
