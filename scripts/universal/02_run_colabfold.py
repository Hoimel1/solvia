#!/usr/bin/env python3
"""
SOLVIA ColabFold Structure Prediction
Runs ColabFold for peptide structure prediction or selects best existing model
"""

import os
import sys
import yaml
import json
import argparse
import subprocess
from pathlib import Path
import numpy as np
import shutil

def read_fasta(path):
    with open(path, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    seq = ''.join(l for l in lines if not l.startswith('>'))
    return seq

def _kyte_doolittle_scores():
    # Kyte-Doolittle hydropathy index
    return {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

def predict_tm_segments(seq: str, window: int = 19, threshold: float = 1.6):
    """Predict TM helices using Kyte-Doolittle sliding window.
    Returns list of (start, end) indices (0-based, inclusive) for putative segments.
    Optimized and with corrected segment ends.
    """
    seq = (seq or '').upper()
    if not seq or len(seq) < window:
        return []
    kd = _kyte_doolittle_scores()
    import numpy as _np
    vals = _np.array([kd.get(ch, 0.0) for ch in seq], dtype=float)
    n = len(vals)
    cs = _np.cumsum(_np.r_[0.0, vals])
    ma = (cs[window:] - cs[:-window]) / float(window)  # length n-window+1
    segs = []
    i = 0
    last = n - window
    while i <= last:
        if ma[i] >= threshold:
            start = i
            j = i + 1
            while j <= last and ma[j] >= threshold:
                j += 1
            end = j + window - 2  # inclusive end index
            segs.append((start, min(end, n - 1)))
            i = j
        else:
            i += 1
    merged = []
    for s, e in segs:
        if not merged or s > merged[-1][1] + 1:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(int(a), int(b)) for a, b in merged]

def analyze_sequence_properties(seq: str):
    seq = (seq or '').upper()
    n = len(seq)
    if n == 0:
        return {'length': 0, 'cysteine_fraction': 0.0, 'has_transmembrane': False, 'tm_segments': []}
    cyst_frac = seq.count('C') / n
    tm_segs = predict_tm_segments(seq)
    return {'length': n, 'cysteine_fraction': cyst_frac, 'has_transmembrane': bool(tm_segs), 'tm_segments': tm_segs}

def _tool_available(name: str) -> bool:
    try:
        return shutil.which(name) is not None
    except Exception:
        return False

def predict_tm_segments_tools(fasta_path: str):
    """Try TMHMM or Phobius if available; return list of (start,end) 0-based inclusive.
    Fallback to empty list if tools unavailable or parsing fails.
    """
    def _parse_segments_from_text(txt: str):
        import re
        segs = []
        for line in (txt or '').splitlines():
            u = line.strip()
            if not u:
                continue
            U = u.upper()
            # 1) TMHMM Topology lines: only 'M' (membrane) segments
            m = re.search(r'TOPOLOGY\s*=\s*(.*)', U)
            if m:
                topo = m.group(1)
                # Find label-range triples like M12-34, I1-10, O11-20 (case-insensitive)
                for lab, a, b in re.findall(r'([IOM])\s*(\d+)\s*(?:[-.]{1,2}|\s+)\s*(\d+)', topo, flags=re.IGNORECASE):
                    if lab.upper() != 'M':
                        continue
                    try:
                        s = max(0, int(a) - 1)
                        e = max(s, int(b) - 1)
                        segs.append((s, e))
                    except Exception:
                        pass
                continue
            # 2) Phobius/TMHMM long output: only explicit TM lines
            for a, b in re.findall(r'(?:TRANSMEM|TMHELIX)\s+(\d+)\s*(?:[-.]{1,2}|\s+)\s*(\d+)', U):
                try:
                    s = max(0, int(a) - 1)
                    e = max(s, int(b) - 1)
                    segs.append((s, e))
                except Exception:
                    pass
        # Merge overlapping/adjacent
        segs.sort()
        merged = []
        for s, e in segs:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        return [(int(a), int(b)) for a, b in merged]

    # Try TMHMM (portable: pipe FASTA on stdin with --short)
    try:
        if _tool_available('tmhmm'):
            with open(fasta_path, 'r') as fh:
                fa = fh.read()
            res = subprocess.run(['tmhmm', '--short'], input=fa, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
            segs = _parse_segments_from_text(res.stdout or '')
            if segs:
                return segs
    except Exception:
        pass
    # Try Phobius (binary name can be 'phobius' or 'phobius.pl')
    try:
        for exe in ('phobius', 'phobius.pl'):
            if _tool_available(exe):
                # Use -long if supported to include TRANSMEM lines
                res = subprocess.run([exe, '-long', fasta_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
                segs = _parse_segments_from_text(res.stdout or '')
                if segs:
                    return segs
                # Fallback: pipe FASTA via stdin
                with open(fasta_path, 'r') as fh:
                    fa = fh.read()
                res = subprocess.run([exe, '-long'], input=fa, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
                segs = _parse_segments_from_text(res.stdout or '')
                if segs:
                    return segs
    except Exception:
        pass
    return []

def analyze_sequence_properties_fasta(fasta_path: str):
    seq = read_fasta(fasta_path)
    props = analyze_sequence_properties(seq)
    try:
        tool_segs = predict_tm_segments_tools(fasta_path)
        if tool_segs:
            props['tm_segments'] = tool_segs
            props['has_transmembrane'] = True
    except Exception:
        pass
    return props

def optimize_colabfold_params(base_cfg: dict, fasta_path: str):
    cfg = dict(base_cfg)
    try:
        props = analyze_sequence_properties_fasta(fasta_path)
        n = props['length']
        if n and n < 50:
            cfg['num_models'] = max(cfg.get('num_models', 5), 10)
            cfg['num_replicates'] = max(cfg.get('num_replicates', 5), 10)
            cfg['msa_mode'] = 'single_sequence'
        if props['cysteine_fraction'] > 0.1:
            cfg['relax'] = True
        if props['has_transmembrane']:
            cfg['num_replicates'] = max(cfg.get('num_replicates', 5), 8)
    except Exception:
        pass
    return cfg

def parse_pdb_plddt(pdb_file):
    """Extract per-residue pLDDT from CA B-factors and CA coordinates."""
    vals = []
    ca_coords = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                if line[12:16].strip() != 'CA':
                    continue
                try:
                    b = float(line[60:66])
                    vals.append(b)
                except Exception:
                    pass
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    ca_coords.append([x, y, z])
                except Exception:
                    pass
    except Exception:
        return None, None
    return (np.array(vals) if vals else None), (np.array(ca_coords) if ca_coords else None)

def kabsch_rmsd(P, Q):
    """Compute RMSD after optimal superposition (Kabsch)."""
    if P is None or Q is None:
        return None
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape or P.shape[0] < 3:
        return None
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(C)  # C = U @ S @ Vt
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0.0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    Prot = Pc @ R
    diff = Prot - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))

def clash_score(pdb_file, softness=0.4):
    """Atom-typed clash score using vdW radii; returns clashes per atom.
    Coordinates are in Å from PDB. softness (Å) tolerates small overlaps.
    """
    coords = []
    res_ids = []
    elems = []
    with open(pdb_file,'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            try:
                resi = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                name = line[12:16].strip()
                elem = line[76:78].strip() if len(line) >= 78 and line[76:78].strip() else name[0]
            except Exception:
                continue
            coords.append((x,y,z)); res_ids.append(resi); elems.append(elem.upper())
    N = len(coords)
    if N < 2:
        return 0.0
    coords = np.array(coords)
    vdw = {'H':1.20,'C':1.70,'N':1.55,'O':1.52,'S':1.80,'P':1.80}
    clashes = 0
    for i in range(N-1):
        for j in range(i+1, N):
            if abs(res_ids[i]-res_ids[j]) <= 2:
                continue
            d = np.linalg.norm(coords[i]-coords[j])
            ri = vdw.get(elems[i], 1.7); rj = vdw.get(elems[j], 1.7)
            thr = max(0.0, ri + rj - softness)
            if d < thr:
                clashes += 1
    return clashes / N

def _dihedral(a, b, c, d):
    # a,b,c,d: 3D vectors
    import numpy as _np
    b0 = a - b
    b1 = c - b
    b2 = d - c
    b1 /= _np.linalg.norm(b1) if _np.linalg.norm(b1) != 0 else 1.0
    v = b0 - _np.dot(b0, b1) * b1
    w = b2 - _np.dot(b2, b1) * b1
    x = _np.dot(v, w)
    y = _np.dot(_np.cross(b1, v), w)
    return float(_np.degrees(_np.arctan2(y, x)))

def _parse_backbone_coords(pdb_file):
    # returns dict resid -> {'N':xyz,'CA':xyz,'C':xyz,'SG':xyz}
    bb = {}
    with open(pdb_file,'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            try:
                atom = line[12:16].strip()
                res = line[17:20].strip()
                resid = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            if resid not in bb:
                bb[resid] = {'res': res}
            if atom in ('N','CA','C','SG'):
                bb[resid][atom] = np.array([x,y,z], dtype=float)
    return bb

def ramachandran_quality(pdb_file):
    bb = _parse_backbone_coords(pdb_file)
    resids = sorted(bb.keys())
    favored = 0; allowed = 0; total = 0
    helix_like = 0; beta_like = 0
    for i in range(1, len(resids)-1):
        r_prev = bb.get(resids[i-1], {})
        r = bb.get(resids[i], {})
        r_next = bb.get(resids[i+1], {})
        if not all(k in r for k in ('N','CA','C')):
            continue
        if 'C' not in r_prev or 'N' not in r_next:
            continue
        try:
            phi = _dihedral(r_prev['C'], r['N'], r['CA'], r['C'])
            psi = _dihedral(r['N'], r['CA'], r['C'], r_next['N'])
        except Exception:
            continue
        total += 1
        resname = (r.get('res') or '').upper()
        # residue-class specific boxes (simplified MolProbity-like regions)
        # General (non-Gly/Pro)
        gen_fav = ((-160 <= phi <= -40 and -80 <= psi <= 50) or (-180 <= phi <= -70 and 90 <= psi <= 180))
        gen_all = gen_fav or ((-180 <= phi <= -20 and -150 <= psi <= 90) or (-180 <= phi <= -20 and 60 <= psi <= 180))
        # Glycine (wider) — allowed but not entire plane
        gly_fav = ((-180 <= phi <= -20 and -100 <= psi <= 100) or (-100 <= phi <= 100 and 90 <= psi <= 180))
        gly_all = (
            (-180 <= phi <=  -20 and -180 <= psi <= 160) or
            ( -100 <= phi <=  100 and  -80 <= psi <=  80) or
            (  100 <= phi <=  180 and   90 <= psi <= 180)
        )
        # Proline (restricted phi)
        pro_fav = (-90 <= phi <= -30 and -80 <= psi <= 50)
        pro_all = pro_fav or (-120 <= phi <= -10 and -120 <= psi <= 90)
        if resname == 'GLY':
            fav = gly_fav; allow = gly_all
        elif resname == 'PRO':
            fav = pro_fav; allow = pro_all
        else:
            fav = gen_fav; allow = gen_all
        if fav:
            favored += 1
        elif allow:
            allowed += 1
        # simple ss proxies
        if (-90 <= phi <= -30 and -80 <= psi <= -20):
            helix_like += 1
        if (-180 <= phi <= -100 and 90 <= psi <= 180):
            beta_like += 1
    frac_favored = (favored/total) if total else 0.0
    frac_allowed = ((favored+allowed)/total) if total else 0.0
    ss = {
        'helix_frac': (helix_like/total) if total else 0.0,
        'beta_frac': (beta_like/total) if total else 0.0
    }
    return {'favored': frac_favored, 'allowed': frac_allowed, 'n_eval': total, 'secondary_structure': ss}

def disulfide_analysis(pdb_file):
    # detect SG-SG within 1.9–2.2 Å; bad if <1.9 or >2.2 (pairs listed up to 2.5 Å)
    bb = _parse_backbone_coords(pdb_file)
    cys_res = [rid for rid, rec in bb.items() if rec.get('res','').upper().startswith('CYS') and 'SG' in rec]
    pairs = []
    bad = 0
    for i in range(len(cys_res)-1):
        for j in range(i+1, len(cys_res)):
            a = bb[cys_res[i]]['SG']; b = bb[cys_res[j]]['SG']
            d = float(np.linalg.norm(a-b))
            if d <= 2.5:
                pairs.append({'res1': int(cys_res[i]), 'res2': int(cys_res[j]), 'distance_A': float(d)})
                if d < 1.9 or d > 2.2:
                    bad += 1
    return {'count': len(pairs), 'bad_count': bad, 'pairs': pairs}

def load_config():
    """Load SOLVIA configuration"""
    config_dir = Path(__file__).parent.parent.parent / "config"
    
    # Try available config files in order of preference
    config_files = [
        "pmf_standard_config.yaml",
        "config.yaml"
    ]
    
    for config_file in config_files:
        config_path = config_dir / config_file
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                
            # Add default ColabFold values if missing
            if 'colabfold' not in config:
                config['colabfold'] = {}
            
            colabfold_defaults = {
                'num_replicates': 5,
                'num_models': 5,
                'msa_mode': 'mmseqs2_uniref_env',
                'pair_mode': 'unpaired_paired',
                'relax': True,
                'min_plddt': 70,
                'selection': {
                    'method': 'weighted',  # weighted | rank
                    'weights': {
                        'plddt_mean': 0.30,
                        'plddt90_frac': 0.20,
                        'plddt70_frac': 0.15,
                        'plddt_lt50_frac': -0.10,
                        'ptm': 0.10,
                        'ensemble_agreement': 0.10,
                        'rama_favored': 0.05,
                        'disulfide_ok': 0.05,
                        'clash_penalty': -0.10
                    }
                }
            }
            
            for key, default_value in colabfold_defaults.items():
                if key not in config['colabfold']:
                    config['colabfold'][key] = default_value
                
            print(f"Loaded configuration from: {config_path}")
            return config
    
    # If no config file found, create a minimal default config
    print("Warning: No config file found, using defaults")
    return {
        'colabfold': {
            'num_replicates': 5,
            'num_models': 5,
            'msa_mode': 'mmseqs2_uniref_env',
            'pair_mode': 'unpaired_paired',
            'relax': True,
            'min_plddt': 70
        }
    }

def load_run_metadata(run_dir):
    """Load run metadata"""
    metadata_path = os.path.join(run_dir, "metadata.yaml")
    with open(metadata_path, 'r') as f:
        return yaml.safe_load(f)

def check_colabfold_available():
    """Check if colabfold_batch is available in PATH (portable)."""
    try:
        return shutil.which("colabfold_batch") is not None
    except Exception:
        return False

def run_colabfold(run_dir):
    """Run ColabFold structure prediction or select best model if already run"""
    config = load_config()
    metadata = load_run_metadata(run_dir)
    
    # Paths
    input_fasta = os.path.join(run_dir, "input", "peptide.fasta")
    output_dir = os.path.join(run_dir, "colabfold")
    
    # Check if ColabFold outputs already exist
    pdb_files = list(Path(output_dir).glob("*_unrelaxed_rank_*.pdb")) if os.path.exists(output_dir) else []
    relaxed_files = list(Path(output_dir).glob("*_relaxed_rank_*.pdb")) if os.path.exists(output_dir) else []
    
    if pdb_files or relaxed_files:
        print(f"Found {len(pdb_files + relaxed_files)} ColabFold models. Selecting best...")
        return select_best_model(output_dir, config)
    
    # Check if we can run ColabFold natively
    colabfold_available = check_colabfold_available()
    
    if not colabfold_available:
        print("")
        print("="*70)
        print("ColabFold not found in PATH.")
        print("="*70)
        print("")
        print("Please run ColabFold using Docker with this command:")
        print("")
        print(f"docker compose run --rm \\")
        print(f"  -e XDG_CACHE_HOME=/cache \\")
        print(f"  -e MPLCONFIGDIR=/cache/mpl \\")
        print(f"  colabfold \\")
        print(f"  /work/{run_dir}/input/peptide.fasta /work/{run_dir}/colabfold \\")
        print(f"  --num-seeds {config['colabfold']['num_replicates']} \\")
        print(f"  --num-models {config['colabfold']['num_models']} \\")
        print(f"  --msa-mode {config['colabfold']['msa_mode']} \\")
        print(f"  --pair-mode {config['colabfold']['pair_mode']}")
        print("")
        print("After ColabFold completes, run this script again to select the best model.")
        print("="*70)
        return None, 0
    
    # Optimize parameters based on sequence properties
    cfg_cf = optimize_colabfold_params(config['colabfold'], input_fasta)

    # Build ColabFold command (native installation)
    cmd = [
        "colabfold_batch",
        input_fasta,
        output_dir,
        "--num-seeds", str(cfg_cf['num_replicates']),
        "--num-models", str(cfg_cf['num_models']),
        "--msa-mode", cfg_cf['msa_mode'],
        "--pair-mode", cfg_cf['pair_mode'],
    ]

    if cfg_cf.get('relax', True):
        cmd.append("--amber")
        cmd.append("--use-gpu-relax")
    
    # Log file
    log_dir = os.path.join(run_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "colabfold.log")
    
    print(f"Running ColabFold for {metadata['peptide_id']}...")
    print(f"Command: {' '.join(cmd)}")
    print(f"This may take a while depending on sequence length and settings...")
    
    # Run ColabFold
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
            print("✓ ColabFold completed successfully")
        except subprocess.CalledProcessError as e:
            log.write(e.stdout if e.stdout else "")
            print(f"✗ ColabFold failed. Check log: {log_file}")
            sys.exit(1)
    
    # Select best model
    return select_best_model(output_dir, config)

def select_best_model_by_bfactor(pdb_files):
    """Select best model based on B-factors when score files are missing"""
    best_pdb = None
    best_plddt = 0
    
    for pdb_file in pdb_files:
        # Read PDB and extract B-factors (pLDDT scores)
        plddt_values = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        b_factor = float(line[60:66])
                        plddt_values.append(b_factor)
                    except:
                        pass
        
        if plddt_values:
            avg_plddt = np.mean(plddt_values)
            print(f"  {os.path.basename(pdb_file)}: pLDDT = {avg_plddt:.1f}")
            if avg_plddt > best_plddt:
                best_plddt = avg_plddt
                best_pdb = str(pdb_file)
    
    return best_pdb, best_plddt

def select_best_model(output_dir, config):
    """Select best model using multi-criteria scientific scoring and export ensemble analysis."""
    # Find all PDB files first
    import glob
    pdb_files = glob.glob(os.path.join(output_dir, "*_unrelaxed_rank_*.pdb"))
    relaxed_pdb_files = glob.glob(os.path.join(output_dir, "*_relaxed_rank_*.pdb"))
    all_pdb_files = (relaxed_pdb_files if relaxed_pdb_files else pdb_files)
    
    if not all_pdb_files:
        print("Error: No PDB files found in", output_dir)
        print("Please run ColabFold first using the Docker command shown above.")
        return None, 0
    
    # Find all score files
    score_files = glob.glob(os.path.join(output_dir, "*_scores_rank_*.json"))
    # Also consider ranking_debug.json produced by (Colab)Fold
    ranking_debug = os.path.join(output_dir, "ranking_debug.json")
    
    # Gather per-model metrics
    models = []
    # Map rank -> score data (unused map removed)
    if not score_files:
        print("Warning: No score files found, extracting pLDDT from B-factors...")
        for pdb_file in all_pdb_files:
            plddt_vec, ca = parse_pdb_plddt(pdb_file)
            plddt_mean = float(np.mean(plddt_vec)) if plddt_vec is not None and plddt_vec.size else 0.0
            models.append({'pdb': pdb_file, 'rank': None, 'plddt_vec': plddt_vec, 'plddt_mean': plddt_mean, 'ptm': None, 'ca': ca})
    else:
        # Extract pLDDT list and ptm if present
        plddt_scores = {}
        for score_file in score_files:
            with open(score_file, 'r') as f:
                data = json.load(f)
            
            # Extract rank from filename
            import re
            match = re.search(r'rank_(\d+)', os.path.basename(score_file))
            if match:
                rank = int(match.group(1))
                model_name = f"rank_{rank:03d}"
                
                # pLDDT vector if available
                plddt_vec = None
                if 'plddt' in data and isinstance(data['plddt'], list):
                    plddt_vec = np.array(data['plddt'], dtype=float)
                elif 'plddts' in data and isinstance(data['plddts'], list):
                    plddt_vec = np.array(data['plddts'], dtype=float)
                plddt_mean = float(np.mean(plddt_vec)) if plddt_vec is not None and plddt_vec.size else (float(data.get('plddt', data.get('plddts', 0.0))))
                ptm = data.get('ptm', data.get('iptm', None))
                plddt_scores[model_name] = plddt_mean
                # map rank to pdb and metrics
                rank_num = rank
                pdb_match = None
                for pdb_file in all_pdb_files:
                    if f"rank_{rank_num:03d}" in pdb_file or f"rank_{rank_num:01d}" in pdb_file:
                        pdb_match = pdb_file; break
                plddt_vec_pdb, ca = parse_pdb_plddt(pdb_match) if pdb_match else (None, None)
                # Prefer JSON pLDDT vector when available; else fall back to parsed PDB B-factors
                plddt_vec_final = plddt_vec if plddt_vec is not None else plddt_vec_pdb
                models.append({'pdb': pdb_match, 'rank': rank_num, 'plddt_vec': plddt_vec_final, 'plddt_mean': plddt_mean, 'ptm': ptm, 'ca': ca})

        # Optional: augment missing means from ranking_debug.json
        try:
            if os.path.exists(ranking_debug):
                with open(ranking_debug, 'r') as f:
                    rd = json.load(f)
                pl_map = rd.get('plddts') or rd.get('plddt') or {}
                for i, m in enumerate(models):
                    if (not m.get('plddt_mean')) and m.get('pdb'):
                        base = os.path.basename(m['pdb'])
                        for k, v in pl_map.items():
                            if str(k) in base:
                                try:
                                    models[i]['plddt_mean'] = float(v)
                                except Exception:
                                    pass
                                break
        except Exception:
            pass

    # Compute ensemble RMSD matrix
    n_models = len(models)
    rmsd_matrix = [[None]*n_models for _ in range(n_models)]
    for i in range(n_models):
        for j in range(i+1, n_models):
            rms = kabsch_rmsd(models[i]['ca'], models[j]['ca'])
            rmsd_matrix[i][j] = rmsd_matrix[j][i] = (float(rms) if rms is not None else None)

    # Compute composite scores (configurable)
    sel_cfg = (config.get('colabfold') or {}).get('selection', {})
    method = str(sel_cfg.get('method', 'weighted')).lower()
    weights = (sel_cfg.get('weights') or {})
    composite = []
    ranks = []
    for i, m in enumerate(models):
        pvec = m['plddt_vec']
        mean_plddt = float(m['plddt_mean'] or 0.0)
        if pvec is None or pvec.size == 0:
            p90 = p70 = low = 0.0
        else:
            p90 = float(np.mean(pvec >= 90.0))
            p70 = float(np.mean(pvec >= 70.0))
            low = float(np.mean(pvec < 50.0))
        ptm = m['ptm'] if m['ptm'] is not None else 0.5
        # clash penalty (lower is better)
        cpen = clash_score(m['pdb']) if m['pdb'] else 0.0
        # ensemble agreement score
        if n_models >= 2:
            vals = [v for v in rmsd_matrix[i] if v is not None]
            mean_r = float(np.mean(vals)) if vals else None
            agree = (max(0.0, 1.0 - min(mean_r, 5.0)/5.0) if mean_r is not None else 0.5)
        else:
            agree = 0.5
        # validation-derived metrics
        rama_favored = 0.0
        disulf_ok = 1.0
        # place holders, adjusted below when validations computed
        metrics = {
            'plddt_mean': (mean_plddt / 100.0),
            'plddt90_frac': p90,
            'plddt70_frac': p70,
            'plddt_lt50_frac': low,
            'ptm': float(ptm) if isinstance(ptm, (int, float)) else 0.5,
            'ensemble_agreement': agree,
            'clash_penalty': min(1.0, cpen)  # clamp
        }
        composite.append(metrics)

    # Structural validation for each model (Ramachandran + disulfide)
    validations = []
    for m in models:
        try:
            rama = ramachandran_quality(m['pdb']) if m['pdb'] else {'favored': 0.0, 'allowed': 0.0, 'n_eval': 0}
            ssb = disulfide_analysis(m['pdb']) if m['pdb'] else {'count': 0, 'bad_count': 0}
        except Exception:
            rama = {'favored': 0.0, 'allowed': 0.0, 'n_eval': 0}
            ssb = {'count': 0, 'bad_count': 0}
        validations.append({'ramachandran': rama, 'disulfides': ssb})

    # Finalize composite scores per config (and collect selection metrics for logging)
    final_scores = []
    selection_metrics = []
    if method == 'rank':
        # Rank-aggregate each metric then sum ranks (lower is better)
        metric_names = ['plddt_mean','plddt90_frac','plddt70_frac','ptm','ensemble_agreement','rama_favored','disulfide_ok']
        arrays = {k: [] for k in metric_names}
        for i, metrics in enumerate(composite):
            arrays['plddt_mean'].append(metrics['plddt_mean'])
            arrays['plddt90_frac'].append(metrics['plddt90_frac'])
            arrays['plddt70_frac'].append(metrics['plddt70_frac'])
            arrays['ptm'].append(metrics['ptm'])
            arrays['ensemble_agreement'].append(metrics['ensemble_agreement'])
            arrays['rama_favored'].append(validations[i]['ramachandran']['favored'])
            arrays['disulfide_ok'].append(1.0 if validations[i]['disulfides']['bad_count'] == 0 else 0.0)
        # compute per-metric ranks (descending desirable)
        ranks_by_metric = {}
        ranksums = [0.0] * len(composite)
        for k, vals in arrays.items():
            order = np.argsort(-np.array(vals, dtype=float))
            rank = np.empty_like(order, dtype=float)
            rank[order] = np.arange(1, len(order) + 1)
            ranks_by_metric[k] = rank.tolist()
            for i in range(len(composite)):
                ranksums[i] += rank[i]
        final_scores = [-r for r in ranksums]
        # build selection metrics per model
        for i, metrics in enumerate(composite):
            mvals = dict(metrics)
            mvals['rama_favored'] = validations[i]['ramachandran']['favored']
            mvals['disulfide_ok'] = 1.0 if validations[i]['disulfides']['bad_count'] == 0 else 0.0
            selection_metrics.append({
                'model': os.path.basename(models[i]['pdb']) if models[i]['pdb'] else None,
                'method': 'rank',
                'metric_values': mvals,
                'metric_ranks': {k: ranks_by_metric[k][i] for k in ranks_by_metric},
                'rank_sum': float(-final_scores[i]),
                'final_score': float(final_scores[i])
            })
    else:
        # Weighted sum per config
        w = {
            'plddt_mean': float(weights.get('plddt_mean', 0.30)),
            'plddt90_frac': float(weights.get('plddt90_frac', 0.20)),
            'plddt70_frac': float(weights.get('plddt70_frac', 0.15)),
            'plddt_lt50_frac': float(weights.get('plddt_lt50_frac', -0.10)),
            'ptm': float(weights.get('ptm', 0.10)),
            'ensemble_agreement': float(weights.get('ensemble_agreement', 0.10)),
            'rama_favored': float(weights.get('rama_favored', 0.05)),
            'disulfide_ok': float(weights.get('disulfide_ok', 0.05)),
            'clash_penalty': float(weights.get('clash_penalty', -0.10)),
        }
        for i, metrics in enumerate(composite):
            score = (
                w['plddt_mean'] * metrics['plddt_mean'] +
                w['plddt90_frac'] * metrics['plddt90_frac'] +
                w['plddt70_frac'] * metrics['plddt70_frac'] +
                w['plddt_lt50_frac'] * metrics['plddt_lt50_frac'] +
                w['ptm'] * metrics['ptm'] +
                w['ensemble_agreement'] * metrics['ensemble_agreement'] +
                w['rama_favored'] * validations[i]['ramachandran']['favored'] +
                w['disulfide_ok'] * (1.0 if validations[i]['disulfides']['bad_count'] == 0 else 0.0) +
                w['clash_penalty'] * metrics['clash_penalty']
            )
            final_scores.append(score)
            mvals = dict(metrics)
            mvals['rama_favored'] = validations[i]['ramachandran']['favored']
            mvals['disulfide_ok'] = 1.0 if validations[i]['disulfides']['bad_count'] == 0 else 0.0
            selection_metrics.append({
                'model': os.path.basename(models[i]['pdb']) if models[i]['pdb'] else None,
                'method': 'weighted',
                'metric_values': mvals,
                'weights': w,
                'final_score': float(score)
            })

    best_idx = int(np.argmax(final_scores)) if final_scores else 0
    best_pdb = models[best_idx]['pdb'] if models else None
    best_plddt = float(models[best_idx]['plddt_mean']) if models else 0.0
    
    # Check minimum pLDDT threshold
    if best_plddt < config['colabfold']['min_plddt']:
        print(f"⚠ Warning: Best pLDDT ({best_plddt:.1f}) is below threshold ({config['colabfold']['min_plddt']})")
        print("Consider using a longer sequence or different modeling approach")
    
    # Save selection metadata
    selection = {
        'best_model': os.path.basename(best_pdb) if best_pdb else None,
        'best_plddt': float(best_plddt),
        'best_pdb': best_pdb,
        'min_plddt_threshold': config['colabfold']['min_plddt'],
        'models': [
            {
                'pdb': os.path.basename(m['pdb']) if m['pdb'] else None,
                'rank': m['rank'],
                'plddt_mean': float(m['plddt_mean'] or 0.0),
                'plddt90_frac': float(np.mean(m['plddt_vec']>=90.0)) if m['plddt_vec'] is not None and m['plddt_vec'].size else None,
                'plddt70_frac': float(np.mean(m['plddt_vec']>=70.0)) if m['plddt_vec'] is not None and m['plddt_vec'].size else None,
                'plddt_lt50_frac': float(np.mean(m['plddt_vec']<50.0)) if m['plddt_vec'] is not None and m['plddt_vec'].size else None,
                'ptm': (float(m['ptm']) if isinstance(m['ptm'], (int,float)) else None),
                'composite_score': float(final_scores[i]),
                'validation': validations[i]
            }
            for i, m in enumerate(models)
        ]
    }
    
    selection_file = os.path.join(output_dir, "model_selection.yaml")
    with open(selection_file, 'w') as f:
        yaml.dump(selection, f, default_flow_style=False)
    
    # Save ensemble RMSD
    try:
        ens = {
            'rmsd_matrix': rmsd_matrix,
            'mean_pairwise_rmsd': float(np.nanmean([v for row in rmsd_matrix for v in (row or []) if v is not None])) if rmsd_matrix else None
        }
        with open(os.path.join(output_dir, 'ensemble_analysis.yaml'), 'w') as f:
            yaml.dump(ens, f, default_flow_style=False)
    except Exception:
        pass

    # Save selection metrics for auditability
    try:
        selm = {
            'method': method,
            'entries': selection_metrics
        }
        with open(os.path.join(output_dir, 'selection_metrics.yaml'), 'w') as f:
            yaml.dump(selm, f, default_flow_style=False)
    except Exception:
        pass

    # Create symbolic link to best model
    if best_pdb and os.path.exists(best_pdb):
        best_link = os.path.join(output_dir, "best_model.pdb")
        if os.path.exists(best_link):
            os.remove(best_link)
        try:
            os.symlink(os.path.basename(best_pdb), best_link)
        except Exception:
            # Fallback: copy if symlink unsupported
            import shutil as _sh
            _sh.copy2(best_pdb, best_link)
    
    print(f"\n✓ Model selection complete")
    print(f"  Best model: {os.path.basename(best_pdb) if best_pdb else 'None'}")
    print(f"  pLDDT: {best_plddt:.1f}")
    print(f"  Selection saved: {selection_file}")
    
    return best_pdb, best_plddt

def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold structure prediction or select best model"
    )
    parser.add_argument(
        "run_dir",
        help="Run directory containing input FASTA"
    )
    
    args = parser.parse_args()
    
    # Check run directory
    if not os.path.exists(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}")
        sys.exit(1)
    
    # Run ColabFold or select model
    best_pdb, best_plddt = run_colabfold(args.run_dir)
    
    if best_pdb:
        print(f"\n✓ Ready for next step: Coarse-graining")
        print(f"  Command: python3 scripts/universal/03_coarse_grain.py {args.run_dir}")

if __name__ == "__main__":
    main()
