#!/usr/bin/env python3
"""
Generate index_leaflets.ndx by splitting PO4/P atoms across the membrane midplane.

Usage:
  python scripts/utils/generate_leaflets_index.py RUN_DIR [--gro path/to/structure.gro] [--out index_leaflets.ndx]

Defaults:
  - GRO defaults to RUN_DIR/equilibration/npt/npt.gro
  - OUT defaults to RUN_DIR/index_leaflets.ndx
"""
import argparse
from pathlib import Path
import sys
import numpy as np


def read_gro_atoms(gro_path: Path):
    atoms = []
    with open(gro_path, 'r') as f:
        lines = f.readlines()
    if len(lines) < 3:
        return atoms
    for line in lines[2:-1]:
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            x = float(line[20:28]); y = float(line[28:36]); z = float(line[36:44])
            atoms.append((resname, atomname, x, y, z))
        except ValueError:
            continue
    return atoms


def get_po4_atoms(gro_atoms):
    out = []
    for i, (_r, a, x, y, z) in enumerate(gro_atoms):
        if a in ("PO4", "P"):
            out.append((i+1, x, y, z))
    return out


def write_leaflets_index(gro_path: Path, out_path: Path):
    atoms = read_gro_atoms(gro_path)
    po4 = get_po4_atoms(atoms)
    if not po4:
        raise SystemExit("No PO4/P atoms found in GRO; cannot build leaflets index.")
    zs = [z for (_i, _x, _y, z) in po4]
    zmed = float(np.median(zs))
    outer = [i for (i, _x, _y, z) in po4 if z >= zmed]
    inner = [i for (i, _x, _y, z) in po4 if z < zmed]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write("[ OuterPO4 ]\n")
        for k, idx in enumerate(outer):
            if k and (k % 15 == 0):
                f.write("\n")
            f.write(f"{idx} ")
        f.write("\n\n[ InnerPO4 ]\n")
        for k, idx in enumerate(inner):
            if k and (k % 15 == 0):
                f.write("\n")
            f.write(f"{idx} ")
        f.write("\n")


def main():
    ap = argparse.ArgumentParser(description="Generate index_leaflets.ndx from GRO")
    ap.add_argument("run_dir", help="Simulation run directory (contains equilibration/")
    ap.add_argument("--gro", dest="gro", help="Path to structure GRO")
    ap.add_argument("--out", dest="out", help="Output ndx path")
    args = ap.parse_args()
    run_dir = Path(args.run_dir)
    gro = Path(args.gro) if args.gro else (run_dir / "equilibration" / "npt" / "npt.gro")
    out = Path(args.out) if args.out else (run_dir / "index_leaflets.ndx")
    if not gro.exists():
        print(f"GRO not found: {gro}", file=sys.stderr)
        sys.exit(2)
    write_leaflets_index(gro, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

