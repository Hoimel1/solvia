#!/usr/bin/env python3
"""Sequence descriptor utilities and lightweight feature engineering tools.

This module provides two main capabilities that downstream analysis scripts can
reuse directly:

* ``compute_sequence_descriptors`` – calculate charge, hydrophobic metrics,
  logP approximation, and helix indicators for individual peptide sequences
  without relying on pre-computed DBAASP metadata.
* ``compute_pca_components`` – perform a principal component analysis on a
  feature matrix so that modelling steps (e.g. ridge regression) can consume
  decorrelated inputs.

The module still exposes a simple CLI (``python scripts/features.py -i
sequences.fasta``) that mirrors the previous behaviour, but the core logic is
split into import-friendly functions so that the PMF analysis pipeline can call
them programmatically.
"""

from __future__ import annotations

import argparse
import multiprocessing
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import peptides

try:  # Optional dependency used for CLI tabular output
    import polars as pl  # type: ignore
except Exception:  # pragma: no cover - polars is optional
    pl = None

__all__ = [
    "compute_sequence_descriptors",
    "compute_pca_components",
    "calculate_features",
]

# Hydrophobic fragment constants (Fauchere & Pliska, Eur. J. Med. Chem. 1983)
# These values are strongly correlated with logP and provide a reproducible,
# training-data-free estimate of peptide lipophilicity.
AA_FRAGMENT_LOGP = {
    "A": 0.31,
    "R": -1.01,
    "N": -0.60,
    "D": -0.77,
    "C": 1.54,
    "Q": -0.22,
    "E": -0.64,
    "G": 0.00,
    "H": 0.13,
    "I": 1.80,
    "L": 1.70,
    "K": -0.99,
    "M": 1.23,
    "F": 1.79,
    "P": 0.72,
    "S": -0.04,
    "T": 0.26,
    "W": 2.25,
    "Y": 0.96,
    "V": 1.22,
}

# Chou-Fasman alpha-helix propensities (normalized to 1.0 being neutral).
CHOU_FASMAN_ALPHA = {
    "A": 1.45,
    "R": 0.79,
    "N": 0.73,
    "D": 1.01,
    "C": 0.77,
    "Q": 1.11,
    "E": 1.51,
    "G": 0.53,
    "H": 1.00,
    "I": 1.08,
    "L": 1.34,
    "K": 1.16,
    "M": 1.20,
    "F": 1.13,
    "P": 0.59,
    "S": 0.79,
    "T": 0.82,
    "W": 1.08,
    "Y": 0.69,
    "V": 1.06,
}

POSITIVE_RESIDUES = {"K", "R", "H"}
NEGATIVE_RESIDUES = {"D", "E"}
AROMATIC_RESIDUES = {"F", "W", "Y"}


def _clean_sequence(seq: str) -> str:
    seq = (seq or "").strip().upper()
    if not seq:
        raise ValueError("Peptide sequence must not be empty")
    return seq


def _safe_hydrophobic_moment(peptide_obj: peptides.Peptide, window: int, angle: int) -> float:
    try:
        return float(peptide_obj.hydrophobic_moment(window=window, angle=angle))
    except Exception:
        return float("nan")


def compute_sequence_descriptors(
    sequence: str,
    *,
    ph: float = 7.4,
    angle: int = 100,
    peptide_obj: peptides.Peptide | None = None,
) -> dict[str, float]:
    """Compute physicochemical descriptors for a peptide sequence.

    Parameters
    ----------
    sequence:
        Amino-acid sequence string (one-letter codes). Case-insensitive.
    ph:
        Solution pH at which to estimate net charge. Defaults to physiological
        pH 7.4.
    angle:
        Rotation angle (degrees) for hydrophobic moment calculation. 100° is
        appropriate for alpha helices.
    peptide_obj:
        Optional pre-instantiated :class:`peptides.Peptide` to avoid repeated
        initialisation when bulk-processing sequences.

    Returns
    -------
    dict
        Mapping of descriptor names to floats suitable for ML ingestion.
    """

    seq = _clean_sequence(sequence)
    pep = peptide_obj or peptides.Peptide(seq)
    length = len(seq)

    charge = float(pep.charge(pH=ph))
    charge_density = charge / length

    hydro = float(pep.hydrophobicity())
    window = max(2, min(length, 11))
    hydromoment = _safe_hydrophobic_moment(pep, window=window, angle=angle)

    logp_total = sum(AA_FRAGMENT_LOGP.get(res, 0.0) for res in seq)
    helix_values = [CHOU_FASMAN_ALPHA.get(res, 1.0) for res in seq]

    descriptors = {
        "sequence_length": float(length),
        "sequence_charge": charge,
        "sequence_charge_density": charge_density,
        "sequence_isoelectric_point": float(pep.isoelectric_point()),
        "sequence_hydrophobicity": hydro,
        "sequence_hydrophobic_moment": hydromoment,
        "sequence_logp": float(logp_total),
        "sequence_logp_per_residue": float(logp_total / length),
        "sequence_helix_propensity": float(np.mean(helix_values)),
        "sequence_helix_fraction": float(np.mean([1.0 if val > 1.0 else 0.0 for val in helix_values])),
        "sequence_positive_fraction": float(sum(res in POSITIVE_RESIDUES for res in seq) / length),
        "sequence_negative_fraction": float(sum(res in NEGATIVE_RESIDUES for res in seq) / length),
        "sequence_aromatic_fraction": float(sum(res in AROMATIC_RESIDUES for res in seq) / length),
    }

    return descriptors


def compute_pca_components(
    matrix: np.ndarray,
    *,
    variance_threshold: float = 0.95,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return decorrelated principal components for a feature matrix.

    Parameters
    ----------
    matrix:
        ``(n_samples, n_features)`` feature matrix. Rows are centred internally.
    variance_threshold:
        Capture at least this cumulative variance fraction. Must be ``0 < x <= 1``.

    Returns
    -------
    scores, loadings, variance_ratio, means
        PCA scores for each sample, column-wise loadings (eigenvectors),
        explained variance ratios per component, and the feature means used to
        centre the data.

    Raises
    ------
    ValueError
        If the input matrix has fewer than two samples or contains NaNs.
    """

    if matrix.ndim != 2:
        raise ValueError("Feature matrix must be 2-dimensional")
    n_samples, n_features = matrix.shape
    if n_samples < 2 or n_features < 1:
        raise ValueError("Need at least two samples and one feature for PCA")
    if not (0.0 < variance_threshold <= 1.0):
        raise ValueError("variance_threshold must be in (0, 1]")
    if np.isnan(matrix).any() or np.isinf(matrix).any():
        raise ValueError("Feature matrix must be finite")

    means = np.mean(matrix, axis=0)
    centred = matrix - means
    cov = np.cov(centred, rowvar=False)

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    total = np.sum(eigvals)
    if total <= 0:
        raise ValueError("Covariance matrix is degenerate")

    variance_ratio = eigvals / total
    cumulative = np.cumsum(variance_ratio)
    n_components = max(1, int(np.searchsorted(cumulative, variance_threshold) + 1))
    eigvecs = eigvecs[:, :n_components]
    variance_ratio = variance_ratio[:n_components]

    scores = centred @ eigvecs
    return scores, eigvecs, variance_ratio, means


def _calc(peptide: str) -> dict[str, float]:
    seq = _clean_sequence(peptide)
    pep = peptides.Peptide(seq)
    descriptors = compute_sequence_descriptors(seq, peptide_obj=pep)

    res: dict[str, float] = {"seq": seq}
    res.update(descriptors)

    try:
        res["aidx"] = float(pep.aliphatic_index())
    except ZeroDivisionError:
        res["aidx"] = 0.0

    res["boman"] = float(pep.boman())
    res["charge"] = descriptors.get("sequence_charge", float("nan"))
    res["hm"] = descriptors.get("sequence_hydrophobic_moment", float("nan"))
    res["hp"] = descriptors.get("sequence_hydrophobicity", float("nan"))
    res["iep"] = descriptors.get("sequence_isoelectric_point", float("nan"))
    res["iidx"] = float(pep.instability_index())
    res["mol"] = float(pep.molecular_weight()) / len(seq)
    res["mz"] = float(pep.mz()) / len(seq)

    for idx, factor in enumerate(pep.atchley_factors(), start=1):
        res[f"AF_{idx}"] = float(factor)

    for idx, component in enumerate(pep.fasgai_vectors(), start=1):
        res[f"FG_{idx}"] = float(component)

    for idx, vhse in enumerate(pep.vhse_scales(), start=1):
        res[f"VHSE_{idx}"] = float(vhse)

    for idx, zscale in enumerate(pep.z_scales(), start=1):
        res[f"ZS_{idx}"] = float(zscale)

    for idx, kidera in enumerate(pep.kidera_factors(), start=1):
        res[f"KF_{idx}"] = float(kidera)

    return res


def calculate_features(sequences: Sequence[str], n_jobs: int = 8) -> "pl.DataFrame":
    """Return a Polars DataFrame with descriptors for many sequences."""

    if pl is None:  # pragma: no cover - CLI only
        raise RuntimeError("polars is required for DataFrame output")

    with multiprocessing.Pool(processes=max(1, int(n_jobs))) as pool:
        features: list[dict[str, float]] = pool.map(_calc, sequences)
    return pl.from_dicts(features)


def _read_fasta(path: Path) -> list[str]:  # pragma: no cover - exercised via CLI
    seqs = []
    current = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
                continue
            current.append(line)
        if current:
            seqs.append("".join(current))
    return seqs


def _load_sequences(path: Path) -> list[str]:  # pragma: no cover - CLI helper
    ext = path.suffix.lower()
    if ext in {".fa", ".faa", ".fasta"}:
        return _read_fasta(path)
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _cli() -> None:  # pragma: no cover - thin wrapper
    parser = argparse.ArgumentParser(description="Compute peptide descriptors")
    parser.add_argument("-i", "--infile", required=True, help="FASTA or newline-separated sequence file")
    parser.add_argument("-o", "--output", help="Optional CSV output path")
    parser.add_argument("-p", "--processes", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--ph", type=float, default=7.4, help="Solution pH for net charge computation")
    args = parser.parse_args()

    sequences = _load_sequences(Path(args.infile))
    if not sequences:
        raise SystemExit("No sequences found")

    if pl is None:
        raise SystemExit("polars is required for CLI tabular output; pip install polars")

    df = calculate_features(sequences, n_jobs=args.processes)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(args.output)
        print(f"Wrote {len(df)} rows to {args.output}")
    else:
        print(df)


if __name__ == "__main__":  # pragma: no cover
    _cli()
