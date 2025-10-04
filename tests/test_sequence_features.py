import sys
from pathlib import Path

import numpy as np
import peptides
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.features import compute_pca_components, compute_sequence_descriptors


def test_compute_sequence_descriptors_matches_peptides_defaults():
    seq = "ACDE"
    desc = compute_sequence_descriptors(seq)
    pep = peptides.Peptide(seq)

    assert desc["sequence_length"] == pytest.approx(len(seq))
    assert desc["sequence_charge"] == pytest.approx(pep.charge(pH=7.4))
    assert desc["sequence_hydrophobicity"] == pytest.approx(pep.hydrophobicity())
    assert desc["sequence_hydrophobic_moment"] == pytest.approx(
        pep.hydrophobic_moment(window=min(len(seq), 11), angle=100)
    )
    assert desc["sequence_logp_per_residue"] == pytest.approx(
        desc["sequence_logp"] / desc["sequence_length"]
    )
    assert 0.0 <= desc["sequence_helix_fraction"] <= 1.0


def test_compute_pca_components_reconstructs_full_variance():
    matrix = np.array(
        [
            [1.0, 0.0, 2.0],
            [0.0, 1.0, 2.0],
            [1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0],
        ],
        dtype=float,
    )

    scores, loadings, variance_ratio, means = compute_pca_components(
        matrix, variance_threshold=1.0
    )

    assert scores.shape[0] == matrix.shape[0]
    assert loadings.shape[0] == matrix.shape[1]
    assert pytest.approx(np.sum(variance_ratio)) == 1.0

    reconstructed = scores @ loadings.T + means
    assert np.allclose(reconstructed, matrix)
