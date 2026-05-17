"""Marginal conservation detection.

Compares ET rankings computed with and without alphabet expansion to identify
positions that are only recognised as conserved when amino acids are grouped by
stereochemical properties — a signature of marginal conservation.

Algorithm
---------
For each alignment column i, compute:

    Δᵢ = score_standard_i - score_expanded_i

A positive Δ means the position scored *lower* (i.e., more conserved) only after
expansion.  Δᵢ is normalised across all positions to produce a Z-score, from
which a one-tailed p-value is derived.  Benjamini-Hochberg FDR correction is
applied across all positions to account for multiple testing.

A position is called *marginally conserved* when Δᵢ > 0 and the BH-corrected
p-value is below the chosen FDR threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import erfc, sqrt, isfinite

import numpy as np

from evotracex.alignment import MSA
from evotracex.subsets import Clade


@dataclass
class MarginalResult:
    position: int
    amino_acid: str
    score_standard: float
    score_expanded: float
    delta: float
    z_score: float
    p_value: float
    adj_p_value: float
    marginal: bool


def detect_marginal_conservation(
    ranking_standard: dict[int, float],
    ranking_expanded: dict[int, float],
    msa: MSA,
    leaf: Clade,
    fdr: float = 0.05,
) -> list[MarginalResult]:
    """Return per-position marginal conservation statistics for one leaf.

    Parameters
    ----------
    ranking_standard:
        ET scores from a run *without* alphabet expansion (lower = more important).
    ranking_expanded:
        ET scores from a run *with* alphabet expansion.
    msa:
        The MSA used for the standard run (supplies sequence content).
    leaf:
        The leaf whose sequence is used to label positions.
    fdr:
        Benjamini-Hochberg FDR threshold for calling a position marginal.

    Returns
    -------
    List of :class:`MarginalResult`, one per non-gap position, sorted by
    *delta* descending (largest improvement from expansion first).
    """
    positions = sorted(ranking_standard)
    leaf_seq = msa.sequences[msa.sequence_indices[leaf.label]]

    deltas = np.array([
        ranking_standard[col] - ranking_expanded[col]
        for col in positions
    ])

    mu = float(np.mean(deltas))
    sigma = float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0

    if sigma < 1e-12:
        z_scores = np.zeros(len(deltas))
    else:
        z_scores = (deltas - mu) / sigma

    p_values = np.array([_z_to_p(float(z)) for z in z_scores])
    adj_p_values = _bh_correct(p_values)

    results: list[MarginalResult] = []
    for i, col in enumerate(positions):
        aa = leaf_seq[col]
        if aa == "-":
            continue
        d = float(deltas[i])
        z = float(z_scores[i])
        p = float(p_values[i])
        adj_p = float(adj_p_values[i])
        results.append(MarginalResult(
            position=col + 1,
            amino_acid=aa,
            score_standard=ranking_standard[col],
            score_expanded=ranking_expanded[col],
            delta=d,
            z_score=z,
            p_value=p,
            adj_p_value=adj_p,
            marginal=(d > 0 and adj_p < fdr),
        ))

    results.sort(key=lambda r: r.delta, reverse=True)
    return results


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _z_to_p(z: float) -> float:
    """One-tailed p-value (upper tail) from a Z-score via the complementary error function."""
    if not isfinite(z):
        return 1.0
    return erfc(z / sqrt(2)) / 2


def _bh_correct(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR correction.

    Returns an array of adjusted p-values in the same order as the input.
    """
    n = len(p_values)
    if n == 0:
        return p_values.copy()
    order = np.argsort(p_values)
    ranked_p = p_values[order]
    # BH adjusted p for rank k (1-indexed): p[k] * n / k
    adjusted = ranked_p * n / np.arange(1, n + 1)
    # Enforce step-up monotonicity: take cumulative minimum from the right
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = np.empty(n)
    result[order] = np.minimum(adjusted, 1.0)
    return result
