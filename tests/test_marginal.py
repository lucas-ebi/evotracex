"""Tests for marginal conservation detection."""

import math
import pytest

from evotracex.marginal import MarginalResult, _bh_correct, _z_to_p, detect_marginal_conservation
import numpy as np


# ---------------------------------------------------------------------------
# Unit tests for statistical helpers
# ---------------------------------------------------------------------------

def test_z_to_p_zero():
    # Z=0: p = 0.5 (one-tailed)
    assert math.isclose(_z_to_p(0.0), 0.5, abs_tol=1e-10)


def test_z_to_p_large_positive():
    # Very large Z → p → 0
    assert _z_to_p(10.0) < 1e-10


def test_z_to_p_negative():
    # Negative Z → p > 0.5 (not significant in upper tail)
    assert _z_to_p(-2.0) > 0.5


def test_z_to_p_non_finite():
    assert _z_to_p(float("inf")) == 1.0
    assert _z_to_p(float("-inf")) == 1.0
    assert _z_to_p(float("nan")) == 1.0


def test_bh_correct_sorted_input():
    # Known example: 4 p-values, n=4
    p = np.array([0.01, 0.04, 0.10, 0.20])
    adj = _bh_correct(p)
    # All adjusted p-values must be >= raw p-values
    assert all(adj[i] >= p[i] - 1e-12 for i in range(len(p)))


def test_bh_correct_monotone():
    # BH adjusted p-values, when sorted by original p-value, should be non-decreasing
    rng = np.random.default_rng(42)
    p = rng.uniform(0, 1, 20)
    adj = _bh_correct(p)
    sorted_adj = adj[np.argsort(p)]
    assert all(sorted_adj[i] <= sorted_adj[i + 1] + 1e-12 for i in range(len(sorted_adj) - 1))


def test_bh_correct_capped_at_one():
    p = np.array([0.9, 0.95, 0.99])
    adj = _bh_correct(p)
    assert all(a <= 1.0 for a in adj)


def test_bh_correct_empty():
    adj = _bh_correct(np.array([]))
    assert len(adj) == 0


# ---------------------------------------------------------------------------
# Integration tests for detect_marginal_conservation
# ---------------------------------------------------------------------------

def _make_rankings(n_cols: int, marginal_cols: list[int], boost: float = 0.5):
    """Synthetic rankings where *marginal_cols* improve by *boost* with expansion."""
    std = {col: 1.0 + col * 0.01 for col in range(n_cols)}
    exp = {col: v - (boost if col in marginal_cols else 0.0) for col, v in std.items()}
    return std, exp


def test_detect_returns_one_result_per_nongap_position(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    std = {col: 1.0 + col * 0.01 for col in range(msa.length)}
    exp = {col: v for col, v in std.items()}

    # Use the first leaf from the MSA fixture
    from evotracex.subsets import Clade
    leaf = Clade(msa, [0], label="seqA")
    results = detect_marginal_conservation(std, exp, msa, leaf)
    # seqA has no gaps in our fixture → should have msa.length results
    assert len(results) == msa.length


def test_detect_sorted_by_delta_descending(fasta_file, tree_file):
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade
    msa = MSA(fasta_file)
    leaf = Clade(msa, list(range(msa.size)), label="seqA")
    std, exp = _make_rankings(msa.length, marginal_cols=[2, 5])
    results = detect_marginal_conservation(std, exp, msa, leaf)
    deltas = [r.delta for r in results]
    assert deltas == sorted(deltas, reverse=True)


def test_detect_marginal_boosted_positions_rank_highest(fasta_file):
    """Boosted columns should have the largest z-scores.

    Note: with only 8 columns, BH correction at FDR=0.05 cannot reach significance
    (Z ≈ 1.62 for 2/8 boosted gives adj_p ≈ 0.42).  We test the detection logic
    directly: the boosted positions must be the top-ranked candidates.
    """
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade

    msa = MSA(fasta_file)
    leaf = Clade(msa, list(range(msa.size)), label="seqA")

    std, exp = _make_rankings(msa.length, marginal_cols=[2, 5], boost=2.0)
    results = detect_marginal_conservation(std, exp, msa, leaf, fdr=0.05)

    # Sorted by delta descending → top 2 should be columns 2 and 5 (positions 3 and 6)
    top_positions = {r.position for r in results[:2]}
    assert 3 in top_positions
    assert 6 in top_positions
    # Their z-scores must be positive (above mean)
    assert all(r.z_score > 0 for r in results[:2])


def test_detect_no_marginal_when_no_boost(fasta_file):
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade

    msa = MSA(fasta_file)
    leaf = Clade(msa, list(range(msa.size)), label="seqA")
    # All deltas = 0 → sigma = 0 → z = 0 → p = 0.5 → adj_p = 0.5 → no marginal calls
    std = {col: 1.5 for col in range(msa.length)}
    exp = {col: 1.5 for col in range(msa.length)}
    results = detect_marginal_conservation(std, exp, msa, leaf, fdr=0.05)
    assert not any(r.marginal for r in results)


def test_detect_delta_values_correct(fasta_file):
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade
    import math as _math

    msa = MSA(fasta_file)
    leaf = Clade(msa, [0], label="seqA")
    std = {col: float(col) + 2.0 for col in range(msa.length)}
    exp = {col: float(col) + 1.0 for col in range(msa.length)}
    results = detect_marginal_conservation(std, exp, msa, leaf)
    for r in results:
        col = r.position - 1
        s, e = std[col], exp[col]
        expected_delta = _math.atan((s - e) / (s + e))
        assert math.isclose(r.delta, expected_delta, abs_tol=1e-10)


def test_detect_p_values_in_range(fasta_file):
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade

    msa = MSA(fasta_file)
    leaf = Clade(msa, list(range(msa.size)), label="seqA")
    std, exp = _make_rankings(msa.length, marginal_cols=[0, 1], boost=0.3)
    results = detect_marginal_conservation(std, exp, msa, leaf)
    for r in results:
        assert 0.0 <= r.p_value <= 1.0
        assert 0.0 <= r.adj_p_value <= 1.0


# ---------------------------------------------------------------------------
# CLI integration tests for --marginal
# ---------------------------------------------------------------------------

def test_cli_marginal_flag_runs(fasta_file, tree_file, capsys):
    import sys
    from evotracex.cli import main
    sys.argv = ["evotracex", fasta_file, tree_file, "--marginal"]
    main()
    out = capsys.readouterr().out
    assert "position" in out
    assert "delta" in out
    assert "marginal" in out


def test_cli_marginal_writes_tsv(fasta_file, tree_file, tmp_path):
    import sys
    from evotracex.cli import main
    prefix = str(tmp_path / "out")
    sys.argv = ["evotracex", fasta_file, tree_file, "--marginal", "--out", prefix]
    main()
    tsv = tmp_path / "out_marginal.tsv"
    assert tsv.exists()
    content = tsv.read_text()
    assert "position\tresidue\tscore_standard\tscore_expanded\tdelta" in content


def test_cli_marginal_ref_flag(fasta_file, tree_file, capsys):
    import sys
    from evotracex.cli import main
    sys.argv = ["evotracex", fasta_file, tree_file, "--marginal", "--ref", "seqA"]
    main()
    out = capsys.readouterr().out
    assert "seqA" in out
    assert "seqB" not in out


def test_cli_marginal_fdr_flag_accepted(fasta_file, tree_file, capsys):
    import sys
    from evotracex.cli import main
    sys.argv = ["evotracex", fasta_file, tree_file, "--marginal", "--fdr", "0.1", "--ref", "seqA"]
    main()
    out = capsys.readouterr().out
    assert "position" in out
