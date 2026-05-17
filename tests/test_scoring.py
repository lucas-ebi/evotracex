import math
import pytest

from evotracex.constants import ALPHABET
from evotracex.scoring import similarity


# ---------------------------------------------------------------------------
# similarity()
# ---------------------------------------------------------------------------

def test_similarity_identical():
    for aa in ALPHABET:
        assert similarity(aa, aa) == 1


def test_similarity_gap_vs_residue_returns_zero():
    assert similarity("-", "A") == 0
    assert similarity("A", "-") == 0


def test_similarity_gap_gap_returns_one():
    # The identity check (aa1 == aa2) fires before the gap check, so gap-gap = 1.
    # This matches the original ET behaviour and is harmless: two groups that both
    # have no coverage at a position "agree" and contribute no distance.
    assert similarity("-", "-") == 1


def test_similarity_no_off_diagonal_above_threshold():
    # With BLOSUM62_THRESHOLD ≈ 5.8, no off-diagonal pair exceeds the threshold.
    # similarity() is effectively an identity check for amino-acid pairs.
    from evotracex.constants import ALPHABET
    for i, aa1 in enumerate(ALPHABET):
        for aa2 in ALPHABET[i + 1:]:
            assert similarity(aa1, aa2) == 0, f"{aa1}-{aa2} should be 0"


def test_similarity_symmetric():
    for i, aa1 in enumerate(ALPHABET):
        for aa2 in ALPHABET[i:]:
            assert similarity(aa1, aa2) == similarity(aa2, aa1)


# ---------------------------------------------------------------------------
# evolutionary_trace()
# ---------------------------------------------------------------------------

def test_et_ranking_conserved_columns_rank_lowest(fasta_file, tree_file):
    from evotracex.alignment import MSA
    from evotracex.scoring import evolutionary_trace
    from evotracex.tree import build_clade_network, build_groups

    msa = MSA(fasta_file)
    graph, root, leaves = build_clade_network(msa, tree_file)
    groups = build_groups(graph, root)

    for leaf in leaves:
        ranking = evolutionary_trace(msa, graph, groups, root, leaf)
        # Columns 0 and 1 (all 'A') must have the lowest scores
        min_score = min(ranking.values())
        assert math.isclose(ranking[0], min_score, abs_tol=1e-8)
        assert math.isclose(ranking[1], min_score, abs_tol=1e-8)


def test_et_ranking_returns_all_columns(fasta_file, tree_file):
    from evotracex.alignment import MSA
    from evotracex.scoring import evolutionary_trace
    from evotracex.tree import build_clade_network, build_groups

    msa = MSA(fasta_file)
    graph, root, leaves = build_clade_network(msa, tree_file)
    groups = build_groups(graph, root)
    leaf = leaves[0]

    ranking = evolutionary_trace(msa, graph, groups, root, leaf)
    assert set(ranking.keys()) == set(range(msa.length))


def test_et_ranking_scores_non_negative(fasta_file, tree_file):
    from evotracex.alignment import MSA
    from evotracex.scoring import evolutionary_trace
    from evotracex.tree import build_clade_network, build_groups

    msa = MSA(fasta_file)
    graph, root, leaves = build_clade_network(msa, tree_file)
    groups = build_groups(graph, root)

    for leaf in leaves:
        ranking = evolutionary_trace(msa, graph, groups, root, leaf)
        assert all(v >= 0.0 for v in ranking.values())


def test_et_ranking_custom_d0(fasta_file, tree_file):
    """d0 parameter must be accepted and affect scores (different from default)."""
    from evotracex.alignment import MSA
    from evotracex.scoring import evolutionary_trace
    from evotracex.tree import build_clade_network, build_groups

    msa = MSA(fasta_file)
    graph, root, leaves = build_clade_network(msa, tree_file)
    groups = build_groups(graph, root)
    leaf = leaves[0]

    r_default = evolutionary_trace(msa, graph, groups, root, leaf, d0=0.05)
    r_large = evolutionary_trace(msa, graph, groups, root, leaf, d0=1.0)
    # With d0=1.0 (wider Gaussian) the weights differ → at least some scores differ
    assert any(
        not math.isclose(r_default[col], r_large[col], rel_tol=1e-6)
        for col in range(msa.length)
    )
