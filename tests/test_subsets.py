import math
import numpy as np
import pytest


def test_weight_full_subset(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    assert math.isclose(full.weight, 1.0, abs_tol=1e-10)


def test_weight_empty_subset(msa):
    from evotracex.subsets import Subset
    empty = Subset(msa, [], label="none")
    assert empty.weight == 0.0


def test_conditional_weight_self(msa):
    from evotracex.subsets import Subset
    sub = Subset(msa, [0, 1, 2], label="abc")
    assert math.isclose(sub.conditional_weight(sub), 1.0, abs_tol=1e-10)


def test_conditional_weight_disjoint(msa):
    from evotracex.subsets import Subset
    a = Subset(msa, [0, 1], label="a")
    b = Subset(msa, [2, 3], label="b")
    assert a.conditional_weight(b) == 0.0


def test_entropy_fully_conserved_column(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    # Column 0 is all 'A' → entropy should be 0
    assert math.isclose(full.entropy(0), 0.0, abs_tol=1e-10)
    assert math.isclose(full.entropy(1), 0.0, abs_tol=1e-10)


def test_entropy_variable_column_positive(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    # Columns 2-7 each have 5 distinct amino acids → entropy > 0
    for col in range(2, 8):
        assert full.entropy(col) > 0.0


def test_entropy_variable_column_value(msa):
    """With 5 equal-weight sequences each having a unique AA, entropy = log(5)."""
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    expected = math.log(5)
    for col in range(2, 8):
        assert math.isclose(full.entropy(col), expected, abs_tol=1e-8)


def test_entropy_below_maximum(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    max_entropy = math.log(20)
    for col in range(msa.length):
        assert full.entropy(col) <= max_entropy + 1e-10


def test_consensus_conserved_column(msa):
    from evotracex.subsets import Subset
    full = Subset(msa, range(msa.size), label="all")
    con = full.consensus
    assert con[0] == "A"
    assert con[1] == "A"


def test_residue_label_single(msa):
    residue = msa.collection[0][0]  # Ala at col 0
    assert residue.label == "Ala1"


def test_residue_label_gap(msa):
    assert msa.gap_content[0].label == "-1"


def test_clade_children_initialised_empty():
    """Clade.children starts as [] so len() never raises TypeError."""
    from evotracex.alignment import MSA
    from evotracex.subsets import Clade
    import tempfile, os
    fasta = ">s1\nAAAA\n>s2\nLLLL\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(fasta)
        path = f.name
    try:
        m = MSA(path)
        c = Clade(m, [0, 1])
        assert c.children == []
        assert len(c.children) == 0
    finally:
        os.unlink(path)
