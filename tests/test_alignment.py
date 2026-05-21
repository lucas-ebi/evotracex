import math
import numpy as np
import pytest


def test_henikoff_weights_sum_to_one(msa):
    assert math.isclose(sum(msa.weights), 1.0, abs_tol=1e-10)


def test_henikoff_weights_equal_for_uniform_msa(msa):
    # All 5 sequences have unique AAs at every column → equal weights
    expected = 1.0 / msa.size
    assert all(math.isclose(w, expected, abs_tol=1e-10) for w in msa.weights)


def test_henikoff_weights_favour_rare_sequences(fasta_file):
    from evotracex.alignment import MSA

    # Three identical + two identical sequences → minority gets higher weight
    fasta = (
        ">s1\nAAAA\n>s2\nAAAA\n>s3\nAAAA\n"
        ">s4\nLLLL\n>s5\nLLLL\n"
    )
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(fasta)
        path = f.name
    try:
        m = MSA(path)
        # seqA-C (majority) should weigh less than seqD-E (minority)
        assert m.weights[0] < m.weights[3]
        assert math.isclose(m.weights[0], m.weights[1], abs_tol=1e-10)
        assert math.isclose(m.weights[3], m.weights[4], abs_tol=1e-10)
    finally:
        os.unlink(path)


def test_msa_shape(msa):
    assert msa.size == 5
    assert msa.length == 8


def test_collection_fully_conserved_column(msa):
    # Column 0 is all 'A' → exactly one Residue object
    col0 = msa.collection[0]
    assert len(col0) == 1
    assert col0[0].amino_acids == ["A"]
    assert len(col0[0].sequence_indices) == 5


def test_collection_variable_column(msa):
    # Column 6: W, F, Y, K, R → five distinct single-residue groups
    col6 = [r for r in msa.collection[6] if len(r.amino_acids) == 1]
    assert len(col6) == 5


def test_collection_expand_adds_groups(msa, msa_plus):
    # Alphabet expansion should add at least one multi-residue group somewhere
    any_multi = any(
        any(len(r.amino_acids) > 1 for r in msa_plus.collection[col])
        for col in range(msa_plus.length)
    )
    assert any_multi


def test_gap_content_no_gaps(msa):
    # Our fixture has no gaps → all gap Residues have empty sequence_indices
    for gr in msa.gap_content:
        assert len(gr.sequence_indices) == 0


def test_sequence_indices_map(msa):
    assert msa.sequence_indices["seqA"] == 0
    assert msa.sequence_indices["seqE"] == 4


def test_missing_file_exits():
    from evotracex.alignment import MSA
    with pytest.raises(SystemExit):
        MSA("/nonexistent/path.fasta")
