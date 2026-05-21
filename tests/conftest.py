"""Shared fixtures: a tiny 5-sequence × 8-column MSA and matching Newick tree.

Column layout
-------------
Col 0-1 : all 'A'  → fully conserved, ET score should be minimal
Col 2-5 : each of the 5 sequences has a distinct amino acid → moderately variable
Col 6-7 : five distinct amino acids → maximally variable (within this alphabet)

With equal Henikoff weights (all sequences contribute equally), the entropy at
cols 0-1 is 0, and at cols 2-7 it is log(5) ≈ 1.609.
"""

import pytest

# 5 sequences × 8 columns, no gaps
FASTA_TEXT = """\
>seqA
AACDEFWD
>seqB
AAGHIKFE
>seqC
AALMNPYH
>seqD
AAQRSTKN
>seqE
AAVWYARQ
"""

# Bifurcating tree with two internal nodes and sensible branch lengths.
# ((seqA,seqB),(seqC,(seqD,seqE)));
NEWICK_TEXT = "((seqA:0.1,seqB:0.1):0.2,(seqC:0.1,(seqD:0.1,seqE:0.1):0.1):0.2);"


@pytest.fixture
def fasta_file(tmp_path):
    p = tmp_path / "test.fasta"
    p.write_text(FASTA_TEXT)
    return str(p)


@pytest.fixture
def tree_file(tmp_path):
    p = tmp_path / "test.nwk"
    p.write_text(NEWICK_TEXT)
    return str(p)


@pytest.fixture
def msa(fasta_file):
    from evotracex.alignment import MSA
    return MSA(fasta_file)


@pytest.fixture
def msa_plus(fasta_file):
    from evotracex.alignment import MSA
    return MSA(fasta_file, expand=True)
