import sys
import pytest


def _run_main(args: list[str]) -> None:
    """Invoke cli.main() with the given sys.argv."""
    from evotracex.cli import main
    sys.argv = ["evotracex"] + args
    main()


def test_cli_runs_without_error(fasta_file, tree_file, capsys):
    _run_main([fasta_file, tree_file])
    out = capsys.readouterr().out
    assert "rank" in out
    assert "residue" in out


def test_cli_output_contains_all_leaves(fasta_file, tree_file, capsys):
    _run_main([fasta_file, tree_file])
    out = capsys.readouterr().out
    for name in ("seqA", "seqB", "seqC", "seqD", "seqE"):
        assert name in out


def test_cli_leaf_flag_restricts_output(fasta_file, tree_file, capsys):
    _run_main([fasta_file, tree_file, "--leaf", "seqA"])
    out = capsys.readouterr().out
    assert "seqA" in out
    assert "seqB" not in out


def test_cli_invalid_leaf_exits(fasta_file, tree_file):
    with pytest.raises(SystemExit):
        _run_main([fasta_file, tree_file, "--leaf", "nonexistent"])


def test_cli_out_writes_tsv(fasta_file, tree_file, tmp_path):
    prefix = str(tmp_path / "results")
    _run_main([fasta_file, tree_file, "--out", prefix])
    tsv = tmp_path / "results_et_ranking.tsv"
    assert tsv.exists()
    lines = tsv.read_text().splitlines()
    headers = [l for l in lines if l.startswith("rank\t")]
    assert headers, "TSV must contain a header row"
    # Should have header + at least one data row per leaf
    assert len(lines) > 5


def test_cli_out_tsv_columns(fasta_file, tree_file, tmp_path):
    prefix = str(tmp_path / "results")
    _run_main([fasta_file, tree_file, "--out", prefix, "--leaf", "seqA"])
    tsv = tmp_path / "results_et_ranking.tsv"
    content = tsv.read_text()
    assert "rank\tresidue\tposition\tscore" in content


def test_cli_plus_aa_flag_accepted(fasta_file, tree_file, capsys):
    _run_main([fasta_file, tree_file, "--plus-aa", "--leaf", "seqA"])
    out = capsys.readouterr().out
    assert "seqA" in out


def test_cli_d0_flag_accepted(fasta_file, tree_file, capsys):
    _run_main([fasta_file, tree_file, "--d0", "0.5", "--leaf", "seqA"])
    out = capsys.readouterr().out
    assert "rank" in out


def test_cli_msa_mismatch_exits(tree_file, tmp_path):
    bad_fasta = tmp_path / "bad.fasta"
    bad_fasta.write_text(">wrong_name\nAAAAAAAA\n")
    with pytest.raises(SystemExit):
        _run_main([str(bad_fasta), tree_file])


def test_cli_missing_msa_exits(tree_file):
    with pytest.raises(SystemExit):
        _run_main(["/no/such/file.fasta", tree_file])
