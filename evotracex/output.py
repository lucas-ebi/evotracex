from __future__ import annotations

import sys
from pathlib import Path

from Bio.SeqUtils import seq3

from evotracex.alignment import MSA
from evotracex.marginal import MarginalResult
from evotracex.subsets import Clade


def write_ranking(
    ranking: dict[int, float],
    msa: MSA,
    leaf: Clade,
    out: Path | None = None,
    ref_pos_map: dict[int, int] | None = None,
) -> None:
    """Write the ET ranking for *leaf* as a TSV table.

    Columns are written in ascending score order (most important first).
    Gap positions (in the leaf sequence) are skipped.  When *ref_pos_map*
    is provided, positions are numbered by the reference sequence and
    columns where the reference has a gap are also skipped.
    Output goes to *out*_et_ranking.tsv when *out* is given, otherwise
    to stdout.
    """
    leaf_seq = msa.sequences[msa.sequence_indices[leaf.label]]

    rows: list[tuple[int, str, int, float]] = []
    rank = 1
    for col, score in sorted(ranking.items(), key=lambda x: x[1]):
        aa = leaf_seq[col]
        if aa == "-":
            continue
        if ref_pos_map is not None and col not in ref_pos_map:
            continue
        position = ref_pos_map[col] if ref_pos_map is not None else col + 1
        rows.append((rank, seq3(aa), position, score))
        rank += 1

    header = "\t".join(["rank", "residue", "position", "score"])
    lines = [header] + [
        f"{r}\t{res}\t{pos}\t{score:.6f}" for r, res, pos, score in rows
    ]
    text = "\n".join(lines) + "\n"

    if out is None:
        sys.stdout.write(f"# {leaf.label}\n")
        sys.stdout.write(text)
    else:
        path = Path(str(out) + "_et_ranking.tsv")
        mode = "a" if path.exists() else "w"
        with path.open(mode) as fh:
            fh.write(f"# {leaf.label}\n")
            fh.write(text)


def write_marginal_ranking(
    results: list[MarginalResult],
    leaf: Clade,
    out: Path | None = None,
    ref_pos_map: dict[int, int] | None = None,
) -> None:
    """Write the marginal conservation table for *leaf* as a TSV.

    Rows are sorted by delta descending (largest improvement from expansion first).
    Only non-gap positions are included.  When *ref_pos_map* is provided,
    positions are numbered by the reference sequence and columns where the
    reference has a gap are omitted.
    Output goes to *out*_marginal.tsv when *out* is given, otherwise to stdout.
    """
    header = "\t".join([
        "position", "residue",
        "score_standard", "score_expanded", "delta",
        "z_score", "p_value", "adj_p_value", "marginal",
    ])
    lines = [header]
    for r in results:
        col = r.position - 1
        if ref_pos_map is not None and col not in ref_pos_map:
            continue
        position = ref_pos_map[col] if ref_pos_map is not None else r.position
        lines.append(
            f"{position}\t{seq3(r.amino_acid)}\t"
            f"{r.score_standard:.6f}\t{r.score_expanded:.6f}\t{r.delta:.6f}\t"
            f"{r.z_score:.4f}\t{r.p_value:.4e}\t{r.adj_p_value:.4e}\t"
            f"{'TRUE' if r.marginal else 'FALSE'}"
        )
    text = "\n".join(lines) + "\n"

    if out is None:
        sys.stdout.write(f"# {leaf.label}\n")
        sys.stdout.write(text)
    else:
        path = Path(str(out) + "_marginal.tsv")
        mode = "a" if path.exists() else "w"
        with path.open(mode) as fh:
            fh.write(f"# {leaf.label}\n")
            fh.write(text)
