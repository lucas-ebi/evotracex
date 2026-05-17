from __future__ import annotations

import sys
from pathlib import Path

from Bio.SeqUtils import seq3

from evotracex.alignment import MSA
from evotracex.subsets import Clade


def write_ranking(
    ranking: dict[int, float],
    msa: MSA,
    leaf: Clade,
    out: Path | None = None,
) -> None:
    """Write the ET ranking for *leaf* as a TSV table.

    Columns are written in ascending score order (most important first).
    Gap positions are skipped. Output goes to *out*_et_ranking.tsv when
    *out* is given, otherwise to stdout.
    """
    leaf_seq = msa.sequences[msa.sequence_indices[leaf.label]]

    rows: list[tuple[int, str, int, float]] = []
    rank = 1
    for col, score in sorted(ranking.items(), key=lambda x: x[1]):
        aa = leaf_seq[col]
        if aa != "-":
            rows.append((rank, seq3(aa), col + 1, score))
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
