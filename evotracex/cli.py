from __future__ import annotations

import argparse
import sys
from pathlib import Path

from evotracex.alignment import MSA
from evotracex.output import write_ranking
from evotracex.scoring import evolutionary_trace
from evotracex.tree import build_clade_network, build_groups


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="evotracex",
        description=(
            "EvoTraceX — Evolutionary Trace with optional alphabet expansion. "
            "Ranks alignment columns by functional importance relative to each leaf."
        ),
    )
    parser.add_argument(
        "msa_file",
        help="Multiple sequence alignment in FASTA format.",
        type=str,
    )
    parser.add_argument(
        "tree_file",
        help="Phylogenetic tree in Newick format (leaf names must match MSA headers).",
        type=str,
    )
    parser.add_argument(
        "-o", "--out",
        help=(
            "Output prefix. Results are written to <out>_et_ranking.tsv. "
            "Defaults to stdout."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "-x", "--plus-aa",
        help="Expand the amino-acid alphabet using stereochemical groupings.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--d0",
        help=(
            "Gaussian bandwidth for subfamily weighting (default: 0.05). "
            "Controls how sharply groups are down-weighted by consensus distance."
        ),
        type=float,
        default=0.05,
        metavar="FLOAT",
    )
    parser.add_argument(
        "--leaf",
        help=(
            "Run analysis for this leaf only (must match an MSA header exactly). "
            "By default all leaves are analysed."
        ),
        type=str,
        default=None,
        metavar="NAME",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    try:
        msa = MSA(args.msa_file, plus_aa=args.plus_aa)
        graph, root, leaves = build_clade_network(msa, args.tree_file)
        groups = build_groups(graph, root)
    except (ValueError, SystemExit) as exc:
        sys.exit(str(exc))

    if args.leaf is not None:
        matching = [lf for lf in leaves if lf.label == args.leaf]
        if not matching:
            sys.exit(
                f"Leaf '{args.leaf}' not found. "
                f"Available: {sorted(lf.label for lf in leaves)}"
            )
        target_leaves = matching
    else:
        target_leaves = leaves

    out = Path(args.out) if args.out else None

    for leaf in target_leaves:
        ranking = evolutionary_trace(
            msa, graph, groups, root, leaf, d0=args.d0
        )
        write_ranking(ranking, msa, leaf, out)


if __name__ == "__main__":
    main()
