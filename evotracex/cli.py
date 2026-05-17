from __future__ import annotations

import argparse
import sys
from pathlib import Path

from evotracex.alignment import MSA
from evotracex.marginal import detect_marginal_conservation
from evotracex.output import write_marginal_ranking, write_ranking
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
            "Output prefix. Results are written to <out>_et_ranking.tsv "
            "(or <out>_marginal.tsv with --marginal). Defaults to stdout."
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
    parser.add_argument(
        "--marginal",
        help=(
            "Detect marginally conserved positions by comparing ET rankings "
            "with and without alphabet expansion. Runs both analyses automatically "
            "and applies Benjamini-Hochberg FDR correction. "
            "Writes <out>_marginal.tsv (or stdout). Implies --plus-aa."
        ),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--fdr",
        help="FDR threshold for marginal conservation calls (default: 0.05).",
        type=float,
        default=0.05,
        metavar="FLOAT",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    out = Path(args.out) if args.out else None

    if args.marginal:
        _run_marginal(args, out)
    else:
        _run_standard(args, out)


def _run_standard(args: argparse.Namespace, out: Path | None) -> None:
    try:
        msa = MSA(args.msa_file, plus_aa=args.plus_aa)
        graph, root, leaves = build_clade_network(msa, args.tree_file)
        groups = build_groups(graph, root)
    except (ValueError, SystemExit) as exc:
        sys.exit(str(exc))

    target_leaves = _select_leaves(leaves, args.leaf)
    for leaf in target_leaves:
        ranking = evolutionary_trace(msa, graph, groups, root, leaf, d0=args.d0)
        write_ranking(ranking, msa, leaf, out)


def _run_marginal(args: argparse.Namespace, out: Path | None) -> None:
    try:
        msa_std = MSA(args.msa_file, plus_aa=False)
        msa_exp = MSA(args.msa_file, plus_aa=True)
        graph_std, root_std, leaves_std = build_clade_network(msa_std, args.tree_file)
        graph_exp, root_exp, leaves_exp = build_clade_network(msa_exp, args.tree_file)
        groups_std = build_groups(graph_std, root_std)
        groups_exp = build_groups(graph_exp, root_exp)
    except (ValueError, SystemExit) as exc:
        sys.exit(str(exc))

    leaves_exp_by_label = {lf.label: lf for lf in leaves_exp}
    target_std = _select_leaves(leaves_std, args.leaf)

    for leaf_std in target_std:
        leaf_exp = leaves_exp_by_label[leaf_std.label]
        ranking_std = evolutionary_trace(
            msa_std, graph_std, groups_std, root_std, leaf_std, d0=args.d0
        )
        ranking_exp = evolutionary_trace(
            msa_exp, graph_exp, groups_exp, root_exp, leaf_exp, d0=args.d0
        )
        results = detect_marginal_conservation(
            ranking_std, ranking_exp, msa_std, leaf_std, fdr=args.fdr
        )
        write_marginal_ranking(results, leaf_std, out)


def _select_leaves(leaves: list, label: str | None) -> list:
    if label is None:
        return leaves
    matching = [lf for lf in leaves if lf.label == label]
    if not matching:
        sys.exit(
            f"Leaf '{label}' not found. "
            f"Available: {sorted(lf.label for lf in leaves)}"
        )
    return matching


if __name__ == "__main__":
    main()
