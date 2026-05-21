from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path
from typing import Callable

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
        "-x", "--expand",
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
        "--ref",
        help=(
            "Reference sequence for analysis and position numbering "
            "(must match an MSA header exactly). "
            "Restricts output to this sequence; alignment columns where it has "
            "a gap are omitted. Defaults to the first sequence in the MSA."
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
            "Writes <out>_marginal.tsv (or stdout). Implies --expand."
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


_SPINNER = list("|/-\\")
_BAR_WIDTH = 28


def _log(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def _make_progress(label: str) -> Callable[[int, int], None]:
    spinner = itertools.cycle(_SPINNER)

    def _callback(current: int, total: int) -> None:
        pct = current / total
        filled = int(_BAR_WIDTH * pct)
        bar = "█" * filled + "░" * (_BAR_WIDTH - filled)
        sys.stderr.write(f"\r  {next(spinner)} {label}  [{bar}] {current}/{total}")
        sys.stderr.flush()

    return _callback


def _end_progress() -> None:
    sys.stderr.write("\n")
    sys.stderr.flush()


def main() -> None:
    args = build_parser().parse_args()
    out = Path(args.out) if args.out else None

    t0 = time.perf_counter()
    if args.marginal:
        _run_marginal(args, out)
    else:
        _run_standard(args, out)
    _log(f"Done in {time.perf_counter() - t0:.1f}s.")


def _build_ref_pos_map(msa: MSA, ref_label: str) -> dict[int, int]:
    if ref_label not in msa.sequence_indices:
        sys.exit(
            f"Reference sequence '{ref_label}' not found. "
            f"Available: {sorted(msa.sequence_indices)}"
        )
    ref_idx = msa.sequence_indices[ref_label]
    ref_seq = msa.sequences[ref_idx]
    pos = 0
    ref_pos_map: dict[int, int] = {}
    for col, aa in enumerate(ref_seq):
        if aa != "-":
            pos += 1
            ref_pos_map[col] = pos
    return ref_pos_map


def _run_standard(args: argparse.Namespace, out: Path | None) -> None:
    _log("Loading MSA and tree…")
    try:
        msa = MSA(args.msa_file, expand=args.expand)
        graph, root, leaves = build_clade_network(msa, args.tree_file)
        groups = build_groups(graph, root)
    except (ValueError, SystemExit) as exc:
        sys.exit(str(exc))
    _log(f"  {msa.size} sequences · {msa.length} columns · {len(leaves)} leaves")

    ref = args.ref or msa.headers[0]
    ref_pos_map = _build_ref_pos_map(msa, ref)
    _log(f"  reference: {ref} ({len(ref_pos_map)} residues / {msa.length} alignment columns)")
    target_leaves = _select_leaves(leaves, ref)
    for leaf in target_leaves:
        _log(f"Running ET for {leaf.label}…")
        ranking = evolutionary_trace(
            msa, graph, groups, root, leaf, d0=args.d0,
            progress=_make_progress("ET"),
        )
        _end_progress()
        write_ranking(ranking, msa, leaf, out, ref_pos_map=ref_pos_map)


def _run_marginal(args: argparse.Namespace, out: Path | None) -> None:
    _log("Loading MSA and tree…")
    try:
        msa_std = MSA(args.msa_file, expand=False)
        msa_exp = MSA(args.msa_file, expand=True)
        graph_std, root_std, leaves_std = build_clade_network(msa_std, args.tree_file)
        graph_exp, root_exp, leaves_exp = build_clade_network(msa_exp, args.tree_file)
        groups_std = build_groups(graph_std, root_std)
        groups_exp = build_groups(graph_exp, root_exp)
    except (ValueError, SystemExit) as exc:
        sys.exit(str(exc))
    _log(f"  {msa_std.size} sequences · {msa_std.length} columns · {len(leaves_std)} leaves")

    ref = args.ref or msa_std.headers[0]
    ref_pos_map = _build_ref_pos_map(msa_std, ref)
    _log(f"  reference: {ref} ({len(ref_pos_map)} residues / {msa_std.length} alignment columns)")
    leaves_exp_by_label = {lf.label: lf for lf in leaves_exp}
    target_std = _select_leaves(leaves_std, ref)

    for leaf_std in target_std:
        _log(f"Running ET for {leaf_std.label}…")
        leaf_exp = leaves_exp_by_label[leaf_std.label]
        ranking_std = evolutionary_trace(
            msa_std, graph_std, groups_std, root_std, leaf_std, d0=args.d0,
            progress=_make_progress("standard"),
        )
        _end_progress()
        ranking_exp = evolutionary_trace(
            msa_exp, graph_exp, groups_exp, root_exp, leaf_exp, d0=args.d0,
            progress=_make_progress("expanded"),
        )
        _end_progress()
        results = detect_marginal_conservation(
            ranking_std, ranking_exp, msa_std, leaf_std, fdr=args.fdr
        )
        write_marginal_ranking(results, leaf_std, out, ref_pos_map=ref_pos_map)


def _select_leaves(leaves: list, ref: str | None) -> list:
    if ref is None:
        return leaves
    matching = [lf for lf in leaves if lf.label == ref]
    if not matching:
        sys.exit(
            f"Reference sequence '{ref}' not found. "
            f"Available: {sorted(lf.label for lf in leaves)}"
        )
    return matching


if __name__ == "__main__":
    main()
