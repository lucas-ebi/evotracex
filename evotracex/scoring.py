from __future__ import annotations

from typing import Callable

import numpy as np
import networkx as nx

from evotracex.constants import ALPHABET, BLOSUM62, BLOSUM62_THRESHOLD
from evotracex.alignment import MSA
from evotracex.subsets import Clade

# Precomputed binary similarity for every AA pair — O(1) lookup instead of
# iterating BLOSUM62 keys on every call.
_SIM: dict[tuple[str, str], int] = {}
for _a in ALPHABET:
    for _b in ALPHABET:
        if _a == _b:
            _SIM[(_a, _b)] = 1
        else:
            try:
                _s = BLOSUM62[(_a, _b)]
            except (KeyError, IndexError):
                _s = BLOSUM62[(_b, _a)]
            _SIM[(_a, _b)] = int(_s >= BLOSUM62_THRESHOLD)


def similarity(aa1: str, aa2: str) -> int:
    """Binary similarity: 1 if identical or BLOSUM62 score >= mean diagonal, else 0."""
    if aa1 == aa2:
        return 1
    if "-" in (aa1, aa2):
        return 0
    return _SIM.get((aa1, aa2), 0)


def evolutionary_trace(
    msa: MSA,
    graph: nx.DiGraph,
    groups: dict[Clade, list[Clade]],
    root: Clade,
    leaf: Clade,
    d0: float = 0.05,
    progress: Callable[[int, int], None] | None = None,
) -> dict[int, float]:
    """Compute ET importance scores for every alignment column relative to *leaf*.

    Returns a dict mapping column index → score (lower = more important).
    """
    path: list[Clade] = nx.shortest_path(graph, root, leaf)
    leaf_idx = msa.sequence_indices[leaf.label]

    consensus_cache: dict[Clade, str] = {
        g: g.consensus
        for node in path if node.children
        for g in groups[node]
    }

    # Weights between group pairs are column-independent: precompute once.
    weighted_groups: list[tuple[Clade, float]] = []
    for node in path:
        if not node.children:
            continue
        host = next(
            (g for g in groups[node] if leaf_idx in g.sequence_indices), None
        )
        if host is None:
            continue
        host_con = consensus_cache[host]
        for g in groups[node]:
            d = 1.0 - float(np.mean(
                [similarity(a, b) for a, b in zip(host_con, consensus_cache[g])]
            ))
            w = float(np.exp(-(d ** 2) / (d0 ** 2)))
            weighted_groups.append((g, w))

    ranking: dict[int, float] = {}
    for col in range(msa.length):
        if progress:
            progress(col + 1, msa.length)
        ranking[col] = 1.0 + sum(w * g.entropy(col) for g, w in weighted_groups)
    return ranking
