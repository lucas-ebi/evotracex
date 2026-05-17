from __future__ import annotations

import numpy as np
import networkx as nx

from evotracex.constants import BLOSUM62, BLOSUM62_THRESHOLD
from evotracex.alignment import MSA
from evotracex.subsets import Clade


def similarity(aa1: str, aa2: str) -> int:
    """Binary similarity: 1 if identical or BLOSUM62 score >= mean diagonal, else 0."""
    if aa1 == aa2:
        return 1
    if "-" in (aa1, aa2):
        return 0
    pair = (aa1, aa2) if (aa1, aa2) in BLOSUM62 else (aa2, aa1)
    return int(BLOSUM62[pair] >= BLOSUM62_THRESHOLD)


def evolutionary_trace(
    msa: MSA,
    graph: nx.DiGraph,
    groups: dict[Clade, list[Clade]],
    root: Clade,
    leaf: Clade,
    d0: float = 0.05,
) -> dict[int, float]:
    """Compute ET importance scores for every alignment column relative to *leaf*.

    Returns a dict mapping column index → score (lower = more important).
    """
    path: list[Clade] = nx.shortest_path(graph, root, leaf)
    leaf_idx = msa.sequence_indices[leaf.label]

    # Pre-compute consensus strings for all groups on this path once.
    # In the original code this was done inside the column loop — O(L × nodes × groups).
    consensus_cache: dict[Clade, str] = {
        g: g.consensus
        for node in path if node.children
        for g in groups[node]
    }

    ranking: dict[int, float] = {}
    for col in range(msa.length):
        score = 0.0
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
                sim_vals = [
                    similarity(a, b)
                    for a, b in zip(host_con, consensus_cache[g])
                ]
                d = 1.0 - float(np.mean(sim_vals))
                w = float(np.exp(-(d ** 2) / (d0 ** 2)))
                score += w * g.entropy(col)
        ranking[col] = 1.0 + score
    return ranking
