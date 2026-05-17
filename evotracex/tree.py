from __future__ import annotations

import networkx as nx
from Bio import Phylo

from evotracex.alignment import MSA
from evotracex.subsets import Clade


def build_clade_network(
    msa: MSA, tree_file: str
) -> tuple[nx.DiGraph, Clade, list[Clade]]:
    """Parse a Newick tree and return (graph, root, leaves) as Clade objects.

    Raises ValueError if the tree leaf names do not exactly match the MSA headers.
    """
    try:
        tree = Phylo.read(tree_file, "newick")
    except FileNotFoundError:
        raise SystemExit(f"Tree file not found: '{tree_file}'")

    tree.ladderize()
    bio_clades = list(tree.find_clades(order="level"))

    # Validate names before building any objects
    tree_leaves = {c.name for c in bio_clades if c.is_terminal()}
    msa_names = set(msa.sequence_indices)
    if tree_leaves != msa_names:
        parts: list[str] = []
        missing = sorted(msa_names - tree_leaves)
        extra = sorted(tree_leaves - msa_names)
        if missing:
            parts.append(f"in MSA but not in tree: {missing}")
        if extra:
            parts.append(f"in tree but not in MSA: {extra}")
        raise ValueError("Tree/MSA name mismatch — " + "; ".join(parts))

    # Build Clade objects
    clade_map: dict = {}
    leaves: list[Clade] = []

    for i, bio_clade in enumerate(bio_clades):
        is_leaf = bio_clade.is_terminal()
        name = bio_clade.name if is_leaf else f"N{i}"
        terminals = [t.name for t in bio_clade.get_terminals()]
        seq_indices = [msa.sequence_indices[t] for t in terminals]
        clade_obj = Clade(msa, seq_indices, bio_clade.branch_length, name)
        clade_map[bio_clade] = clade_obj
        if is_leaf:
            leaves.append(clade_obj)

    # Connect in a directed graph (parent → child)
    graph: nx.DiGraph = nx.DiGraph()
    for bio_clade in bio_clades:
        parent = clade_map[bio_clade]
        for child_bio in bio_clade.clades:
            graph.add_edge(parent, clade_map[child_bio])

    # Populate Clade.children so callers can check node type without the graph
    for node in graph.nodes():
        node.children = list(graph.successors(node))

    root = next(n for n in graph.nodes() if graph.in_degree(n) == 0)
    return graph, root, leaves


def build_groups(
    graph: nx.DiGraph, root: Clade
) -> dict[Clade, list[Clade]]:
    """Determine the set of active subfamilies for each internal node.

    For a node at depth *d*, the subfamilies are the internal nodes that are
    roots (in-degree 0) of the sub-graph obtained by removing all nodes at
    depth < *d*.  This replicates the original progressive-pruning algorithm
    while computing depths in O(V) instead of O(V²).
    """
    depths = _compute_depths(graph, root)
    internal: set[Clade] = {n for n in graph.nodes() if n.children}

    aux = graph.copy()
    groups: dict[Clade, list[Clade]] = {}

    for node, _ in sorted(depths.items(), key=lambda x: x[1]):
        d = depths[node]
        for n in list(aux.nodes()):
            if depths[n] < d:
                aux.remove_node(n)
        if node in internal:
            groups[node] = [
                n for n in aux.nodes()
                if aux.in_degree(n) == 0 and n in internal
            ]

    return groups


def _compute_depths(graph: nx.DiGraph, root: Clade) -> dict[Clade, float]:
    """Accumulate branch-length depths via a single DFS from root — O(V)."""
    depths: dict[Clade, float] = {root: 0.0}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in graph.successors(node):
            depths[child] = depths[node] + (child.branch_length or 0.0)
            stack.append(child)
    return depths
