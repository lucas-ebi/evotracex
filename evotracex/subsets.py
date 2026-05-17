from __future__ import annotations

import math
from collections.abc import Iterable
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from Bio.SeqUtils import seq3

if TYPE_CHECKING:
    from evotracex.alignment import MSA


class Subset:
    """A subset of sequences from an MSA, identified by sequence indices."""

    def __init__(
        self,
        msa: MSA,
        sequence_indices: Iterable[int],
        label: str | None = None,
    ) -> None:
        self.msa = msa
        self.sequence_indices: frozenset[int] = frozenset(sequence_indices)
        self.label = label

    def __repr__(self) -> str:
        return self.label or ""

    @cached_property
    def weight(self) -> float:
        """Sum of Henikoff weights for sequences in this subset."""
        return float(sum(self.msa.weights[i] for i in self.sequence_indices))

    def conditional_weight(self, given: Subset) -> float:
        """Fraction of *given*'s weight that overlaps with this subset."""
        shared = self.sequence_indices & given.sequence_indices
        return float(sum(self.msa.weights[i] for i in shared)) / given.weight

    @cached_property
    def consensus(self) -> str:
        """Majority amino acid at each column, '-' where the subset has no coverage."""
        result: list[str] = []
        for col in range(self.msa.length):
            candidates = [
                r for r in self.msa.collection[col]
                if len(r.amino_acids) == 1
                and not r.sequence_indices.isdisjoint(self.sequence_indices)
            ]
            if candidates:
                best = max(candidates, key=lambda r: r.conditional_weight(self))
                result.append(best.amino_acids[0])
            else:
                result.append("-")
        return "".join(result)

    def entropy(self, position: int) -> float:
        """Gap-corrected Shannon entropy at *position* for sequences in this subset.

        When msa.plus_aa is True, a greedy set cover selects the minimal set of
        non-overlapping stereochemical groups that spans all amino acids present,
        matching the original X-ET alphabet-expansion behaviour.
        """
        K = [
            r for r in self.msa.collection[position]
            if not r.sequence_indices.isdisjoint(self.sequence_indices)
        ]
        if not K:
            return math.log(20)

        if self.msa.plus_aa:
            content = {
                self.msa.sequences[i][position] for i in self.sequence_indices
            }
            K = sorted(K, key=lambda r: r.conditional_weight(self), reverse=True)
            selected: list[Residue] = []
            covered: set[str] = set()
            for r in K:
                if not selected:
                    selected.append(r)
                    covered.update(r.amino_acids)
                elif covered == content:
                    break
                elif set(r.amino_acids).isdisjoint(covered):
                    selected.append(r)
                    covered.update(r.amino_acids)
            K = selected

        g = self.msa.gap_content[position].conditional_weight(self)
        non_gap = 1.0 - g
        S = (non_gap * np.log(non_gap / 20) if non_gap > 0 else 0.0) + np.log(20)
        for r in K:
            p = r.conditional_weight(self)
            if p > 0:
                S -= p * np.log(p)
        return float(S)


class Residue(Subset):
    """A subset defined by having a specific amino acid (or group) at a position."""

    def __init__(
        self,
        msa: MSA,
        amino_acids: list[str],
        position: int,
        sequence_indices: Iterable[int],
        stereochemistry: list[str] | None = None,
        label: str | None = None,
    ) -> None:
        super().__init__(msa, sequence_indices, label)
        self.stereochemistry = stereochemistry
        self.amino_acids = list(amino_acids)
        self.position = position
        self.label = self._make_label(label)

    def _make_label(self, provided: str | None) -> str:
        if provided is not None:
            return provided
        if len(self.amino_acids) == 1:
            return seq3(self.amino_acids[0]) + str(self.position + 1)
        self.amino_acids = sorted(self.amino_acids)
        names = [
            "Similar" if "Similar" in s else s
            for s in (self.stereochemistry or [])
        ]
        names = sorted(set(names))
        feature = ", ".join(names) if len(names) > 1 else (names[0] if names else "")
        middle = ", ".join(seq3(a) for a in self.amino_acids[:-1])
        return f"{feature} ({middle} or {seq3(self.amino_acids[-1])}) {self.position + 1}"


class Clade(Subset):
    """A subset corresponding to a clade (node) in the phylogenetic tree."""

    def __init__(
        self,
        msa: MSA,
        sequence_indices: Iterable[int],
        branch_length: float | None = None,
        label: str | None = None,
    ) -> None:
        super().__init__(msa, sequence_indices, label)
        self.branch_length = branch_length
        self.children: list[Clade] = []
