from __future__ import annotations

import numpy as np
from Bio import SeqIO

from evotracex.constants import ALPHABET, STEREOCHEMISTRY
from evotracex.subsets import Residue


class MSA:
    """Multiple Sequence Alignment with Henikoff weights and residue collection."""

    def __init__(self, msa_file: str, plus_aa: bool = False) -> None:
        self.plus_aa = plus_aa
        records = self._parse(msa_file)
        self.headers, self.sequences = self._read(records)
        self.size, self.length = self.sequences.shape
        self.weights: np.ndarray = self._henikoff_weights()
        self.sequence_indices: dict[str, int] = {
            h: i for i, h in enumerate(self.headers)
        }
        self.collection: dict[int, list[Residue]] = self._build_collection()
        self.gap_content: list[Residue] = self._count_gaps()

    def _parse(self, msa_file: str) -> list:
        try:
            return list(SeqIO.parse(msa_file, "fasta"))
        except FileNotFoundError:
            raise SystemExit(f"File not found: '{msa_file}'")

    def _read(self, records: list) -> tuple[list[str], np.ndarray]:
        headers = [r.id for r in records]
        sequences = np.array([list(str(r.seq)) for r in records])
        return headers, sequences

    def _henikoff_weights(self) -> np.ndarray:
        """Vectorised Henikoff (1994) sequence weighting."""
        weights = np.zeros(self.size)
        for col in range(self.length):
            column = self.sequences[:, col]
            _, inverse, counts = np.unique(
                column, return_inverse=True, return_counts=True
            )
            k = len(counts)
            weights += 1.0 / (k * counts[inverse])
        return weights / self.length

    def _build_collection(self) -> dict[int, list[Residue]]:
        collection: dict[int, list[Residue]] = {}

        for col in range(self.length):
            collection[col] = []
            for aa in ALPHABET:
                indices = [n for n in range(self.size) if self.sequences[n][col] == aa]
                if indices:
                    collection[col].append(Residue(self, [aa], col, indices))

        if self.plus_aa:
            for col in range(self.length):
                col_chars = set(self.sequences[:, col])

                # Map each stereo category → (amino_acids_present, sequence_indices)
                temp: dict[str, tuple[list[str], list[int]]] = {}
                for name, members in STEREOCHEMISTRY.items():
                    if col_chars & set(members):
                        temp[name] = ([], [])

                for n in range(self.size):
                    aa = self.sequences[n][col]
                    for name in list(temp):
                        if aa in STEREOCHEMISTRY[name]:
                            aas, idxs = temp[name]
                            if aa not in aas:
                                aas.append(aa)
                            idxs.append(n)

                # Freeze to tuples so we can use them as dict keys
                frozen: dict[str, tuple[tuple[str, ...], tuple[int, ...]]] = {
                    name: (tuple(aas), tuple(idxs))
                    for name, (aas, idxs) in temp.items()
                }

                # Group stereochemistry names by identical (aa, idx) fingerprint
                by_fingerprint: dict[
                    tuple[tuple[str, ...], tuple[int, ...]], list[str]
                ] = {}
                for name, fp in frozen.items():
                    by_fingerprint.setdefault(fp, []).append(name)

                for (aa_tuple, idx_tuple), names in by_fingerprint.items():
                    if len(aa_tuple) > 1:
                        # Normalise "Similar (X or Y)" → "Similar"
                        norm = [
                            "Similar" if "Similar" in n else n for n in names
                        ]
                        collection[col].append(
                            Residue(self, list(aa_tuple), col, list(idx_tuple), norm)
                        )

        return collection

    def _count_gaps(self) -> list[Residue]:
        result: list[Residue] = []
        for col in range(self.length):
            indices = [n for n in range(self.size) if self.sequences[n][col] == "-"]
            result.append(
                Residue(self, ["-"], col, indices, label=f"-{col + 1}")
            )
        return result
