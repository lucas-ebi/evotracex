# EvoTraceX

EvoTraceX implements **Evolutionary Trace (ET) with optional alphabet expansion** — a method for identifying functionally important residue positions in a protein family.

Given a multiple sequence alignment (MSA) and a phylogenetic tree, it ranks each alignment column by how well its variation mirrors the tree's branching pattern. Residues that are conserved within subfamilies but differ between them are highlighted as likely determinants of functional specificity.

The alphabet expansion (`--plus-aa`) groups amino acids by stereochemical properties (e.g., Aliphatic, Aromatic, Hydrophilic), boosting sensitivity at positions where chemical character is conserved even when the exact residue is not.

## Installation

Requires Python ≥ 3.10.

```bash
pip install -e ".[dev]"   # editable install with test dependencies
```

Dependencies: `biopython>=1.80`, `numpy>=1.24`, `networkx>=3.0`.

## Usage

```bash
evotracex <msa.fasta> <tree.nwk> [OPTIONS]
```

### Options

| Flag | Description |
| --- | --- |
| `-o/--out PREFIX` | Write results to `PREFIX_et_ranking.tsv` instead of stdout |
| `-x/--plus-aa` | Enable alphabet expansion (stereochemical groupings) |
| `--d0 FLOAT` | Gaussian bandwidth for subfamily weighting (default: `0.05`) |
| `--leaf NAME` | Analyse one named leaf only (default: all leaves) |

### Example

```bash
evotracex alignment.fasta tree.nwk --plus-aa --d0 0.1 --out results
# writes results_et_ranking.tsv
```

### Output format

Tab-separated, one block per leaf (preceded by a `# leaf_name` comment line):

```text
# seqA
rank    residue    position    score
1       Ala        12          1.000000
2       Leu        45          1.023411
...
```

Lower score = more important. Gap positions are omitted.

## Running the tests

```bash
pytest tests/ -v
```

## Scientific background

The method is based on the Evolutionary Trace approach (Lichtarge et al. 1996), which uses a phylogenetic tree to partition sequences into nested subfamilies and identifies positions whose conservation pattern mirrors the tree topology. The alphabet expansion is described in the X-ET extension (Carrijo de Oliveira 2018), which replaces single-residue conservation with stereochemical-group conservation to capture marginally conserved sites.

## License

GNU General Public License v3. See [LICENSE](LICENSE).
