# EvoTraceX

EvoTraceX implements **Evolutionary Trace (ET) with optional alphabet expansion** — a method for identifying functionally important residue positions in a protein family.

Given a multiple sequence alignment (MSA) and a phylogenetic tree, it ranks each alignment column by how well its variation mirrors the tree's branching pattern. Residues that are conserved within subfamilies but differ between them are highlighted as likely determinants of functional specificity.

The alphabet expansion (`--plus-aa`) groups amino acids by stereochemical properties (e.g., Aliphatic, Aromatic, Hydrophilic), boosting sensitivity at positions where chemical character is conserved even when the exact residue is not.

The marginal conservation mode (`--marginal`) automatically compares standard and expanded rankings to pinpoint positions that only emerge as conserved under the expanded alphabet — the defining signature of marginal conservation.

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

| Flag | Default | Description |
| --- | --- | --- |
| `-o/--out PREFIX` | stdout | Write results to `PREFIX_et_ranking.tsv` (or `PREFIX_marginal.tsv`) |
| `-x/--plus-aa` | off | Enable alphabet expansion (stereochemical groupings) |
| `--d0 FLOAT` | `0.05` | Gaussian bandwidth for subfamily weighting |
| `--ref NAME` | first sequence | Reference sequence for analysis and position numbering (must match an MSA header exactly). Restricts output to this sequence; alignment columns where it has a gap are omitted. |
| `--marginal` | off | Detect marginally conserved positions (see below) |
| `--fdr FLOAT` | `0.05` | Benjamini-Hochberg FDR threshold for marginal calls |

### Standard ET ranking

Ranks every alignment column by functional importance relative to the reference sequence.
Positions are numbered by residues in the reference (gaps excluded).

```bash
evotracex alignment.fasta tree.nwk --ref seqA --plus-aa --d0 0.1 --out results
# → results_et_ranking.tsv
```

Output — lower score = more important:

```text
# seqA
rank    residue    position    score
1       Ala        12          1.000000
2       Leu        45          1.023411
...
```

### Marginal conservation detection

Runs ET twice (with and without alphabet expansion) and compares the rankings
position by position. Positions whose score drops significantly when the expanded
alphabet is applied are flagged as marginally conserved.

```bash
evotracex alignment.fasta tree.nwk --ref seqA --marginal --fdr 0.05 --out results
# → results_marginal.tsv
```

Output:

```text
# seqA
position    residue    score_standard    score_expanded    delta    z_score    p_value    adj_p_value    marginal
23          Val        1.412600          1.031200          0.381400 2.3112     0.0104     0.0312         TRUE
45          Leu        1.389100          1.201300          0.187800 1.0431     0.1484     0.3710         FALSE
...
```

| Column | Meaning |
| --- | --- |
| `delta` | `score_standard − score_expanded` — positive means expansion helped |
| `z_score` | Standardised delta across all positions for this leaf |
| `p_value` | One-tailed p-value (upper tail) from the Z-score |
| `adj_p_value` | Benjamini-Hochberg FDR corrected p-value |
| `marginal` | `TRUE` if `delta > 0` and `adj_p_value < --fdr` |

> **Note on alignment size.** BH correction is conservative with very few columns (<~30). On small datasets, ranked z-scores are a reliable guide even when no position crosses the `marginal = TRUE` threshold.

## Running the tests

```bash
pytest tests/ -v
```

## Scientific background

The method is based on the Evolutionary Trace approach (Lichtarge et al. 1996), which uses a phylogenetic tree to partition sequences into nested subfamilies and identifies positions whose conservation pattern mirrors the tree topology. The alphabet expansion is described in the X-ET extension (Carrijo de Oliveira 2018), which replaces single-residue conservation with stereochemical-group conservation to capture marginally conserved sites. The marginal conservation detection introduced here formalises this comparison using a Z-score framework with Benjamini-Hochberg multiple-testing correction.

## License

GNU General Public License v3. See [LICENSE](LICENSE).
