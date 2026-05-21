# EvoTraceX

EvoTraceX implements **Evolutionary Trace (ET) with optional alphabet expansion** — a method for identifying functionally important residue positions in a protein family.

Given a multiple sequence alignment (MSA) and a phylogenetic tree, it ranks each alignment column by how well its variation mirrors the tree's branching pattern. Residues that are conserved within subfamilies but differ between them are highlighted as likely determinants of functional specificity. Sequence redundancy bias is corrected by Henikoff & Henikoff (1994) positional weighting; column conservation within each subfamily is quantified by a gap-corrected Shannon entropy.

The alphabet expansion (`--expand`) groups amino acids by stereochemical properties (Taylor 1986; e.g., Aliphatic, Aromatic, Hydrophilic), boosting sensitivity at positions where physicochemical property is conserved even when the exact residue is not. In terms of the entropy formula, enabling `--expand` replaces individual residues with stereochemical groups as the symbols over which $\tilde{p}_i$ is defined.

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
| `-x/--expand` | off | Enable alphabet expansion (stereochemical groupings) |
| `--d0 FLOAT` | `0.05` | Gaussian bandwidth for subfamily weighting |
| `--ref NAME` | first sequence | Reference sequence for analysis and position numbering (must match an MSA header exactly). Restricts output to this sequence; alignment columns where it has a gap are omitted. |
| `--marginal` | off | Detect marginally conserved positions (see below) |
| `--fdr FLOAT` | `0.05` | Benjamini-Hochberg FDR threshold for marginal calls |

### Standard ET ranking

Ranks every alignment column by functional importance relative to the reference sequence.
Positions are numbered by residues in the reference (gaps excluded).

```bash
evotracex alignment.fasta tree.nwk --ref seqA --d0 0.1 --out results
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

### Expanded ET ranking

Same as above but residue conservation is measured at the level of stereochemical groups (Taylor 1986), boosting sensitivity at physicochemically constrained positions.

```bash
evotracex alignment.fasta tree.nwk --ref seqA --expand --d0 0.1 --out results
# → results_et_ranking.tsv
```

Output format is identical; scores will differ at positions where group-level conservation is stronger than residue-level conservation.

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
| `delta` | `arctan((score_standard − score_expanded) / (score_standard + score_expanded))` — angle below the score diagonal in radians; positive means expansion helped, scale-invariant |
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

The method is based on the Evolutionary Trace approach (Lichtarge et al. 1996), which uses a phylogenetic tree to partition sequences into nested subfamilies and identifies positions whose conservation pattern mirrors the tree topology. The alphabet expansion (Taylor 1986) extends its functionality by replacing single-residue conservation with stereochemical-group conservation to capture marginally conserved sites. The marginal conservation detection introduced here formalises this comparison using a Z-score framework with Benjamini-Hochberg multiple-testing correction.

**Sequence weighting (Henikoff & Henikoff 1994).** To reduce bias from over-represented sequence clusters, EvoTraceX applies positional weights. For sequence $i$ and alignment column $j$, let $k(j)$ be the number of distinct residue types in column $j$ and $c(i,j)$ the count of the residue carried by sequence $i$ there. The weight is:

```math
w_i = \frac{1}{L} \sum_j \frac{1}{k(j)\cdot c(i,j)},
```

where $L$ is the alignment length and the outer $1/L$ normalises so that $\sum_i w_i = 1$. Rare residues at a column (small $c$, large $k$) receive a higher contribution, upweighting divergent sequences and downweighting redundant clusters without requiring an explicit identity threshold.

**Gap-corrected Shannon entropy.** Without a gap correction, a column where most sequences carry a gap can appear spuriously conserved: the few observed residues may happen to agree, yielding a low raw entropy that ranks the position as functionally important. To prevent this and ensure that gap-rich columns are penalised as variable rather than conserved, EvoTraceX adapts the Shannon entropy to account for gap content.

Let $g$ be the total Henikoff weight of gap-bearing sequences at a position, and $p_i$ the raw weight fraction of residue group $i$ relative to the full subset (so $\sum_i p_i = 1 - g$). Starting from Shannon's definition, the entropy of the actual column is:

```math
H^\text{actual} = -\sum_i p_i \ln p_i.
```

A **reference column** is then constructed with the same gap fraction $g$ but with the non-gap weight spread uniformly over all 20 amino acids ($p_i^\text{ref} = (1-g)/20$). Its entropy follows directly from Shannon's definition:

```math
H^\text{ref} = -\sum_{i=1}^{20} \frac{1-g}{20}\ln\frac{1-g}{20} = -(1-g)\ln\frac{1-g}{20}.
```

The gap-corrected score is defined by shifting $H^\text{actual}$ so that a column matching the reference maps to exactly $\ln 20$ — the maximum uncertainty — regardless of gap content:

```math
S = H^\text{actual} + \bigl(\ln 20 - H^\text{ref}\bigr)
  = -\sum_i p_i \ln p_i + (1-g)\ln\frac{1-g}{20} + \ln 20.
```

Substituting $p_i = (1-g)\,\tilde{p}_i$, where $\tilde{p}_i = p_i/(1-g)$ is the residue distribution conditioned on non-gap sequences ($\sum_i \tilde{p}_i = 1$), gives:

```math
-\sum_i p_i \ln p_i = -(1-g)\ln(1-g) + (1-g)\,H(\tilde{p}_i),
```

and substituting $H^\text{ref}$ gives $(1-g)\ln\frac{1-g}{20} = (1-g)\ln(1-g) - (1-g)\ln 20$, so the $(1-g)\ln(1-g)$ terms cancel and the expression simplifies to:

```math
S = g\ln 20 + (1-g)\,H(\tilde{p}_i),
\qquad H(\tilde{p}_i) = -\sum_i \tilde{p}_i \ln \tilde{p}_i.
```

Fully conserved, gap-free columns give $S = 0$; a fully gapped column gives $S = \ln 20$, the entropy of a uniform distribution over all 20 amino acids.

### References

- Lichtarge O, Bourne HR, Cohen FE. (1996). An evolutionary trace method defines binding surfaces common to protein families. *J Mol Biol*, 257(2), 342–358. <https://doi.org/10.1006/jmbi.1996.0167>
- Henikoff S, Henikoff JG. (1994). Position-based sequence weights. *J Mol Biol*, 243(4), 574–578. <https://doi.org/10.1016/0022-2836(94)90032-9>
- Taylor WR. (1986). The classification of amino acid conservation. *J Theor Biol*, 119(2), 205–218. <https://doi.org/10.1016/S0022-5193(86)80075-3>

## License

GNU General Public License v3. See [LICENSE](LICENSE).
