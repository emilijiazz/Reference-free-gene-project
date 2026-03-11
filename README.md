# Reference-free-gene project
# assemble_kmers:  k-mer Graph Assembly in R, built on Bruijn graph-based approach

## Overview

`assemble_kmers` is an R function for assembling sequence contigs from small sets of k-mers. It supports linear and branching paths, single or multi-cohort enrichment/delta filtering, customizable overlap requirements, tie-breaking rules, and robust handling of real-world dataset irregularities. The algorithm handles arbitrarily long chains of overlapping k-mers in a single run, without requiring secondary merging rounds.

---

## Features

- **Branching logic:**  
  Automatically resolves branchpoints in the k-mer graph, selecting paths based on enrichment and delta counts.
- **Customizable overlap:**  
  Default overlap is k-1, but you can specify any `min_overlap`.
- **Single & multi-cohort support:**  
  Dataframes with `enrichment`, `delta` (single cohort) or `enrichment_x/enrichment_y`, `delta_x/delta_y` (multi-cohort) are supported. Branch decisions use minimum values across cohorts.
- **Tie-breaker options:**  
  When enrichment and delta are indistinguishable, set `tie_action = "pick_lexical"` to resolve alphabetically.
- **One-pass robust merging:**  
  Long chains of overlapping k-mers are merged into a single contig in one step — no need for a second merge round.
- **Isolated k-mer handling:**  
  K-mers with no overlaps are returned as single contigs.
- **Edge-case robustness:**  
  Handles missing or extreme values in data without crashing.
- **Verbose output:**  
  Set `verbose = TRUE` to print branch and node selection decisions.

---

## Installation

Simply copy `assemble_kmers_finalversion.R` and `print_contig_info` function to your R project.  
Requires the `igraph` package.

---

## Function Usage

```r
source("assemble_kmers.R")  # Load the function

contigs <- assemble_kmers(
  kmers = kmer_count_df$kmer,       # vector of k-mers
  kmer_count_df = kmer_count_df,    # dataframe with enrichment/delta info
  enrichment_threshold = 1.25,      # user-specified thresholds
  delta_threshold = 200,
  tie_action = "stop",              # or "pick_lexical"
  verbose = TRUE                    # optional: see branching details
)

print_contig_info(contigs)          # neatly print all assembled contigs
```

- For single cohort, ensure `kmer_count_df` has `enrichment` and `delta` columns.
- For two cohorts, use `enrichment_x`, `enrichment_y`, `delta_x`, and `delta_y`.

---

## Example Tests

See `test_assemble_kmers_full.R` for a comprehensive suite covering:

- Simple linear assembly (all k-mers merged in sequence)
- Branchpoints resolved by enrichment/delta (single and multi-cohort)
- Tie-breaking lexically when enrichment/delta are equal
- Isolated k-mers (no overlaps)
- Handling of missing/extreme values
- Maximal merging with strict kminus1 setting
- Verbose/debug output for transparency
- Long redundant contig test (multi-step chain merge)

**Example: Long linear contig**
```r
kmers_long <- c("AAAAG", "AAAGG", "AAGGG", "AGGGG", "GGGGT", "GGGTG")
result_long <- assemble_kmers(kmers = kmers_long, min_overlap = 4)
print_contig_info(result_long)
```
Output: a single contig merging all k-mers, demonstrating redundancy and robustness.

---

## Scientific Summary

- Handles branching
- Selects optimal paths by enrichment/delta
- Supports single and multi-cohort datasets
- Tie-breaks resolve ambiguities
- Merges any chain of overlapping k-mers in one step (no extra merging needed)
- Robust to edge-cases and real data irregularities

---



## Citation

If you use or adapt this function, please cite this repository and acknowledge the original algorithm. For questions or improvements, open an issue or pull request.

---

## License

MIT License.
