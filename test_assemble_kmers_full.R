# test_assemble_kmers_full.R
# Robust testing for assemble_kmers function
# Covers: linear assembly, branchpoints (single/multi-cohort), tie-breaks, no-overlap, edge-cases, strict overlap, verbose/debug

source("assemble_kmers_finalversion.R")  # Load your function. Adjust filename if needed.

cat("\n====== TEST 1: Simple linear assembly ======\n")
kmers1 <- c("ATCG", "TCGA", "CGAT")
result1 <- assemble_kmers(kmers1, min_overlap = 3)
print_contig_info(result1)
cat("Expected: Single contig, all k-mers used in sequence.\n\n")

cat("====== TEST 2: Branchpoint (single cohort with enrichment/delta) ======\n")
kmers2 <- c("AAAAT", "AAATG", "AAACC", "AATCC")
kmer_count_df2 <- data.frame(
  kmer        = kmers2,
  enrichment  = c(10, 20, 5, 1),
  delta       = c(100, 200, 50, 10),
  stringsAsFactors = FALSE
)
result2 <- assemble_kmers(
  kmers = kmer_count_df2$kmer,
  min_overlap = 4,
  kmer_count_df = kmer_count_df2,
  enrichment_threshold = 1.5,
  delta_threshold = 100,
  tie_action = "stop"
)
print_contig_info(result2)
cat("Expected: Branch resolved by enrichment/delta threshold; selects AAATG.\n\n")

cat("====== TEST 3: Branchpoint (multi-cohort) ======\n")
kmers3 <- c("AAAAT", "AAATG", "AAACC", "AATCC")
kmer_count_df3 <- data.frame(
  kmer        = kmers3,
  enrichment_x= c(10, 20, 5, 1),
  enrichment_y= c(15, 30, 4, 1),
  delta_x     = c(100, 200, 50, 10),
  delta_y     = c(110, 250, 45, 10),
  stringsAsFactors = FALSE
)
result3 <- assemble_kmers(
  kmers = kmer_count_df3$kmer,
  min_overlap = 4,
  kmer_count_df = kmer_count_df3,
  enrichment_threshold = 1.5,
  delta_threshold = 100,
  tie_action = "stop"
)
print_contig_info(result3)
cat("Expected: Branch resolved using minimum enrichment/delta across cohorts; selects AAATG.\n\n")

cat("====== TEST 4: Branchpoint (tie-break, lexical order) ======\n")
kmers4 <- c("ATCGA", "CGAAC", "CGAAG")
kmer_count_df4 <- data.frame(
  kmer        = kmers4,
  enrichment  = c(5, 5, 5),
  delta       = c(100, 100, 100),
  stringsAsFactors = FALSE
)
result4 <- assemble_kmers(
  kmers = kmer_count_df4$kmer,
  min_overlap = 3,
  kmer_count_df = kmer_count_df4,
  enrichment_threshold = 1.0,
  delta_threshold = 0,
  tie_action = "pick_lexical"
)
print_contig_info(result4)
cat("Expected: Branchpoint with identical enrichment/delta; picks CGAAC (lexical sort).\n\n")

# Extra explicit lexical tiebreak test
cat("\n====== TEST 4B: Lexical tiebreak at branch with identical enrichment/delta ======\n")
kmers <- c("ATCGA", "CGAAG", "CGAAC")
kmer_count_df <- data.frame(
  kmer        = kmers,
  enrichment  = c(5, 5, 5),
  delta       = c(100, 100, 100),
  stringsAsFactors = FALSE
)
result <- assemble_kmers(
  kmers = kmer_count_df$kmer,
  min_overlap = 3,
  kmer_count_df = kmer_count_df,
  enrichment_threshold = 1.0,
  delta_threshold = 0,
  tie_action = "pick_lexical",
  verbose = TRUE
)
print_contig_info(result)
cat("Expected: Picks CGAAC (alphabetically before CGAAG) after ATCGA.\n\n")

cat("====== TEST 5: No overlap (isolated k-mers) ======\n")
kmers5 <- c("AAAA", "TTTT", "GGGG")
result5 <- assemble_kmers(kmers5, min_overlap = 3)
print_contig_info(result5)
cat("Expected: Each k-mer returned as an isolated contig.\n\n")

cat("====== TEST 6: Single cohort with varying thresholds ======\n")
kmers6 <- c("AAATT", "AATTC", "ATTCC")
kmer_count_df6 <- data.frame(
  kmer        = kmers6,
  enrichment  = c(2, 1.3, 1.1),
  delta       = c(400, 180, 90),
  stringsAsFactors = FALSE
)
result6 <- assemble_kmers(
  kmers = kmer_count_df6$kmer,
  min_overlap = 4,
  kmer_count_df = kmer_count_df6,
  enrichment_threshold = 1.25,
  delta_threshold = 200,
  tie_action = "stop"
)
print_contig_info(result6)
cat("Expected: Only k-mers meeting threshold are merged; others isolated.\n\n")

cat("====== TEST 7: Multi-cohort (edge-case: missing/extreme values) ======\n")
kmers7 <- c("AAGTC", "AGTCA", "GTCAA")
kmer_count_df7 <- data.frame(
  kmer        = kmers7,
  enrichment_x= c(10, NA, 1),
  enrichment_y= c(9, 2, NA),
  delta_x     = c(100, NA, 1),
  delta_y     = c(90, 2, NA),
  stringsAsFactors = FALSE
)
result7 <- assemble_kmers(
  kmers = kmer_count_df7$kmer,
  min_overlap = 4,
  kmer_count_df = kmer_count_df7,
  enrichment_threshold = 1.5,
  delta_threshold = 10,
  tie_action = "stop"
)
print_contig_info(result7)
cat("Expected: Handles missing/extreme values without crash; merges only valid k-mers.\n\n")

cat("====== TEST 8: Strict kminus1 (k-mer length - 1 overlap) ======\n")
kmers8 <- c("ATCGG", "TCGGA", "CGGAT")
result8 <- assemble_kmers(
  kmers = kmers8,
  strict_kminus1 = TRUE
)
print_contig_info(result8)
cat("Expected: min_overlap = k-1; maximal merging.\n\n")

cat("====== TEST 9: Verbose & branch debug ======\n")
kmers9 <- c("ATCGA", "CGAAC", "CGAAG")
kmer_count_df9 <- data.frame(
  kmer        = kmers9,
  enrichment  = c(5, 5, 5),
  delta       = c(100, 100, 100),
  stringsAsFactors = FALSE
)
result9 <- assemble_kmers(
  kmers = kmer_count_df9$kmer,
  min_overlap = 3,
  kmer_count_df = kmer_count_df9,
  enrichment_threshold = 1.0,
  delta_threshold = 0,
  tie_action = "pick_lexical",
  verbose = TRUE
)
print_contig_info(result9)
cat("Expected: Branchpoints and node selection printed to console.\n\n")

cat("====== TEST 10: Long linear contig (redundancy/full merge) ======\n")
kmers_long <- c(
  "AAAAG",     # 5-mer
  "AAAGG",
  "AAGGG",
  "AGGGG",
  "GGGGT",
  "GGGTG"
)
# All overlap on last 4 letters, so one long contig possible

result_long <- assemble_kmers(
  kmers = kmers_long,
  min_overlap = 4
)
print_contig_info(result_long)
cat("Expected: One long contig, all k-mers merged in sequence. No fragments.\n\n")
cat("\n====== ALL TESTS COMPLETE ======\n")