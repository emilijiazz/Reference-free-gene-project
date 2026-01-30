#!/usr/bin/env Rscript
# Advanced example demonstrating the assemble_kmers function

source("assemble_kmers.R")

cat("=== Advanced K-mer Assembly Examples ===\n\n")

# Example 1: Simulating genomic DNA fragment assembly
cat("Example 1: Genomic DNA Fragment Assembly\n")
cat("-------------------------------------------\n")
dna_kmers <- c(
  "ATGCGTA",
  "CGTAGCA",
  "AGCATGC",
  "TGCGCTA"
)
cat("Input DNA k-mers:\n")
print(dna_kmers)
cat("\nAssembling with min_overlap=4:\n")
result1 <- assemble_kmers(dna_kmers, min_overlap = 4)
cat("Result:", result1, "\n\n")

# Example 2: Different overlap thresholds
cat("Example 2: Effect of Different min_overlap Values\n")
cat("---------------------------------------------------\n")
test_kmers <- c("ABCDEF", "CDEFGH", "EFGHIJ")
cat("Input k-mers:", paste(test_kmers, collapse = ", "), "\n\n")

for (overlap in c(2, 3, 4, 5)) {
  result <- assemble_kmers(test_kmers, min_overlap = overlap)
  cat(sprintf("min_overlap=%d: %s\n", overlap, paste(result, collapse = ", ")))
}
cat("\n")

# Example 3: Multiple independent contigs
cat("Example 3: Multiple Independent Contigs\n")
cat("----------------------------------------\n")
multi_kmers <- c(
  "AAATTT",  # First contig start
  "ATTTGG",  # First contig middle
  "TTTGGG",  # First contig end
  "CCCGGG",  # Second contig start
  "CGGGAA",  # Second contig end
  "XXXXXX"   # Isolated k-mer
)
cat("Input k-mers:\n")
print(multi_kmers)
cat("\nAssembling with min_overlap=3:\n")
result3 <- assemble_kmers(multi_kmers, min_overlap = 3)
cat("Result contigs:\n")
for (i in seq_along(result3)) {
  cat(sprintf("Contig %d: %s\n", i, result3[i]))
}
cat("\n")

# Example 4: Short k-mers with varying overlap
cat("Example 4: Short K-mers Assembly\n")
cat("---------------------------------\n")
short_kmers <- c("ATCG", "TCGA", "CGAT", "GATC")
cat("Input k-mers:", paste(short_kmers, collapse = ", "), "\n")
result4 <- assemble_kmers(short_kmers, min_overlap = 3)
cat("Result:", paste(result4, collapse = ", "), "\n\n")

# Example 5: Handling no overlap case
cat("Example 5: No Overlap Detection\n")
cat("--------------------------------\n")
no_overlap_kmers <- c("AAAA", "TTTT", "GGGG", "CCCC")
cat("Input k-mers:", paste(no_overlap_kmers, collapse = ", "), "\n")
result5 <- assemble_kmers(no_overlap_kmers, min_overlap = 3)
cat("Result (should return original k-mers):\n")
print(result5)
cat("\n")

cat("=== All examples completed successfully! ===\n")
