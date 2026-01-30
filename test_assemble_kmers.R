# Test script for assemble_kmers function
source("assemble_kmers.R")

# Test 1: Simple linear assembly
cat("Test 1: Simple linear assembly\n")
kmers1 <- c("ATCG", "TCGA", "CGAT")
result1 <- assemble_kmers(kmers1, min_overlap = 3)
cat("Input k-mers:", paste(kmers1, collapse = ", "), "\n")
cat("Result:", paste(result1, collapse = ", "), "\n")
cat("Expected: ATCGAT (or similar based on overlap)\n\n")

# Test 2: More complex example with multiple k-mers
cat("Test 2: Multiple overlapping k-mers\n")
kmers2 <- c("ABCDE", "CDEFG", "EFGHI")
result2 <- assemble_kmers(kmers2, min_overlap = 3)
cat("Input k-mers:", paste(kmers2, collapse = ", "), "\n")
cat("Result:", paste(result2, collapse = ", "), "\n")
cat("Expected: ABCDEFGHI\n\n")

# Test 3: No overlap case
cat("Test 3: No overlap (min_overlap too high)\n")
kmers3 <- c("AAAA", "TTTT", "GGGG")
result3 <- assemble_kmers(kmers3, min_overlap = 3)
cat("Input k-mers:", paste(kmers3, collapse = ", "), "\n")
cat("Result:", paste(result3, collapse = ", "), "\n")
cat("Expected: Original k-mers returned as no overlap\n\n")

# Test 4: DNA-like sequence
cat("Test 4: DNA-like sequence\n")
kmers4 <- c("ATGCGA", "GCGATT", "GATTCC", "ATTCCA")
result4 <- assemble_kmers(kmers4, min_overlap = 4)
cat("Input k-mers:", paste(kmers4, collapse = ", "), "\n")
cat("Result:", paste(result4, collapse = ", "), "\n\n")

cat("All tests completed!\n")
