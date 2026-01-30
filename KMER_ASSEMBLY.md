# K-mer Assembly Function

## Description

This R script provides a function `assemble_kmers()` that assembles k-mers into contigs using a graph-based approach. The function builds a de Bruijn graph from the given list of k-mers and assembles them into contigs while allowing for adjustable control over the `min_overlap` parameter.

## Function Signature

```R
assemble_kmers(kmers, min_overlap = 8)
```

### Parameters

- **kmers**: A vector of character strings representing k-mers to be assembled
- **min_overlap** (default = 8): The minimum number of overlapping characters required between k-mers to create an edge in the de Bruijn graph

### Returns

A vector of character strings representing the assembled contigs.

## How It Works

1. **Graph Construction**: Creates a de Bruijn graph where:
   - Each k-mer is a node
   - Edges are created between k-mers that have sufficient overlap (≥ min_overlap)
   - The overlap is determined by comparing the suffix of one k-mer with the prefix of another

2. **Graph Traversal**: Traverses the graph to assemble contigs:
   - Starts from each unvisited node
   - Follows edges to neighboring nodes
   - Merges k-mers by removing overlapping regions
   - Avoids cycles by tracking visited nodes

3. **Contig Assembly**: Returns all assembled contigs as a vector

## Usage Example

```R
source("assemble_kmers.R")

# Example 1: Simple linear assembly
kmers1 <- c("ATCG", "TCGA", "CGAT")
result1 <- assemble_kmers(kmers1, min_overlap = 3)
# Result: "ATCGAT"

# Example 2: Multiple overlapping k-mers
kmers2 <- c("ABCDE", "CDEFG", "EFGHI")
result2 <- assemble_kmers(kmers2, min_overlap = 3)
# Result: "ABCDEFGHI"

# Example 3: DNA sequence assembly
kmers3 <- c("ATGCGA", "GCGATT", "GATTCC", "ATTCCA")
result3 <- assemble_kmers(kmers3, min_overlap = 4)
# Result: "ATGCGATTCCA"
```

## Requirements

- R (version 4.0 or higher)
- igraph package

### Installing igraph

```R
install.packages("igraph")
```

Or on Ubuntu/Debian:
```bash
sudo apt install r-cran-igraph
```

## Testing

Run the test script to verify the function works correctly:

```bash
Rscript test_assemble_kmers.R
```

## Key Features

- **Flexible Overlap Control**: Adjustable `min_overlap` parameter allows fine-tuning of assembly stringency
- **Cycle Detection**: Avoids infinite loops by tracking visited nodes
- **Self-loop Prevention**: Automatically excludes edges from a k-mer to itself
- **Robust Error Handling**: Returns original k-mers if no overlaps are found
- **Greedy Assembly**: Uses a greedy approach for contig extension

## Limitations

- The function uses a greedy approach and takes the first neighbor when multiple neighbors exist
- Does not handle complex branching structures optimally
- Performance may degrade with very large k-mer sets due to nested loops in edge construction

## Code Fixes

The original incomplete code had the following issues that were fixed:

1. **Incomplete line**: The line `shared_overlap <- sum(substring(start_node, -min_overlap, nchar(start_value) returned0..` was incomplete and had multiple syntax errors
2. **Missing logic**: The contig assembly logic was incomplete
3. **Variable name errors**: References to undefined variables like `start_value`
4. **Missing cycle detection**: No mechanism to prevent infinite loops
5. **Incomplete traversal logic**: The graph traversal and contig extension logic was not implemented

All these issues have been resolved in the current implementation.
