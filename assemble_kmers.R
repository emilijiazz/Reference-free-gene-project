# Function to assemble k-mers into contigs
assemble_kmers <- function(kmers, min_overlap = 8) {
  library(igraph)

  # Create adjacency list for the de Bruijn graph
  edges <- data.frame(source = character(0), target = character(0), stringsAsFactors = FALSE)
  for (kmer in kmers) {
    for (other_kmer in kmers) {
      # Skip self-loops
      if (kmer == other_kmer) {
        next
      }
      # Find overlap between current k-mer and another k-mer
      for (overlap_size in seq(nchar(kmer) - 1, min_overlap, by = -1)) {
        if (substr(kmer, nchar(kmer) - overlap_size + 1, nchar(kmer)) == substr(other_kmer, 1, overlap_size)) {
          edges <- rbind(edges, data.frame(source = kmer, target = other_kmer, stringsAsFactors = FALSE))
          break
        }
      }
    }
  }

  # Handle case with no edges (no overlaps found)
  if (nrow(edges) == 0) {
    return(kmers)
  }

  # Create de Bruijn graph
  g <- graph_from_data_frame(edges, directed = TRUE)

  # Traverse the graph to find contigs
  contigs <- list()
  visited <- rep(FALSE, vcount(g))
  names(visited) <- V(g)$name

  assemble_contig <- function(start_node) {
    contig <- start_node
    path_kmers <- c(start_node)
    current_node <- start_node

    while (TRUE) {
      neighbors <- neighbors(g, current_node, mode = "out")
      if (length(neighbors) == 0) {
        break
      }
      next_node <- V(g)$name[neighbors[1]]  # Take the first neighbor
      
      # Check if we've already visited this node (avoid cycles)
      if (visited[next_node]) {
        break
      }
      
      # Find the overlap size between current_node and next_node
      overlap_found <- FALSE
      for (overlap_size in seq(nchar(current_node) - 1, min_overlap, by = -1)) {
        if (substr(current_node, nchar(current_node) - overlap_size + 1, nchar(current_node)) == 
            substr(next_node, 1, overlap_size)) {
          # Append the non-overlapping part of next_node to contig
          contig <- paste0(contig, substr(next_node, overlap_size + 1, nchar(next_node)))
          overlap_found <- TRUE
          break
        }
      }
      
      if (!overlap_found) {
        break
      }
      
      path_kmers <- c(path_kmers, next_node)
      visited[next_node] <<- TRUE
      current_node <- next_node
    }

    return(contig)
  }

  # Find all contigs by traversing from each unvisited node
  for (node_name in V(g)$name) {
    if (!visited[node_name]) {
      visited[node_name] <- TRUE
      contig <- assemble_contig(node_name)
      contigs <- c(contigs, contig)
    }
  }

  return(unlist(contigs))
}
