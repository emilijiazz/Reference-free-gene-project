assemble_kmers <- function(
    kmers,
    min_overlap = NULL,
    kmer_count_df = NULL,
    strict_kminus1 = FALSE,
    enrichment_threshold = 1.15,
    delta_threshold = 10000,
    tie_action = c("stop", "pick_lexical"),
    verbose = FALSE
) {
  suppressPackageStartupMessages(library(igraph))
  tie_action <- match.arg(tie_action)
  kmers <- unique(kmers)
  if (is.null(min_overlap)) min_overlap <- 8
  
  if (strict_kminus1) {
    lens <- unique(nchar(kmers))
    if (length(lens) == 1) min_overlap <- lens[1] - 1
    else warning("strict_kminus1=TRUE but varying k-mer length")
  }
  
  edge_list <- vector("list", length(kmers) * length(kmers))
  idx <- 1L
  for (kmer in kmers) {
    for (other_kmer in kmers) {
      if (kmer == other_kmer) next
      max_ol <- min(nchar(kmer), nchar(other_kmer)) - 1
      for (overlap_size in seq(max_ol, min_overlap, by = -1)) {
        if (substr(kmer, nchar(kmer) - overlap_size + 1, nchar(kmer)) ==
            substr(other_kmer, 1, overlap_size)) {
          edge_list[[idx]] <- list(source = kmer, target = other_kmer, overlap = overlap_size)
          idx <- idx + 1L
          break
        }
      }
    }
  }
  edge_list <- edge_list[seq_len(idx - 1L)]
  edges <- if (length(edge_list)) {
    data.frame(
      source  = vapply(edge_list, `[[`, character(1), "source"),
      target  = vapply(edge_list, `[[`, character(1), "target"),
      overlap = as.integer(vapply(edge_list, `[[`, numeric(1), "overlap")),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(source = character(0), target = character(0), overlap = integer(0), stringsAsFactors = FALSE)
  }
  
  contigs_info <- list()
  if (nrow(edges) == 0) {
    for (i in seq_along(kmers)) {
      contigs_info[[i]] <- list(contig = kmers[i], kmers_used = kmers[i], num_kmers = 1, is_isolated = TRUE)
    }
    return(contigs_info)
  }
  
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = kmers)
  E(g)$overlap <- edges$overlap
  
  kmers_in_graph <- V(g)$name
  isolated_kmers <- setdiff(kmers, kmers_in_graph)
  visited <- rep(FALSE, vcount(g)); names(visited) <- V(g)$name
  
  pick_next_node <- function(current_node) {
    neigh <- neighbors(g, current_node, mode = "out")
    if (length(neigh) == 0) return(NA_character_)
    neigh_names <- V(g)$name[neigh]
    neigh_names <- neigh_names[!visited[neigh_names]]
    if (length(neigh_names) == 0) return(NA_character_)
    ov <- sapply(neigh_names, function(nn) {
      eid <- get.edge.ids(g, c(current_node, nn))
      E(g)$overlap[eid]
    })
    best_ov <- max(ov)
    cand <- neigh_names[ov == best_ov]
    if (length(cand) == 1) return(cand)
    if (isTRUE(verbose)) {
      cat(sprintf("\n[BRANCH] current=%s best_overlap=%d candidates=%s\n",
                  current_node, best_ov, paste(cand, collapse = ", ")))
    }
    # Single cohort: enrichment & delta
    if (!is.null(kmer_count_df) && all(c("enrichment", "delta") %in% names(kmer_count_df))) {
      enrichment <- setNames(kmer_count_df$enrichment, kmer_count_df$kmer)
      delta <- setNames(kmer_count_df$delta, kmer_count_df$kmer)
      cand_enrichment <- enrichment[cand]
      cand_delta <- delta[cand]
      ord <- order(-cand_enrichment, cand)
      cand <- cand[ord]
      top <- cand[1]; second <- cand[2]
      enrichment_ratio <- if (enrichment[second] == 0) Inf else (enrichment[top] / enrichment[second])
      delta_diff <- delta[top] - delta[second]
      ratio_ok <- enrichment_ratio >= enrichment_threshold
      diff_ok <- delta_diff >= delta_threshold
      if ((ratio_ok && diff_ok)) {
        if (isTRUE(verbose)) cat(sprintf("[BRANCH] PICK %s\n", top))
        return(top)
      }
      if (tie_action == "pick_lexical") return(sort(cand)[1])
      return(NA_character_)
    }
    # Multi-cohort: enrichment_x/y and delta_x/y
    if (!is.null(kmer_count_df) &&
        all(c("enrichment_x", "enrichment_y", "delta_x", "delta_y") %in% names(kmer_count_df))) {
      enrichment <- setNames(pmin(kmer_count_df$enrichment_x, kmer_count_df$enrichment_y), kmer_count_df$kmer)
      delta <- setNames(pmin(kmer_count_df$delta_x, kmer_count_df$delta_y), kmer_count_df$kmer)
      cand_enrichment <- enrichment[cand]
      cand_delta <- delta[cand]
      ord <- order(-cand_enrichment, cand)
      cand <- cand[ord]
      top <- cand[1]; second <- cand[2]
      enrichment_ratio <- if (enrichment[second] == 0) Inf else (enrichment[top] / enrichment[second])
      delta_diff <- delta[top] - delta[second]
      ratio_ok <- enrichment_ratio >= enrichment_threshold
      diff_ok <- delta_diff >= delta_threshold
      if ((ratio_ok && diff_ok)) {
        if (isTRUE(verbose)) cat(sprintf("[BRANCH] PICK %s\n", top))
        return(top)
      }
      if (tie_action == "pick_lexical") return(sort(cand)[1])
      return(NA_character_)
    }
    # No counts: pure overlap/lexical
    if (tie_action == "pick_lexical") return(sort(cand)[1])
    if (isTRUE(verbose)) cat("[BRANCH] no counts -> STOP\n")
    NA_character_
  }
  
  assemble_contig <- function(start_node) {
    contig <- start_node
    path_kmers <- c(start_node)
    current_node <- start_node
    while (TRUE) {
      next_node <- pick_next_node(current_node)
      if (is.na(next_node)) break
      eid <- get.edge.ids(g, c(current_node, next_node))
      overlap_size <- E(g)$overlap[eid]
      contig <- paste0(contig, substr(next_node, overlap_size + 1, nchar(next_node)))
      path_kmers <- c(path_kmers, next_node)
      visited[next_node] <<- TRUE
      current_node <- next_node
    }
    list(contig = contig, kmers_used = path_kmers)
  }
  
  contig_counter <- 1
  for (node_name in V(g)$name) {
    if (!visited[node_name]) {
      visited[node_name] <- TRUE
      result <- assemble_contig(node_name)
      contigs_info[[contig_counter]] <- list(
        contig = result$contig,
        kmers_used = result$kmers_used,
        num_kmers = length(result$kmers_used),
        is_isolated = length(result$kmers_used) == 1
      )
      contig_counter <- contig_counter + 1
    }
  }
  for (isolated_kmer in isolated_kmers) {
    contigs_info[[contig_counter]] <- list(
      contig = isolated_kmer,
      kmers_used = isolated_kmer,
      num_kmers = 1,
      is_isolated = TRUE
    )
    contig_counter <- contig_counter + 1
  }
  contigs_info
}

print_contig_info <- function(contigs_info) {
  for (i in seq_along(contigs_info)) {
    cat(sprintf("\n=== Contig %d ===\n", i))
    cat(sprintf("Sequence: %s\n", contigs_info[[i]]$contig))
    cat(sprintf("Length: %d bp\n", nchar(contigs_info[[i]]$contig)))
    cat(sprintf("Number of k-mers: %d\n", contigs_info[[i]]$num_kmers))
    cat(sprintf("Isolated: %s\n", contigs_info[[i]]$is_isolated))
    cat(sprintf("K-mers used: %s\n", paste(contigs_info[[i]]$kmers_used, collapse = ", ")))
  }
}