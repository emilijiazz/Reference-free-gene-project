library(readr)
source("assemble_kmers_finalversion.R")

df <- read_tsv("intersected_top_kmers.tsv")
df_selected <- df[df$top1 & df$top2, ]

kmer_count_df <- data.frame(
  kmer          = df_selected$kmer,
  enrichment_x  = df_selected$enrichment.x,
  delta_x       = df_selected$disease_count.x - df_selected$healthy_count.x,
  enrichment_y  = df_selected$enrichment.y,
  delta_y       = df_selected$disease_count.y - df_selected$healthy_count.y,
  stringsAsFactors = FALSE
)

contigs <- assemble_kmers(
  kmers = kmer_count_df$kmer,
  kmer_count_df = kmer_count_df,
  enrichment_threshold = 1.25,
  delta_threshold = 200,
  tie_action = "stop"
)

print_contig_info(contigs)