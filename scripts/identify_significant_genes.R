# Snakemake script for identifying significant genes

suppressPackageStartupMessages({
  library(tidyverse)
})

# Get Snakemake parameters
contrasts_file <- snakemake@input[["contrasts"]]
summary_file <- snakemake@input[["summary"]]
output_file <- snakemake@output[["sig_genes"]]
p_threshold <- snakemake@params[["p_threshold"]]

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Identifying genes with significant contrasts...")

# Read data
contrasts <- read_csv(contrasts_file, show_col_types = FALSE)
summary <- read_csv(summary_file, show_col_types = FALSE)

# Filter for significant contrasts
log_msg(paste("Using p-value threshold:", p_threshold))

sig_contrasts <- contrasts %>%
  filter(p < p_threshold)

# Get unique genes with significant contrasts
genes_sig <- unique(sig_contrasts$GOI)

log_msg(paste("Found", length(genes_sig), "genes with significant contrasts"))

# Additional statistics
contrast_counts <- sig_contrasts %>%
  group_by(GOI) %>%
  summarise(
    n_sig_contrasts = n(),
    min_p = min(p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_contrasts))

log_msg("Top 10 genes by number of significant contrasts:")
print(head(contrast_counts, 10))

# Save significant genes
writeLines(genes_sig, output_file)
log_msg(paste("Significant genes saved to:", output_file))

log_msg("Identification complete!")
