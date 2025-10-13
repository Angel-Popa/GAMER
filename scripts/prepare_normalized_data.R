# Snakemake script for data preparation and normalization

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
})

# Get Snakemake parameters
metadata_file <- snakemake@input[["metadata"]]
counts_file <- snakemake@input[["counts"]]
output_normalized <- snakemake@output[["normalized"]]
output_genes <- snakemake@output[["genes_list"]]
min_samples <- snakemake@params[["min_samples"]]
report_dir <- snakemake@params[["report_dir"]]

# Get config parameters
config <- snakemake@config
lines_to_include <- config$lines_to_include
metadata_cols <- config$metadata_columns

# Create report directory if needed
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Starting data preparation and normalization...")

# Import metadata
log_msg(paste("Reading metadata from:", metadata_file))
metadata <- read_csv(metadata_file, col_names = TRUE, show_col_types = FALSE)

# Import count data
log_msg(paste("Reading count matrix from:", counts_file))
counts <- read_csv(counts_file, col_names = TRUE, show_col_types = FALSE)
counts_df <- as.data.frame(counts)
rownames(counts_df) <- counts_df[[1]]  # First column is gene ID
counts_df <- counts_df[, -1]  # Remove the gene_id column

log_msg(paste("Count matrix dimensions:", nrow(counts_df), "genes x", ncol(counts_df), "samples"))

# Order metadata to match count columns
sample_col <- metadata_cols$sample_id
metadata_ordered <- metadata[match(colnames(counts_df), metadata[[sample_col]]), ]

# Filter lowly expressed genes
log_msg("Filtering lowly expressed genes...")
keep <- rowSums(counts_df > 10) >= min_samples
filtered_counts <- counts_df[keep, ]
log_msg(paste("Retained", sum(keep), "of", nrow(counts_df), "genes after filtering"))

# Calculate TMM normalization factors
log_msg("Calculating TMM normalization factors...")
tmm <- calcNormFactors(filtered_counts, method = "TMM")
N <- colSums(filtered_counts)
tmm.counts <- N * tmm / exp(mean(log(N * tmm)))

# Normalize counts
log_msg("Normalizing counts...")
normalized_counts <- filtered_counts
for (sample in colnames(normalized_counts)) {
  normalized_counts[, sample] <- normalized_counts[, sample] / tmm.counts[sample]
}

# Prepare long format data for modeling
log_msg("Converting to long format...")
line_col <- metadata_cols$line
time_col <- metadata_cols$timepoint

normalized_counts_long <- normalized_counts %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = -c("gene_id"), names_to = sample_col, values_to = "expression") %>% 
  left_join(metadata_ordered, by = sample_col) %>% 
  mutate(!!sym(line_col) := as.factor(!!sym(line_col)))

# Filter for specific lines if specified
if (!is.null(lines_to_include)) {
  log_msg(paste("Filtering for lines:", paste(lines_to_include, collapse = ", ")))
  normalized_counts_long <- normalized_counts_long %>%
    filter(!!sym(line_col) %in% lines_to_include) %>%
    droplevels()
}

# Rename columns to standard names for downstream analysis
normalized_counts_long <- normalized_counts_long %>%
  rename(Line = !!sym(line_col),
         DPA = !!sym(time_col),
         Sample = !!sym(sample_col))

log_msg(paste("Final dataset:", 
              length(unique(normalized_counts_long$gene_id)), "genes,",
              length(unique(normalized_counts_long$Sample)), "samples,",
              length(unique(normalized_counts_long$Line)), "lines"))

# Save normalized data
log_msg(paste("Saving normalized data to:", output_normalized))
save(normalized_counts_long, file = output_normalized)

# Extract and save gene list
gene_ids <- unique(normalized_counts_long$gene_id)
log_msg(paste("Saving gene list to:", output_genes))
writeLines(gene_ids, output_genes)

log_msg("Data preparation complete!")
