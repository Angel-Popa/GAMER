# Snakemake script for creating gene chunks

suppressPackageStartupMessages({
  library(tidyverse)
})

# Get Snakemake parameters
genes_file <- snakemake@input[["genes_list"]]
chunks_dir <- snakemake@output[["chunks_dir"]]
chunk_list_file <- snakemake@output[["chunk_list"]]
genes_per_chunk <- snakemake@params[["genes_per_chunk"]]

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Creating gene chunks for parallel processing...")

# Read gene IDs
gene_ids <- readLines(genes_file)
total_genes <- length(gene_ids)

log_msg(paste("Total genes to process:", total_genes))

# Calculate number of chunks
num_chunks <- ceiling(total_genes / genes_per_chunk)
log_msg(paste("Creating", num_chunks, "chunks with ~", genes_per_chunk, "genes each"))

# Create output directory
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)

# Split genes into chunks and write to files
chunk_files <- character(num_chunks)

for (i in 1:num_chunks) {
  # Calculate indices for this chunk
  start_idx <- (i - 1) * genes_per_chunk + 1
  end_idx <- min(i * genes_per_chunk, total_genes)
  
  # Extract chunk
  gene_chunk <- gene_ids[start_idx:end_idx]
  
  # Create chunk number with leading zeros
  chunk_num <- sprintf("%03d", i)
  
  # Create filename
  filename <- file.path(chunks_dir, paste0("genes_chunk_", chunk_num, ".txt"))
  chunk_files[i] <- filename
  
  # Write chunk to file
  writeLines(gene_chunk, filename)
  
  if (i %% 100 == 0 || i == num_chunks) {
    log_msg(paste("Created", i, "of", num_chunks, "chunks"))
  }
}

# Write list of chunk files
writeLines(chunk_files, chunk_list_file)

log_msg(paste("Chunk creation complete! Created", num_chunks, "chunk files"))
log_msg(paste("Chunk list saved to:", chunk_list_file))
