# Snakemake script for combining GAM results from all chunks

suppressPackageStartupMessages({
  library(tidyverse)
})

# Get Snakemake parameters
output_summary <- snakemake@output[["summary"]]
output_contrasts <- snakemake@output[["contrasts"]]
results_dir <- snakemake@params[["results_dir"]]

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Combining GAM results from all chunks...")

# Find all chunk directories
chunk_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
chunk_dirs <- chunk_dirs[grepl("GAM_chunk_", chunk_dirs)]

log_msg(paste("Found", length(chunk_dirs), "chunk directories"))

# Combine summary files
log_msg("Combining summary statistics...")
all_summaries <- list()

for (chunk_dir in chunk_dirs) {
  summary_files <- list.files(chunk_dir, pattern = "_summary.csv$", full.names = TRUE)
  
  if (length(summary_files) > 0) {
    chunk_summaries <- lapply(summary_files, function(f) {
      tryCatch(
        read_csv(f, show_col_types = FALSE),
        error = function(e) NULL
      )
    })
    
    # Remove NULL entries (failed reads)
    chunk_summaries <- chunk_summaries[!sapply(chunk_summaries, is.null)]
    
    if (length(chunk_summaries) > 0) {
      all_summaries <- c(all_summaries, chunk_summaries)
    }
  }
}

if (length(all_summaries) > 0) {
  combined_summary <- bind_rows(all_summaries)
  write_csv(combined_summary, output_summary)
  log_msg(paste("Combined", nrow(combined_summary), "gene summaries"))
  
  # Print summary statistics
  log_msg("=== Overall Analysis Summary ===")
  log_msg(paste("Total genes analyzed:", nrow(combined_summary)))
  log_msg(paste("Successful analyses:", sum(combined_summary$Status == "Analysis completed", na.rm = TRUE)))
  log_msg(paste("Genes with significant contrasts:", 
                sum(combined_summary$Significant_Contrasts > 0, na.rm = TRUE)))
  log_msg(paste("Mean R-squared:", 
                round(mean(combined_summary$R_squared, na.rm = TRUE), 3)))
  log_msg(paste("Mean deviance explained:", 
                round(mean(combined_summary$Deviance_Explained, na.rm = TRUE) * 100, 2), "%"))
} else {
  log_msg("WARNING: No summary files found!")
  # Create empty file
  combined_summary <- data.frame(
    Gene = character(),
    AIC = numeric(),
    Deviance_Explained = numeric(),
    R_squared = numeric(),
    Significant_Contrasts = integer(),
    Total_Contrasts = integer(),
    Status = character()
  )
  write_csv(combined_summary, output_summary)
}

# Combine contrast files
log_msg("Combining contrast results...")
all_contrasts <- list()

for (chunk_dir in chunk_dirs) {
  contrast_files <- list.files(chunk_dir, pattern = "_contrasts.csv$", full.names = TRUE)
  
  if (length(contrast_files) > 0) {
    chunk_contrasts <- lapply(contrast_files, function(f) {
      tryCatch(
        read_csv(f, show_col_types = FALSE),
        error = function(e) NULL
      )
    })
    
    # Remove NULL entries
    chunk_contrasts <- chunk_contrasts[!sapply(chunk_contrasts, is.null)]
    
    if (length(chunk_contrasts) > 0) {
      all_contrasts <- c(all_contrasts, chunk_contrasts)
    }
  }
}

if (length(all_contrasts) > 0) {
  combined_contrasts <- bind_rows(all_contrasts)
  write_csv(combined_contrasts, output_contrasts)
}

log_msg("Results combination complete!")
  log_msg(paste("Combined", nrow(combined_contrasts), "contrast results"))
} else {
  log_msg("WARNING: No contrast files found!")
  # Create empty file with expected columns
  combined_contrasts <- data.frame(
    GOI = character(),
    DPA = numeric(),
    Level1 = character(),
    Level2 = character(),
    Difference = numeric(),
    CI_low = numeric(),
    CI_high = numeric(),
    SE = numeric(),
    p = numeric()
  )
  write_csv(combined_contrasts, output_contrasts)
