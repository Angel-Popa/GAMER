# Snakemake script for running GAM analysis on a gene chunk

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(MASS)
  library(modelbased)
  library(parallel)
})

# Get Snakemake parameters
normalized_file <- snakemake@input[["normalized"]]
gene_chunk_file <- snakemake@input[["gene_chunk"]]
results_dir <- snakemake@output[["results"]]
alpha <- snakemake@params[["alpha"]]
n_cores <- snakemake@threads

# Get GAM parameters from config
config <- snakemake@config
gam_params <- config$gam_params
k <- gam_params$k
bs <- gam_params$bs

# Create results directory
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Starting GAM analysis for gene chunk...")

# Load normalized data
log_msg(paste("Loading normalized data from:", normalized_file))
load(normalized_file)

# Read gene list for this chunk
genes_in_chunk <- readLines(gene_chunk_file)
log_msg(paste("Processing", length(genes_in_chunk), "genes"))

# Function to analyze a single gene
analyze_gene <- function(gene_id) {
  tryCatch({
    # Filter data for this gene
    gene_data <- normalized_counts_long %>% 
      filter(gene_id == !!gene_id)
    
    # Check if data exists
    if (nrow(gene_data) == 0) {
      return(list(
        gene_id = gene_id,
        status = "No data",
        summary = data.frame(
          Gene = gene_id,
          AIC = NA,
          Deviance_Explained = NA,
          R_squared = NA,
          Significant_Contrasts = NA,
          Total_Contrasts = NA,
          Status = "No data available"
        )
      ))
    }
    
    # Fit GAM model
    gam_model <- mgcv::gam(
      expression ~ s(DPA, k = k, bs = bs, by = Line), 
      data = gene_data, 
      family = nb(link = "log")
    )
    
    # Extract model results
    estimate_relation <- estimate_relation(gam_model, length = 100, exponentiate = TRUE)
    parameters <- parameters::parameters(gam_model)
    
    # Perform contrasts analysis
    contrasts <- estimate_contrasts(
      gam_model,
      contrast = "Line",
      by = "DPA",
      length = 100,
      p_adjust = "BH",
      backend = "emmeans"
    )
    
    # Add gene identifier
    contrasts$GOI <- gene_id
    
    # Create summary statistics
    summary_stats <- data.frame(
      Gene = gene_id,
      AIC = AIC(gam_model),
      Deviance_Explained = summary(gam_model)$dev.expl,
      R_squared = summary(gam_model)$r.sq,
      Significant_Contrasts = sum(contrasts$p < alpha, na.rm = TRUE),
      Total_Contrasts = nrow(contrasts),
      Status = "Analysis completed",
      stringsAsFactors = FALSE
    )
    
    return(list(
      gene_id = gene_id,
      status = "Success",
      model = gam_model,
      estimate_relation = estimate_relation,
      parameters = parameters,
      contrasts = contrasts,
      summary = summary_stats
    ))
    
  }, error = function(e) {
    return(list(
      gene_id = gene_id,
      status = "Error",
      error_message = as.character(e),
      summary = data.frame(
        Gene = gene_id,
        AIC = NA,
        Deviance_Explained = NA,
        R_squared = NA,
        Significant_Contrasts = NA,
        Total_Contrasts = NA,
        Status = paste("Error:", as.character(e))
      )
    ))
  })
}

# Process genes in parallel
log_msg(paste("Running GAM analysis using", n_cores, "cores..."))
results <- mclapply(genes_in_chunk, analyze_gene, mc.cores = n_cores)

# Save individual results
log_msg("Saving individual gene results...")
for (result in results) {
  gene_id <- result$gene_id
  
  # Create GAM_results list structure
  GAM_results <- list()
  GAM_results[[gene_id]] <- result
  
  # Save results
  output_file <- file.path(results_dir, paste0(gene_id, "_GAM_results.RData"))
  save(GAM_results, file = output_file)
  
  # Save summary
  summary_file <- file.path(results_dir, paste0(gene_id, "_summary.csv"))
  write.csv(result$summary, summary_file, row.names = FALSE)
  
  # Save contrasts if available
  if (!is.null(result$contrasts)) {
    contrasts_file <- file.path(results_dir, paste0(gene_id, "_contrasts.csv"))
    write.csv(result$contrasts, contrasts_file, row.names = FALSE)
  }
}

# Summary statistics
n_success <- sum(sapply(results, function(x) x$status == "Success"))
n_error <- sum(sapply(results, function(x) x$status == "Error"))
n_nodata <- sum(sapply(results, function(x) x$status == "No data"))

log_msg("=== Chunk Analysis Summary ===")
log_msg(paste("Total genes:", length(genes_in_chunk)))
log_msg(paste("Successful:", n_success))
log_msg(paste("Errors:", n_error))
log_msg(paste("No data:", n_nodata))
log_msg("Chunk analysis complete!")
