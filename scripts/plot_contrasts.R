# Snakemake script for plotting contrast results

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Get Snakemake parameters
sig_genes_file <- snakemake@input[["sig_genes"]]
summary_file <- snakemake@input[["summary"]]
results_dir <- snakemake@params[["results_dir"]]
plot_dir <- snakemake@params[["plot_dir"]]
p_threshold <- snakemake@params[["p_threshold"]]
sig_threshold <- snakemake@params[["sig_threshold"]]

# Get plotting parameters from config
config <- snakemake@config
plot_params <- config$plot_params
base_size <- plot_params$base_size
width <- plot_params$width
height <- plot_params$height
dpi <- plot_params$dpi

# Create plot directory
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Log function
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_msg("Generating contrast plots for significant genes...")

# Read significant genes
genes_sig <- readLines(sig_genes_file)
log_msg(paste("Plotting", length(genes_sig), "significant genes"))

# Function to plot contrasts for a single gene
plot_gene_contrasts <- function(gene_id, 
                                results_dir,
                                p_threshold,
                                sig_threshold,
                                base_size) {
  
  tryCatch({
    # Find all chunk directories and search for the gene
    chunk_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
    chunk_dirs <- chunk_dirs[grepl("GAM_chunk_", chunk_dirs)]
    
    # Look for the gene's results file
    result_file <- NULL
    for (chunk_dir in chunk_dirs) {
      potential_file <- file.path(chunk_dir, paste0(gene_id, "_GAM_results.RData"))
      if (file.exists(potential_file)) {
        result_file <- potential_file
        break
      }
    }
    
    if (is.null(result_file)) {
      warning(paste("Results file not found for gene:", gene_id))
      return(NULL)
    }
    
    # Load the GAM results
    load(result_file)
    
    # Extract and process contrast data
    contrast_sig <- as_tibble(GAM_results[[gene_id]]$contrasts) %>%
      mutate(
        DPA = floor(DPA),
        sig_symbol = ifelse(p < sig_threshold, "*", ""),
        Contrast = paste0(Level1, " vs ", Level2)
      ) %>% 
      filter(p < p_threshold)
    
    if (nrow(contrast_sig) == 0) {
      warning(paste("No significant contrasts for gene:", gene_id))
      return(NULL)
    }
    
    # Create the plot
    plot <- ggplot(contrast_sig, aes(x = DPA, y = Difference)) +
      geom_ribbon(aes(fill = Contrast, ymin = CI_low, ymax = CI_high), alpha = 0.2) +
      geom_pointrange(aes(colour = Contrast, ymin = CI_low, ymax = CI_high)) + 
      geom_line(aes(colour = Contrast), linewidth = 1) +
      geom_text(aes(label = sig_symbol), vjust = -0.5, size = 5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw(base_size = base_size) +
      ylab("Difference") + 
      xlab("DPA (Days Post-Anthesis)") +
      facet_wrap(~Contrast) + 
      ggtitle(paste0("Significant contrast differences for ", gene_id)) + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.position = "none"
      )
    
    return(plot)
    
  }, error = function(e) {
    warning(paste("Error plotting gene", gene_id, ":", e$message))
    return(NULL)
  })
}

# Generate plots for all significant genes
log_msg("Creating individual contrast plots...")

successful_plots <- 0
failed_plots <- 0

for (i in seq_along(genes_sig)) {
  gene_id <- genes_sig[i]
  
  plot <- plot_gene_contrasts(
    gene_id = gene_id,
    results_dir = results_dir,
    p_threshold = p_threshold,
    sig_threshold = sig_threshold,
    base_size = base_size
  )
  
  if (!is.null(plot)) {
    # Save plot
    plot_file <- file.path(plot_dir, paste0(gene_id, "_contrasts.png"))
    ggsave(
      filename = plot_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi
    )
    successful_plots <- successful_plots + 1
  } else {
    failed_plots <- failed_plots + 1
  }
  
  # Progress update every 100 genes
  if (i %% 100 == 0) {
    log_msg(paste("Processed", i, "of", length(genes_sig), "genes"))
  }
}

log_msg("=== Plotting Summary ===")
log_msg(paste("Total genes:", length(genes_sig)))
log_msg(paste("Successful plots:", successful_plots))
log_msg(paste("Failed plots:", failed_plots))
log_msg(paste("Plots saved to:", plot_dir))

log_msg("Contrast plotting complete!")
