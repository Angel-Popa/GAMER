"""
Snakemake pipeline for GAM (Generalized Additive Model) analysis
of gene expression time series data
"""

import os
import math
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

configfile: "config/config.yaml"

# Extract configuration parameters
METADATA = config["metadata"]
COUNT_MATRIX = config["count_matrix"]
OUTPUT_DIR = config["output_dir"]
MIN_SAMPLES = config.get("min_samples", 12)
NORMALIZATION = config.get("normalization", "TMM")
ALPHA = config.get("alpha", 0.05)
P_THRESHOLD = config.get("p_threshold", 0.001)
SIG_THRESHOLD = config.get("sig_threshold", 0.05)
GENES_PER_CHUNK = config.get("genes_per_chunk", 20)
PARALLEL_JOBS = config.get("parallel_jobs", 4)

# Create output directories
RESULTS_DIR = f"{OUTPUT_DIR}/results"
REPORT_DIR = f"{OUTPUT_DIR}/report"
LOG_DIR = f"{OUTPUT_DIR}/logs"

# ============================================================================
# Helper Functions
# ============================================================================

def get_num_chunks(wildcards):
    """Calculate number of gene chunks dynamically"""
    genes_file = f"{REPORT_DIR}/genes_ids.txt"
    if os.path.exists(genes_file):
        with open(genes_file) as f:
            num_genes = sum(1 for _ in f)
        return math.ceil(num_genes / GENES_PER_CHUNK)
    return 1

# ============================================================================
# Target Rules
# ============================================================================

rule all:
    input:
        # Final outputs
        f"{REPORT_DIR}/GAM_summary.csv",
        f"{REPORT_DIR}/GAM_contrasts.csv",
        f"{REPORT_DIR}/GAM_significant_genes.txt",
        f"{REPORT_DIR}/contrast_plots/plots_complete.flag"

# ============================================================================
# Data Preparation
# ============================================================================

rule prepare_normalized_data:
    """
    Normalize count data using TMM and prepare for GAM analysis
    """
    input:
        metadata = METADATA,
        counts = COUNT_MATRIX
    output:
        normalized = f"{REPORT_DIR}/normalized_counts.Rdata",
        genes_list = f"{REPORT_DIR}/genes_ids.txt"
    params:
        min_samples = MIN_SAMPLES,
        report_dir = REPORT_DIR
    log:
        f"{LOG_DIR}/prepare_normalized_data.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/prepare_normalized_data.R"

rule create_gene_chunks:
    """
    Split gene list into chunks for parallel processing
    """
    input:
        genes_list = f"{REPORT_DIR}/genes_ids.txt"
    output:
        chunks_dir = directory(f"{REPORT_DIR}/gene_chunks"),
        chunk_list = f"{REPORT_DIR}/chunk_files.txt"
    params:
        genes_per_chunk = GENES_PER_CHUNK
    log:
        f"{LOG_DIR}/create_gene_chunks.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/create_gene_chunks.R"

# ============================================================================
# GAM Analysis
# ============================================================================

rule run_gam_analysis:
    """
    Run GAM analysis on a chunk of genes
    """
    input:
        normalized = f"{REPORT_DIR}/normalized_counts.Rdata",
        gene_chunk = f"{REPORT_DIR}/gene_chunks/genes_chunk_{{chunk}}.txt"
    output:
        results = directory(f"{RESULTS_DIR}/GAM_chunk_{{chunk}}"),
        done = touch(f"{RESULTS_DIR}/GAM_chunk_{{chunk}}.done")
    params:
        alpha = ALPHA,
        normalization = NORMALIZATION,
        parallel_jobs = PARALLEL_JOBS
    threads: PARALLEL_JOBS
    log:
        f"{LOG_DIR}/gam_analysis_chunk_{{chunk}}.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/run_gam_chunk.R"

checkpoint get_chunks:
    """
    Checkpoint to dynamically determine number of chunks
    """
    input:
        chunk_list = f"{REPORT_DIR}/chunk_files.txt"
    output:
        flag = touch(f"{REPORT_DIR}/chunks_ready.flag")
    run:
        pass

def aggregate_gam_chunks(wildcards):
    """
    Aggregate function to get all chunk outputs
    """
    checkpoint_output = checkpoints.get_chunks.get(**wildcards).output[0]
    chunk_list = f"{REPORT_DIR}/chunk_files.txt"

    with open(chunk_list) as f:
        chunks = [group.strip().split('_')[-1].replace('.txt', '')
                 for group in f if group.strip()]

    return expand(f"{RESULTS_DIR}/GAM_chunk_{{chunk}}.done", chunk=chunks)

rule combine_gam_results:
    """
    Combine all GAM results from chunks into summary files
    """
    input:
        aggregate_gam_chunks
    output:
        summary = f"{REPORT_DIR}/GAM_summary.csv",
        contrasts = f"{REPORT_DIR}/GAM_contrasts.csv"
    params:
        results_dir = RESULTS_DIR
    log:
        f"{LOG_DIR}/combine_gam_results.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/combine_gam_results.R"

# ============================================================================
# Visualization
# ============================================================================

rule identify_significant_genes:
    """
    Identify genes with significant contrasts
    """
    input:
        contrasts = f"{REPORT_DIR}/GAM_contrasts.csv",
        summary = f"{REPORT_DIR}/GAM_summary.csv"
    output:
        sig_genes = f"{REPORT_DIR}/GAM_significant_genes.txt"
    params:
        p_threshold = P_THRESHOLD
    log:
        f"{LOG_DIR}/identify_significant_genes.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/identify_significant_genes.R"

rule plot_contrasts:
    """
    Generate contrast plots for significant genes
    """
    input:
        sig_genes = f"{REPORT_DIR}/GAM_significant_genes.txt",
        summary = f"{REPORT_DIR}/GAM_summary.csv"
    output:
        flag = touch(f"{REPORT_DIR}/contrast_plots/plots_complete.flag")
    params:
        results_dir = RESULTS_DIR,
        plot_dir = f"{REPORT_DIR}/contrast_plots",
        p_threshold = P_THRESHOLD,
        sig_threshold = SIG_THRESHOLD
    log:
        f"{LOG_DIR}/plot_contrasts.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/plot_contrasts.R"

# ============================================================================
# Quality Control and Reporting
# ============================================================================

rule generate_qc_report:
    """
    Generate quality control report
    """
    input:
        summary = f"{REPORT_DIR}/GAM_summary.csv",
        contrasts = f"{REPORT_DIR}/GAM_contrasts.csv",
        sig_genes = f"{REPORT_DIR}/GAM_significant_genes.txt"
    output:
        report = f"{REPORT_DIR}/QC_report.html"
    params:
        output_dir = OUTPUT_DIR
    log:
        f"{LOG_DIR}/qc_report.log"
    conda:
        "envs/r_gam.yaml"
    script:
        "scripts/generate_qc_report.Rmd"
