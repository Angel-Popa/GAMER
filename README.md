# GAMER: Generalized Additive Model Expression Resolver Pipeline

A flexible Snakemake pipeline for Generalized Additive Model (GAM) analysis of gene expression time-series data.

## Overview

This pipeline performs:
1. **Data normalization** using TMM (Trimmed Mean of M-values)
2. **GAM modeling** to analyze gene expression patterns over time across different lines
3. **Contrast analysis** to identify significant differences between lines
4. **Visualization** of significant contrasts

## Features

- **Modular design**: Easy to adapt to different datasets
- **Parallel processing**: Processes genes in chunks for efficiency
- **Configurable**: All parameters controlled via YAML config file
- **Reproducible**: Conda environment management
- **Scalable**: Works with small to large datasets

## Directory Structure

```
.
├── config/
│   └── config.yaml           # Main configuration file
├── data/
│   ├── samples_metadata.csv  # Sample metadata (user-provided)
│   └── manual_gene_count_matrix.csv  # Count matrix (user-provided)
├── envs/
│   └── r_gam.yaml            # Conda environment specification
├── scripts/
│   ├── prepare_normalized_data.R
│   ├── create_gene_chunks.R
│   ├── run_gam_chunk.R
│   ├── combine_gam_results.R
│   ├── identify_significant_genes.R
│   └── plot_contrasts.R
├── Snakefile                  # Main workflow definition
└── README.md
```

## Requirements

- Snakemake (≥7.0)
- Conda/Mamba
- R (≥4.3)

## Installation

1. **Install Snakemake** (if not already installed):
```bash
conda install -c conda-forge -c bioconda snakemake
```

2. **Clone or download this pipeline**:
```bash
git clone <your-repo-url>
cd gam-pipeline
```

## Input Data Format

### 1. Sample Metadata (`samples_metadata.csv`)

Required columns (can be configured in `config.yaml`):
- **Sample**: Unique sample identifier
- **Line**: Line/genotype identifier
- **DPA**: Time point (e.g., Days Post-Anthesis)

Example:
```csv
Sample,Line,DPA,Other_Info
Sample1,1,10,batch1
Sample2,1,15,batch1
Sample3,2,10,batch1
Sample4,2,15,batch2
```

### 2. Count Matrix (`manual_gene_count_matrix.csv`)

Format:
- First column: gene IDs
- Remaining columns: sample counts (column names must match Sample IDs in metadata)

Example:
```csv
gene_id,Sample1,Sample2,Sample3,Sample4
Gene001,150,200,180,220
Gene002,50,45,52,48
Gene003,1200,1500,1300,1600
```

## Configuration

Edit `config/config.yaml` to customize the analysis:

### Key Parameters

```yaml
# Input files
metadata: "data/samples_metadata.csv"
count_matrix: "data/manual_gene_count_matrix.csv"

# Output directory
output_dir: "output"

# Statistical thresholds
alpha: 0.05              # Significance for contrasts
p_threshold: 0.001       # Filter for plotting
sig_threshold: 0.05      # Mark significance in plots

# Processing
genes_per_chunk: 20      # Genes per parallel chunk
parallel_jobs: 4         # Cores per chunk

# Metadata columns (customize to your data)
metadata_columns:
  sample_id: "Sample"
  line: "Line"
  timepoint: "DPA"

# Filter specific lines (or set to null for all)
lines_to_include: ["1", "2", "3", "4", "5", "6"]
```

## Usage

### Basic Execution

1. **Prepare your data** in the `data/` directory

2. **Edit configuration** in `config/config.yaml`

3. **Dry run** to check the workflow:
```bash
snakemake --configfile config/config.yaml --dry-run
```

4. **Run the pipeline**:
```bash
snakemake --configfile config/config.yaml --cores 8 --use-conda
```

### Advanced Options

**Run on a cluster** (SLURM example):
```bash
snakemake --configfile config/config.yaml \
  --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb} --time=02:00:00" \
  --jobs 20 \
  --use-conda
```

**Generate workflow diagram**:
```bash
snakemake --configfile config/config.yaml --dag | dot -Tpng > workflow.png
```

**Run specific steps**:
```bash
# Only data preparation
snakemake --configfile config/config.yaml output/report/normalized_counts.Rdata --use-conda

# Only GAM analysis (no plotting)
snakemake --configfile config/config.yaml output/report/GAM_summary.csv --use-conda
```

## Output Structure

```
output/
├── report/
│   ├── normalized_counts.Rdata       # Normalized expression data
│   ├── genes_ids.txt                 # List of all genes
│   ├── gene_chunks/                  # Gene chunks for parallel processing
│   ├── GAM_summary.csv               # Summary statistics for all genes
│   ├── GAM_contrasts.csv             # All contrast results
│   ├── GAM_significant_genes.txt     # Genes with significant contrasts
│   └── contrast_plots/               # Plots for significant genes
│       ├── Gene001_contrasts.png
│       └── ...
├── results/
│   ├── GAM_chunk_001/                # Results from each chunk
│   │   ├── Gene001_GAM_results.RData
│   │   ├── Gene001_summary.csv
│   │   ├── Gene001_contrasts.csv
│   │   └── ...
│   └── ...
└── logs/                             # Log files for each step
```

## Key Output Files

1. **GAM_summary.csv**: Summary statistics for each gene
   - AIC, R-squared, Deviance Explained
   - Number of significant contrasts

2. **GAM_contrasts.csv**: Detailed contrast results
   - Pairwise line comparisons at each time point
   - Effect sizes, confidence intervals, p-values

3. **GAM_significant_genes.txt**: List of genes with significant contrasts (p < threshold)

4. **Contrast plots**: Visual representation of significant differences between lines over time

## Customization for Different Datasets

### Example 1: Different Time Variable

If your time variable is "Hours" instead of "DPA":

```yaml
metadata_columns:
  sample_id: "SampleID"
  line: "Genotype"
  timepoint: "Hours"
```

### Example 2: All Lines Included

To analyze all lines in your dataset:

```yaml
lines_to_include: null
```

### Example 3: More Stringent Filtering

```yaml
alpha: 0.01
p_threshold: 0.0001
min_samples: 20
```

### Example 4: Larger Chunks for Fewer Jobs

```yaml
genes_per_chunk: 100
parallel_jobs: 8
```

## Troubleshooting

### Issue: "No genes pass filtering"
- Check `min_samples` parameter
- Verify count matrix has sufficient expression
- Review filtering criteria in logs

### Issue: "Memory errors during GAM fitting"
- Reduce `parallel_jobs` in config
- Reduce `genes_per_chunk`
- Increase memory allocation if using cluster

### Issue: "Missing contrasts package"
- Ensure conda environment is activated
- Try: `conda env update -f envs/r_gam.yaml`

### Issue: "Metadata/count matrix mismatch"
- Verify column names in count matrix match Sample IDs in metadata
- Check `metadata_columns` configuration matches your actual column names

## Citation

If you use this pipeline, please cite:
- **mgcv**: Wood, S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition)
- **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). Bioinformatics 26(1): 139-140
- **modelbased**: Makowski et al. (2020). CRAN

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review log files in `output/logs/`
3. Open an issue on GitHub

## License

MIT License - feel free to modify and adapt for your needs!
