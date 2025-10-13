# Quick Start Guide

Get the GAM pipeline running in 5 minutes!

## Step 1: Setup Environment

```bash
# Create conda environment with Snakemake
conda create -n snakemake -c conda-forge -c bioconda snakemake

# Activate environment
conda activate snakemake
```

## Step 2: Prepare Your Data

### Required Files:

1. **Sample metadata** (`data/samples_metadata.csv`)
   - Must include columns for: Sample ID, Line, Time point
   - See `data/samples_metadata_example.csv` for format

2. **Count matrix** (`data/manual_gene_count_matrix.csv`)
   - First column: gene IDs
   - Other columns: sample counts (names must match metadata)

### Example structure:

```
your_project/
├── data/
│   ├── samples_metadata.csv      ← Your metadata
│   └── gene_counts.csv            ← Your count matrix
```

## Step 3: Configure Pipeline

Edit `config/config.yaml`:

```yaml
# Point to your data files
metadata: "data/samples_metadata.csv"
count_matrix: "data/gene_counts.csv"

# Adjust column names to match YOUR data
metadata_columns:
  sample_id: "Sample"      ← Your sample ID column name
  line: "Line"             ← Your line/genotype column name
  timepoint: "DPA"         ← Your time variable column name

# Optional: Specify which lines to analyze
lines_to_include: ["1", "2", "3"]  # or null for all
```

## Step 4: Test Run

```bash
# Check if everything is configured correctly
./run_pipeline.sh --dry-run

# Or using snakemake directly
snakemake --configfile config/config.yaml --dry-run
```

If you see a list of jobs that will be executed, you're ready to go!

## Step 5: Run Pipeline

### Option A: Simple local run
```bash
./run_pipeline.sh --cores 8
```

### Option B: Using snakemake directly
```bash
snakemake --configfile config/config.yaml --cores 8 --use-conda
```

### Option C: On a SLURM cluster
```bash
./run_pipeline.sh --cluster --cores 20
```

## Step 6: Check Results

After completion, find results in:

```
output/
├── report/
│   ├── GAM_summary.csv              ← All gene statistics
│   ├── GAM_contrasts.csv            ← All contrast results
│   ├── GAM_significant_genes.txt    ← Significant genes list
│   └── contrast_plots/              ← Visualization plots
└── logs/                            ← Check if errors occur
```

## Common First-Time Issues

### Issue: "Metadata column not found"
**Solution**: Make sure `metadata_columns` in config matches your actual column names

### Issue: "No samples match between metadata and counts"
**Solution**: Check that sample IDs in count matrix column names exactly match Sample column in metadata

### Issue: "
