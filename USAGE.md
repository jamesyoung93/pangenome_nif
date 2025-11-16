# Pangenome Analysis Pipeline - Usage Guide

## Quick Start

To run the consolidated pangenome analysis pipeline on the HPC cluster:

### 1. Submit the SLURM job

You can submit the consolidated pipeline with either of these wrappers:

```bash
sbatch run_consolidated_pipeline.sh
```

or the new script that directly calls the all-in-one Python driver:

```bash
sbatch run_pipeline_consolidated.sh [path/to/input.csv]
```

If you omit the optional input path, `run_pipeline_consolidated.sh` defaults to
`nif_hdk_hits_enriched_with_quality_checkm.csv` in the repository directory.

Both wrappers will:
- Load all required modules (Python, MMseqs2, NCBI datasets)
- Set up the Python environment
- Run all 8 pipeline steps automatically
- Generate comprehensive output files

### 2. Monitor the job

```bash
# Check job status
squeue -u $USER

# View output in real-time (replace JOBID with your job ID)
tail -f pangenome_consolidated_JOBID.out
```

### 3. Check results

After completion, check the output file:
```bash
cat pangenome_consolidated_JOBID.out
```

## Important Notes

### File Location
**Always run from the repository directory!** The submission script will automatically:
- Change to the script directory
- Locate all required scripts (`pangenome_pipeline_consolidated.py`, `02_download_proteins.py`)
- Set up the environment properly

### Required Files
Before submitting, ensure these files are in your repository directory:
- `pangenome_pipeline_consolidated.py` - Main pipeline script
- `02_download_proteins.py` - Protein download script
- `nif_hdk_hits_enriched_with_quality_checkm.csv` - Input data
- `run_consolidated_pipeline.sh` - SLURM submission script

### NCBI Dataset Tools
The submission script will attempt to load NCBI dataset tools using common module names:
- `ncbi-datasets`
- `datasets`
- `ncbi`

If none are available, the pipeline will fall back to direct HTTPS downloads from NCBI FTP servers.

### Email Configuration
The pipeline is pre-configured with the email: `james.young@jacks.sdstate.edu`

This email is sent in the User-Agent header to NCBI as recommended by their API guidelines.

## Pipeline Steps

The consolidated pipeline runs these steps automatically:

1. **Data Preparation** - Load and label genomes as diazotroph/non-diazotroph
2. **Download Proteins** - Download protein sequences from NCBI for all genomes
3. **Gene Family Clustering** - Cluster proteins into gene families using MMseqs2
4. **Classification** - Train ML models with genus-level cross-validation
5. **Feature Directionality** - Analyze which features favor diazotrophs vs non-diazotrophs
6. **Swiss-Prot Annotation** - Merge functional annotations
7. **Gene Family Expansion** - Add detailed gene family information
8. **Narrative Table** - Create final interpretable results

## Customizing Resources

Edit `run_consolidated_pipeline.sh` to adjust:

```bash
#SBATCH --cpus-per-task=32    # Number of CPU cores
#SBATCH --mem=256G            # Memory allocation
#SBATCH --time=24:00:00       # Maximum run time
#SBATCH --partition=quickq    # SLURM partition
```

## Troubleshooting

### Download script not found
If you see:
```
WARNING: 02_download_proteins.py not found in any of the following locations
```

**Solution**: Make sure you're running from the repository directory and both Python scripts are present.

### NCBI download failures
If protein downloads fail:
1. Check `downloads/download_log.tsv` for details
2. Verify internet connectivity from compute nodes
3. Try running from a login node or data transfer node (DTN)
4. The pipeline will continue with available proteins

### Module not found
If required modules are missing:
```bash
module avail python
module avail mmseqs
module avail ncbi
```

Contact your HPC administrator if modules are unavailable.

## Advanced Usage

### Running individual steps
```bash
python3 pangenome_pipeline_consolidated.py \
    --input nif_hdk_hits_enriched_with_quality_checkm.csv \
    --threads 32 \
    --steps 3-5 \
    --non-interactive
```

### Interactive mode
For local testing (not on HPC):
```bash
python3 pangenome_pipeline_consolidated.py \
    --input nif_hdk_hits_enriched_with_quality_checkm.csv \
    --threads 8
```

## Output Files

Key output files include:
- `complete_genomes_labeled.csv` - Labeled genome dataset
- `downloads/proteins/*.faa` - Protein sequences
- `gene_family_matrix.csv` - Gene family presence/absence matrix
- `classification_summary.csv` - Model performance metrics
- `feature_directionality_full.csv` - Feature analysis results

See the main output file for complete details.
