# SDSU Innovator HPC Setup and Execution Instructions

## Consolidated Pangenome Analysis Pipeline

This document provides step-by-step instructions for running the consolidated pangenome analysis pipeline on the SDSU Innovator HPC system from scratch.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Initial Setup](#initial-setup)
3. [Uploading Data](#uploading-data)
4. [Running the Pipeline](#running-the-pipeline)
5. [Monitoring the Job](#monitoring-the-job)
6. [Retrieving Results](#retrieving-results)
7. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Files

Before starting, ensure you have these files ready:

1. **Input data file**: `nif_hdk_hits_enriched_with_quality_checkm.csv`
2. **Main pipeline script**: `pangenome_pipeline_consolidated.py`
3. **SLURM submission script**: `run_consolidated_pipeline.slurm`

### HPC Account

- Active account on SDSU Innovator HPC
- SSH access configured
- Basic familiarity with Linux command line

---

## Initial Setup

### Step 1: Connect to SDSU Innovator HPC

Open a terminal and connect via SSH:

```bash
ssh your_username@innovator.sdsu.edu
```

Replace `your_username` with your actual SDSU username.

### Step 2: Create a Working Directory

Once logged in, create a dedicated directory for this analysis:

```bash
cd $HOME
mkdir -p pangenome_nif_analysis
cd pangenome_nif_analysis
```

### Step 3: Load required modules

Load the base environment used by the SLURM scripts. On Innovator the NCBI
tools module is published as **`ncbi-blast/2.14.1`** (not `ncbi`):

```bash
module purge
module load StdEnv
module load python/3.12
module load mmseqs2/15-6f452
# For NCBI utilities
module keyword ncbi        # shows available NCBI-related modules
module load ncbi-blast/2.14.1
```

If your site uses a different module name, run `module keyword ncbi` and load
the closest match. The pipeline will fall back to direct HTTPS downloads if no
NCBI module is available.

### Step 4: Create and activate your Python environment (required)

The submission scripts no longer create a virtual environment or install packages. Activate your own venv/conda environment **before** running the SLURM scripts so the required packages are available.

Example using `venv`:
```bash
module load python/3.12   # or your preferred Python 3.9+
python3 -m venv ~/pangenome_venv
source ~/pangenome_venv/bin/activate
pip install -U pip pandas numpy scipy scikit-learn matplotlib seaborn xgboost
```

Example using conda:
```bash
conda create -n pangenome python=3.11 -y
conda activate pangenome
conda install -y pandas numpy scipy scikit-learn matplotlib seaborn xgboost
```

---

## Uploading Data

### Step 5: Transfer Files to HPC

From your **local machine** (in a new terminal window), upload the required files:

```bash
# Upload input data
scp nif_hdk_hits_enriched_with_quality_checkm.csv your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/

# Upload pipeline scripts
scp pangenome_pipeline_consolidated.py your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/
scp run_consolidated_pipeline.slurm your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/

# If you have the original scripts (optional for annotation steps):
scp 06_merge_annotations_selected.py your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/
scp 07_expand_gene_families.py your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/
scp 08_build_narrative_table.py your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/
```

### Alternative: Git Clone

If your files are in a Git repository:

```bash
cd ~/pangenome_nif_analysis
git clone https://github.com/your_repo/pangenome_nif.git .
```

---

## Running the Pipeline

### Step 6: Verify Files

Back in your SSH session on the HPC, verify all files are present:

```bash
cd ~/pangenome_nif_analysis
ls -lh
```

You should see:
- `nif_hdk_hits_enriched_with_quality_checkm.csv`
- `pangenome_pipeline_consolidated.py`
- `run_consolidated_pipeline.slurm`

### Step 7: Make Scripts Executable

```bash
chmod +x pangenome_pipeline_consolidated.py
chmod +x run_consolidated_pipeline.slurm
```

### Step 8: Review SLURM Settings (Optional)

The SLURM script is configured for the **bigmem** partition with:
- **CPUs**: 32 cores
- **Memory**: 256 GB
- **Time limit**: 48 hours
- **Partition**: bigmem

To modify these settings, edit `run_consolidated_pipeline.slurm`:

```bash
nano run_consolidated_pipeline.slurm
```

Adjust these lines as needed:
```bash
#SBATCH --cpus-per-task=32      # Number of CPU cores
#SBATCH --mem=256G               # Memory allocation
#SBATCH --time=48:00:00          # Maximum runtime (HH:MM:SS)
#SBATCH --partition=bigmem       # Partition name
```

Save with `Ctrl+O`, `Enter`, then exit with `Ctrl+X`.

### Step 9: Submit the Job

Submit the pipeline to the SLURM scheduler:

```bash
# If using a conda environment, export its name so the batch job activates it
CONDA_ENV_NAME=pangenome2 sbatch run_consolidated_pipeline.slurm

# Without conda, omit CONDA_ENV_NAME (the job will use the current Python modules)
# sbatch run_consolidated_pipeline.slurm
```

You should see output like:
```
Submitted batch job 1234567
```

Make note of the job ID (e.g., `1234567`).

---

## Monitoring the Job

### Check Job Status

Monitor your job status:

```bash
squeue -u your_username
```

Output will show:
```
JOBID   PARTITION   NAME                    USER    ST  TIME  NODES
1234567 bigmem      pangenome_consolidated  user    R   2:30  1
```

Status codes:
- **PD**: Pending (waiting for resources)
- **R**: Running
- **CG**: Completing
- **CD**: Completed

### View Real-Time Output

Monitor the output log as the job runs:

```bash
tail -f pangenome_consolidated_1234567.out
```

Press `Ctrl+C` to stop viewing.

### Check for Errors

If you encounter issues, check the error log:

```bash
tail -f pangenome_consolidated_1234567.err
```

### Detailed Job Information

Get detailed information about your job:

```bash
scontrol show job 1234567
```

### Cancel a Job (if needed)

If you need to cancel the job:

```bash
scancel 1234567
```

---

## Retrieving Results

### Step 10: Check Job Completion

Once the job completes (disappears from `squeue`), verify it finished successfully:

```bash
tail -50 pangenome_consolidated_1234567.out
```

Look for:
```
====================================================================
PIPELINE COMPLETED SUCCESSFULLY!
====================================================================
```

### Step 11: Review Results

List all output files:

```bash
ls -lh *.csv *.png *.tsv
```

Key output files include:
- `complete_genomes_labeled.csv`
- `gene_family_matrix.csv`
- `classification_summary.csv`
- `roc_curves.png`
- `feature_importance_*.csv`
- `feature_directionality_full.csv`
- `feature_effects.png`
- `narrative_top_features.csv` (if step 8 completed)

### Step 10: Download Results Archive

The pipeline automatically creates a compressed archive. Check for:

```bash
ls -lh results_*.tar.gz
```

Download the archive to your local machine (from local terminal):

```bash
scp your_username@innovator.sdsu.edu:~/pangenome_nif_analysis/results_*.tar.gz ./
```

### Step 11: Extract Results Locally

On your local machine:

```bash
tar -xzf results_20241113.tar.gz
cd results_20241113
```

---

## Expected Runtime

Based on typical datasets:

| Step | Description | Estimated Time |
|------|-------------|----------------|
| 1    | Data preparation | 5-10 minutes |
| 2    | Protein download | 2-8 hours |
| 3    | Gene family clustering | 1-3 hours |
| 4    | Classification | 30-90 minutes |
| 5    | Feature directionality | 5-15 minutes |
| 6-8  | Annotation and tables | 10-30 minutes |

**Total estimated time**: 4-12 hours (depending on dataset size and download speeds)

---

## Pipeline Steps Overview

The consolidated pipeline runs 8 sequential steps:

1. **Data Preparation**: Filters for complete genomes and labels diazotrophs
2. **Protein Download**: Downloads protein sequences from NCBI
3. **Gene Family Creation**: Clusters proteins using MMseqs2 (80% identity)
4. **Classification**: Trains ML models with genus-level cross-validation
5. **Feature Directionality**: Analyzes positive/negative associations
6. **Merge Annotations**: Adds Swiss-Prot annotations (optional)
7. **Expand Gene Families**: Generates detailed member lists (optional)
8. **Build Narrative Table**: Creates final annotated feature table

---

## Troubleshooting

### Issue: Job Fails Immediately

**Check error log**:
```bash
cat pangenome_consolidated_1234567.err
```

**Common causes**:
- Missing input file
- Insufficient memory
- Module load failures

**Solutions**:
- Verify input file exists: `ls -lh nif_hdk_hits_enriched_with_quality_checkm.csv`
- Check available modules: `module avail`
- Request more memory in SLURM script

### Issue: Job Pending for Long Time

**Check queue status**:
```bash
squeue -p bigmem
```

**If bigmem partition is busy**:
- Wait for resources to become available
- Try a different partition (may need to reduce memory request)
- Contact HPC support for queue status

### Issue: Protein Download Fails

**Symptoms**: Step 2 reports many failed downloads

**Solutions**:
1. Download manually or on login node (better internet access)
2. Use existing protein files if available
3. Run with `--skip-download` flag and provide proteins separately

### Issue: MMseqs2 Clustering Runs Out of Memory

**Solutions**:
- Increase memory allocation in SLURM script: `#SBATCH --mem=512G`
- Reduce coverage parameter: `--min-genomes 20`
- Filter input to fewer genomes

### Issue: XGBoost Not Available

**This is expected** - the pipeline will continue without XGBoost and use:
- Random Forest
- Gradient Boosting
- Logistic Regression

### Issue: Steps 6-8 Skip or Fail

**This is normal** if:
- Swiss-Prot annotation files are not available
- Optional steps can be run separately later if needed

---

## Advanced Options

### Run Specific Steps Only

Modify the Python script call in the SLURM script:

```bash
# Run only steps 1-5 (skip annotation)
python3 pangenome_pipeline_consolidated.py \
    --input "$INPUT_FILE" \
    --threads $SLURM_CPUS_PER_TASK \
    --start-step 1 \
    --end-step 5 \
    --non-interactive
```

### Skip Protein Download

If proteins are already downloaded:

```bash
python3 pangenome_pipeline_consolidated.py \
    --input "$INPUT_FILE" \
    --threads $SLURM_CPUS_PER_TASK \
    --skip-download \
    --non-interactive
```

### Adjust Parameters

```bash
python3 pangenome_pipeline_consolidated.py \
    --input "$INPUT_FILE" \
    --threads $SLURM_CPUS_PER_TASK \
    --min-genomes 20 \      # Lower threshold for small datasets
    --non-interactive
```

### Non-Interactive Mode (Recommended for SLURM)

When running in SLURM batch jobs, use the `--non-interactive` flag to prevent the pipeline from prompting for user input. This flag is already included in the provided SLURM script, but can be added manually if needed:

```bash
python3 pangenome_pipeline_consolidated.py \
    --input "$INPUT_FILE" \
    --threads $SLURM_CPUS_PER_TASK \
    --non-interactive        # No user prompts, auto-continue on warnings
```

This flag:
- Prevents `EOFError` failures when running in batch mode
- Automatically continues even when protein files are limited
- Is recommended for all SLURM/batch job submissions

---

## Getting Help

### HPC-Specific Issues

Contact SDSU HPC Support:
- Email: hpc-support@sdsu.edu
- Check HPC documentation: https://innovator.sdsu.edu/documentation

### Pipeline-Specific Issues

Check the following files for diagnostic information:
1. `pangenome_consolidated_JOBID.out` - Standard output log
2. `pangenome_consolidated_JOBID.err` - Error log
3. `download_log.tsv` - Protein download status (if step 2 ran)

### Resource Monitoring

Check your disk quota:
```bash
df -h $HOME
```

Check partition limits:
```bash
sinfo -p bigmem
```

---

## Best Practices

1. **Test First**: Consider testing with a small subset of genomes first
2. **Monitor Resources**: Check memory and CPU usage during runs
3. **Backup Results**: Copy important results to long-term storage
4. **Document Parameters**: Keep notes on parameters used for reproducibility
5. **Archive Intermediate Files**: Save clustering results for reuse

---

## Quick Reference Commands

```bash
# Submit job
sbatch run_consolidated_pipeline.slurm

# Check job status
squeue -u $USER

# View live output
tail -f pangenome_consolidated_*.out

# Check for errors
tail -f pangenome_consolidated_*.err

# Cancel job
scancel JOBID

# Check disk space
df -h .

# Download results
scp user@innovator.sdsu.edu:~/pangenome_nif_analysis/results_*.tar.gz ./
```

---

## Success Indicators

Your pipeline has completed successfully when you see:

1. **Job status**: Completed (not in `squeue` output)
2. **Output log** shows: "PIPELINE COMPLETED SUCCESSFULLY!"
3. **Results archive** created: `results_YYYYMMDD.tar.gz`
4. **Key files present**:
   - `classification_summary.csv`
   - `roc_curves.png`
   - `feature_directionality_full.csv`

---

## Next Steps After Completion

1. **Review Results**: Examine `SUMMARY.txt` in results archive
2. **Visualize**: Open PNG files to view plots
3. **Analyze Features**: Review `narrative_top_features.csv` for biological insights
4. **Validate**: Check model performance metrics in `classification_summary.csv`
5. **Iterate**: Adjust parameters and rerun if needed

---

## Contact Information

For questions or issues specific to this pipeline, please contact:
- James Young: james.young@jacks.sdstate.edu

For SDSU HPC system issues:
- HPC Support: hpc-support@sdsu.edu

---

**Document Version**: 1.0
**Last Updated**: November 2024
**Pipeline Version**: Consolidated v1.0
