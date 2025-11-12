# Diazotroph Classification Pipeline

A complete machine learning pipeline for classifying cyanobacteria as diazotrophs or non-diazotrophs based on gene family presence/absence patterns.

## Overview

This pipeline performs the following steps:
1. Filters for ONLY "Complete Genome" assembly level and labels them as diazotrophs (all nifH, D, K e-values < 1e-50)
2. Downloads protein sequences from NCBI (see troubleshooting if this fails)
3. Clusters proteins into gene families using MMseqs2 (80% identity, 80% coverage)
4. Creates presence/absence matrix and filters for families with ≥40 genomes (or 10% of available)
5. Removes highly pure gene families (>90% or <10% diazotroph)
6. Runs multiple classification models with genus-level cross-validation
7. Analyzes feature importance and directional effects

**IMPORTANT:** Protein downloads often fail on HPC compute nodes due to network
restrictions. See "Troubleshooting Downloads" section below for solutions.

## Requirements

### HPC Modules
- python/3.12
- mmseqs2/15-6f452

### Python Packages (auto-installed)
- pandas
- numpy
- scipy
- scikit-learn
- matplotlib
- seaborn
- xgboost

## Input File

**nif_hdk_hits_enriched_with_quality_checkm.csv**
- CSV file containing genome metadata and nif gene annotations
- Required columns:
  - assembly_accession
  - nifH_best_evalue, nifD_best_evalue, nifK_best_evalue
  - organism_name (for genus extraction)
  - assembly_level
  - checkm_completeness, checkm_contamination

## Usage

### Quick Start

1. Upload your data and pipeline to the HPC:
```bash
# Transfer files to HPC
scp -r diazotroph_classification/ your_username@hpc_address:~/
scp nif_hdk_hits_enriched_with_quality_checkm.csv your_username@hpc_address:~/diazotroph_classification/
```

2. Submit the job:
```bash
cd ~/diazotroph_classification
sbatch run_pipeline.sh
```

3. Monitor progress:
```bash
# Check job status
squeue -u your_username

# Watch output in real-time
tail -f diazotroph_classification_*.out
```

4. Retrieve results when complete:
```bash
# Download results archive
scp your_username@hpc_address:~/diazotroph_classification/results.tar.gz .

# Extract locally
tar -xzf results.tar.gz
```

### Running Individual Steps

You can also run steps individually for testing or debugging:

```bash
# Step 1: Prepare data
python3 01_prepare_data.py nif_hdk_hits_enriched_with_quality_checkm.csv

# Step 2: Download proteins (can take several hours)
python3 02_download_proteins.py

# Step 3: Create gene families (specify number of threads)
python3 03_create_gene_families.py 16

# Step 4: Classification
python3 04_classify.py

# Step 5: Feature directionality
python3 05_analyze_feature_directionality.py
```

## Output Files

### Primary Results

**classification_summary.csv**
- Performance metrics for all models (accuracy, precision, recall, F1, ROC-AUC)
- Mean ± standard deviation across genus-level CV folds

**roc_curves.png**
- ROC curves for all classification models
- Shows AUC with standard deviation

**feature_directionality_summary.csv**
- Top 50 most important features
- Direction (positive = diazotroph-associated, negative = non-diazotroph-associated)
- Effect sizes and statistical significance

**feature_effects.png**
- Visualization of feature importance with directionality
- Green = diazotroph-associated, Red = non-diazotroph-associated

### Intermediate Files

**complete_genomes_labeled.csv**
- Filtered dataset of complete genomes with diazotroph labels

**gene_family_matrix.csv**
- Presence/absence matrix for all gene families (≥40 genomes)
- Rows = genomes, Columns = gene families

**gene_family_purity_stats.csv**
- Statistics for each gene family
- Identifies pure vs. informative families

**feature_importance_[MODEL].csv**
- Feature importance scores for each classification model

**fold_metrics_[MODEL].csv**
- Per-fold performance metrics for each model

## Classification Models

The pipeline runs four classification models:

1. **Random Forest** (n=100 trees, max_depth=10)
2. **Gradient Boosting** (n=100 estimators, learning_rate=0.1)
3. **Logistic Regression** (L2 regularization)
4. **XGBoost** (n=100 estimators, max_depth=5)

All models use:
- **Genus-level cross-validation** (5 folds)
  - Ensures all genomes from the same genus are in the same fold
  - Tests model generalization to new genera
- Binary classification (diazotroph vs. non-diazotroph)
- Feature filtering (removing >90% or <10% pure families)

## Understanding the Results

### Classification Performance

Check `classification_summary.csv` for model comparison:
- **Accuracy**: Overall correctness
- **Precision**: When model predicts diazotroph, how often is it correct?
- **Recall**: What fraction of diazotrophs are correctly identified?
- **F1**: Harmonic mean of precision and recall
- **ROC-AUC**: Area under ROC curve (discrimination ability)

Values closer to 1.0 are better. Standard deviations indicate consistency across genera.

### Feature Importance

Top features in `feature_directionality_summary.csv` show:
- **Positive direction**: Gene family present more often in diazotrophs
- **Negative direction**: Gene family present more often in non-diazotrophs
- **Effect size**: Change in diazotroph rate when feature is present

### Genus-Level Cross-Validation

This approach:
- Tests if the model can predict diazotrophy for new, unseen genera
- More realistic evaluation than random splitting
- Lower performance than standard CV, but more reliable for generalization

## Troubleshooting

### Troubleshooting Downloads

**Problem:** All protein downloads fail with "Failed (download error)"

**Cause:** HPC compute nodes often cannot access external networks

**Solutions (in order of preference):**

1. **Run from login node** (MOST COMMON SOLUTION)
   ```bash
   # SSH to HPC
   cd ~/diazotroph_classification
   # Run download from login node, not via sbatch
   python3 02_download_proteins.py
   ```

2. **Use data transfer node** (if available)
   ```bash
   ssh username@data-transfer.hpc.edu
   cd ~/diazotroph_classification
   python3 02_download_proteins.py
   ```

3. **Manual download methods**
   See DOWNLOAD_TROUBLESHOOTING.txt for detailed instructions including:
   - NCBI datasets CLI batch download
   - wget/curl scripts
   - Local download + transfer

4. **Continue with partial data**
   If you have 100+ genomes downloaded, you can proceed:
   ```bash
   ls protein_sequences/*.faa | wc -l  # Check count
   python3 03_create_gene_families.py 16
   ```
   The pipeline will automatically adjust to available data.

### Other Issues

### Job runs out of memory
- Increase `#SBATCH --mem=` in run_pipeline.sh (try 128G)
- Or reduce dataset size in Step 1

### Protein downloads fail
- Check internet connectivity from compute node
- Some genomes may not have protein sequences available
- Failed downloads are logged in `failed_downloads.txt`
- Pipeline will continue with available genomes

### MMseqs2 clustering is slow
- Increase `#SBATCH --cpus-per-task=` for more threads
- Default 16 threads should complete in 1-2 hours for ~900 genomes

### No XGBoost model
- XGBoost will be skipped if not available
- Other models will still run successfully

## Expected Runtime

On a typical HPC node:
- Step 1 (Data prep): ~1 minute
- Step 2 (Download): Variable - see download troubleshooting
  - From login node: 2-4 hours (good network)
  - From compute node: Often fails due to network restrictions
- Step 3 (Clustering): 30-90 minutes (with 16 threads, depends on genome count)
- Step 4 (Classification): 5-15 minutes
- Step 5 (Analysis): 1-2 minutes

**Total: 3-6 hours** (assuming successful downloads)

**NOTE:** If downloads fail in the SLURM job, run Step 2 separately from
the login node before resubmitting the rest of the pipeline.

## Citation

If you use this pipeline, please cite:
- MMseqs2: Steinegger & Söding (2017) Nature Biotechnology
- scikit-learn: Pedregosa et al. (2011) JMLR
- XGBoost: Chen & Guestrin (2016) KDD

## Contact

For questions or issues, please contact the pipeline developer.

## Version History

- v1.0 (2025-11): Initial release
  - Complete pipeline with genus-level CV
  - Multiple classification models
  - Feature directionality analysis
