#!/bin/bash
#SBATCH --job-name=diazotroph_classification
#SBATCH --output=diazotroph_classification_%j.out
#SBATCH --error=diazotroph_classification_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=bigmem


# Diazotroph Classification Pipeline
# This script runs the complete analysis pipeline for classifying
# cyanobacteria as diazotrophs or non-diazotrophs based on gene content

echo "======================================================================"
echo "DIAZOTROPH CLASSIFICATION PIPELINE"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "======================================================================"

# Load required modules
module purge
module load StdEnv
module load python/3.12
module load mmseqs2/15-6f452

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Install required Python packages
echo ""
echo "Installing/updating Python packages..."
pip install --upgrade pip
pip install pandas numpy scipy scikit-learn matplotlib seaborn xgboost

echo ""
echo "======================================================================"
echo "STEP 1: Data Preparation"
echo "======================================================================"
python3 01_prepare_data.py nif_hdk_hits_enriched_with_quality_checkm.csv

if [ $? -ne 0 ]; then
    echo "Error in Step 1"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 2: Download Protein Sequences"
echo "======================================================================"
echo "Note: This step may take several hours..."
python3 02_download_proteins.py

if [ $? -ne 0 ]; then
    echo "Error in Step 2"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 3: Create Gene Families"
echo "======================================================================"
echo "Note: This step may take 1-2 hours..."
python3 03_create_gene_families.py $SLURM_CPUS_PER_TASK

if [ $? -ne 0 ]; then
    echo "Error in Step 3"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 4: Classification with Genus Cross-Validation"
echo "======================================================================"
python3 04_classify.py

if [ $? -ne 0 ]; then
    echo "Error in Step 4"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 5: Feature Directionality Analysis"
echo "======================================================================"
python3 05_analyze_feature_directionality.py

if [ $? -ne 0 ]; then
    echo "Error in Step 5"
    exit 1
fi

echo ""
echo "======================================================================"
echo "CREATING RESULTS ARCHIVE"
echo "======================================================================"

# Create results directory
mkdir -p results

# Copy all important outputs
cp complete_genomes_labeled.csv results/
cp gene_family_matrix.csv results/
cp gene_family_info.csv results/
cp gene_family_purity_stats.csv results/
cp classification_summary.csv results/
cp roc_curves.png results/
cp feature_importance_*.csv results/
cp feature_importance_plot.png results/
cp feature_directionality_*.csv results/
cp feature_effects.png results/
cp fold_metrics_*.csv results/
cp *.out results/ 2>/dev/null
cp *.err results/ 2>/dev/null

# Create summary report
cat > results/ANALYSIS_SUMMARY.txt << EOF
====================================================================
DIAZOTROPH CLASSIFICATION ANALYSIS SUMMARY
====================================================================
Job ID: $SLURM_JOB_ID
Completion time: $(date)

PIPELINE STEPS COMPLETED:
1. Data preparation and labeling ✓
2. Protein sequence download ✓
3. Gene family clustering (MMseqs2, 80% identity) ✓
4. Classification with genus-level cross-validation ✓
5. Feature directionality analysis ✓

KEY OUTPUT FILES:
- complete_genomes_labeled.csv: Labeled genome dataset
- gene_family_matrix.csv: Presence/absence matrix
- classification_summary.csv: Model performance metrics
- roc_curves.png: ROC curves for all models
- feature_importance_*.csv: Feature importance per model
- feature_directionality_full.csv: Complete feature analysis
- feature_effects.png: Visualization of feature effects

See README.txt for detailed description of all outputs.
====================================================================
EOF

# Create archive
echo "Creating results.tar.gz..."
tar -czf results.tar.gz results/

echo ""
echo "======================================================================"
echo "PIPELINE COMPLETE!"
echo "======================================================================"
echo "End time: $(date)"
echo ""
echo "Results saved in: results.tar.gz"
echo "Extract with: tar -xzf results.tar.gz"
echo ""
echo "To view summary: cat results/ANALYSIS_SUMMARY.txt"
echo "======================================================================"
