#!/bin/bash
#SBATCH --job-name=pangenome_consolidated
#SBATCH --output=pangenome_consolidated_%j.out
#SBATCH --error=pangenome_consolidated_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=bigmem

# ============================================================================
# Consolidated Pangenome Analysis Pipeline - SLURM Submission Script
# For SDSU Innovator HPC (bigmem partition)
# ============================================================================

echo "======================================================================"
echo "CONSOLIDATED PANGENOME ANALYSIS PIPELINE"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 256GB"
echo "Start time: $(date)"
echo "======================================================================"

# Load required modules (Python environment must already be active)
echo ""
echo "Loading modules..."
module purge
module load StdEnv
module load mmseqs2/15-6f452

# Show loaded modules
echo ""
echo "Loaded modules:"
module list

echo ""
echo "Using pre-activated Python environment (venv or conda)"

# Optional conda environment activation (set CONDA_ENV_NAME before submitting)
if [ -n "${CONDA_ENV_NAME:-}" ]; then
    echo "Activating conda environment: ${CONDA_ENV_NAME}"
    if ! command -v conda >/dev/null 2>&1; then
        if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
            . "$HOME/miniconda3/etc/profile.d/conda.sh"
        elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
            . "$HOME/anaconda3/etc/profile.d/conda.sh"
        elif [ -f "$HOME/.conda/etc/profile.d/conda.sh" ]; then
            . "$HOME/.conda/etc/profile.d/conda.sh"
        else
            echo "ERROR: conda command not found. Load your conda module or add conda.sh to your PATH."
            exit 1
        fi
    fi

    if ! conda activate "$CONDA_ENV_NAME"; then
        echo "ERROR: Failed to activate conda environment '$CONDA_ENV_NAME'." \
             "Ensure the environment exists on the compute node."
        exit 1
    fi
fi

python3 --version

echo "Checking required Python packages..."
python3 - <<'PY'
import importlib.util, sys
required = [
    "pandas",
    "numpy",
    "scipy",
    "sklearn",
    "matplotlib",
    "seaborn",
    "xgboost",
]
missing = [pkg for pkg in required if importlib.util.find_spec(pkg) is None]
if missing:
    print("ERROR: Missing Python packages: " + ", ".join(missing))
    print("Install them in your environment before submitting this job.")
    sys.exit(1)
print("All required Python packages found.")
PY

# Check for required input file
echo ""
echo "Checking for input file..."
INPUT_FILE="nif_hdk_hits_enriched_with_quality_checkm.csv"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "Please ensure the input CSV file is in the current directory"
    exit 1
fi
echo "Input file found: $INPUT_FILE"

# Check for required scripts
echo ""
echo "Checking for required scripts..."
MAIN_SCRIPT="pangenome_pipeline_consolidated.py"
if [ ! -f "$MAIN_SCRIPT" ]; then
    echo "ERROR: Main pipeline script not found: $MAIN_SCRIPT"
    exit 1
fi
echo "Main script found: $MAIN_SCRIPT"

# Make script executable
chmod +x "$MAIN_SCRIPT"

# Display disk space
echo ""
echo "Disk space:"
df -h .

# Run the consolidated pipeline
echo ""
echo "======================================================================"
echo "STARTING PIPELINE EXECUTION"
echo "======================================================================"
echo ""

# Run with all CPUs available
python3 "$MAIN_SCRIPT" \
    --input "$INPUT_FILE" \
    --threads $SLURM_CPUS_PER_TASK \
    --min-genomes 40 \
    --non-interactive

PIPELINE_EXIT_CODE=$?

echo ""
echo "======================================================================"
if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo "PIPELINE COMPLETED SUCCESSFULLY!"
else
    echo "PIPELINE FAILED WITH EXIT CODE: $PIPELINE_EXIT_CODE"
fi
echo "======================================================================"
echo "End time: $(date)"
echo ""

# Create results summary
if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo "Creating results archive..."

    # Create results directory
    mkdir -p results_$(date +%Y%m%d)

    # Copy all important outputs
    cp complete_genomes_labeled.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp complete_genomes_with_proteins.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp gene_family_matrix.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp gene_family_info.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp gene_family_purity_stats.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp classification_summary.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp roc_curves.png results_$(date +%Y%m%d)/ 2>/dev/null
    cp feature_importance_*.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp feature_importance_plot.png results_$(date +%Y%m%d)/ 2>/dev/null
    cp feature_directionality_*.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp feature_effects.png results_$(date +%Y%m%d)/ 2>/dev/null
    cp fold_metrics_*.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp narrative_top_features.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp llm_context_features.csv results_$(date +%Y%m%d)/ 2>/dev/null
    cp gf_*.csv results_$(date +%Y%m%d)/ 2>/dev/null

    # Copy log files
    cp pangenome_consolidated_${SLURM_JOB_ID}.out results_$(date +%Y%m%d)/ 2>/dev/null
    cp pangenome_consolidated_${SLURM_JOB_ID}.err results_$(date +%Y%m%d)/ 2>/dev/null

    # Create summary report
    cat > results_$(date +%Y%m%d)/SUMMARY.txt << EOF
====================================================================
PANGENOME ANALYSIS SUMMARY
====================================================================
Job ID: $SLURM_JOB_ID
Node: $SLURM_NODELIST
Partition: $SLURM_JOB_PARTITION
CPUs: $SLURM_CPUS_PER_TASK
Memory: 256GB
Start time: $(date)
Completion time: $(date)

PIPELINE STEPS COMPLETED:
1. Data preparation and genome labeling ✓
2. Protein sequence download ✓
3. Gene family clustering (MMseqs2, 80% identity) ✓
4. Machine learning classification with genus-level CV ✓
5. Feature directionality analysis ✓
6. Annotation merging (optional) ✓
7. Gene family expansion (optional) ✓
8. Narrative table generation ✓

KEY OUTPUT FILES:
- complete_genomes_labeled.csv: Labeled genome dataset
- gene_family_matrix.csv: Gene family presence/absence matrix
- classification_summary.csv: Model performance metrics
- roc_curves.png: ROC curves for all models
- feature_importance_*.csv: Feature importance per model
- feature_directionality_full.csv: Complete feature analysis
- feature_effects.png: Visualization of feature effects
- narrative_top_features.csv: Top features with annotations

For detailed information about outputs, see the main README.
====================================================================
EOF

    # Create tar archive
    echo "Creating compressed archive..."
    tar -czf results_$(date +%Y%m%d).tar.gz results_$(date +%Y%m%d)/
    echo ""
    echo "Results saved in: results_$(date +%Y%m%d).tar.gz"
    echo "Extract with: tar -xzf results_$(date +%Y%m%d).tar.gz"
fi

echo ""
echo "======================================================================"
echo "JOB COMPLETE"
echo "======================================================================"

exit $PIPELINE_EXIT_CODE
