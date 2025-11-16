#!/bin/bash
#SBATCH --job-name=pangenome_consolidated
#SBATCH --output=pangenome_consolidated_%j.out
#SBATCH --error=pangenome_consolidated_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=quickq

# Consolidated Pangenome Analysis Pipeline SLURM Submission Script
# This script runs the complete pangenome analysis pipeline

echo "======================================================================"
echo "CONSOLIDATED PANGENOME ANALYSIS PIPELINE"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start time: $(date)"
echo "======================================================================"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Script directory: $SCRIPT_DIR"

# Change to script directory to ensure all files are accessible
cd "$SCRIPT_DIR"
echo "Working directory: $(pwd)"
echo ""

# Load required modules
echo "Loading modules..."
module purge
module load StdEnv
module load python/3.12
module load mmseqs2/15-6f452

# Load NCBI dataset tools if available
# Try different possible module names
if module avail ncbi-datasets 2>&1 | grep -q ncbi-datasets; then
    module load ncbi-datasets
    echo "Loaded NCBI datasets module"
elif module avail datasets 2>&1 | grep -q datasets; then
    module load datasets
    echo "Loaded datasets module"
elif module avail ncbi-blast 2>&1 | grep -q ncbi-blast; then
    module load ncbi-blast
    echo "Loaded ncbi-blast module (provides NCBI tools)"
elif module avail ncbi 2>&1 | grep -q ncbi; then
    module load ncbi
    echo "Loaded NCBI module"
else
    echo "WARNING: NCBI datasets module not found."
    echo "The pipeline will try to use direct HTTPS download as fallback."
    echo "Try 'module keyword ncbi' to see available NCBI-related modules on this cluster."
fi

echo ""
echo "Loaded modules:"
module list
echo ""

# Set up Python environment
echo "Setting up Python environment..."

# Create Python virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate
echo ""

# Install/upgrade required Python packages
echo "Python version:"
python3 --version
echo ""

echo "Installing Python packages..."
pip install --upgrade pip --quiet
pip install --quiet pandas numpy scipy scikit-learn matplotlib seaborn xgboost biopython

if [ $? -eq 0 ]; then
    echo "Python packages installed successfully"
else
    echo "WARNING: Some Python packages may have failed to install"
fi
echo ""

# Check for input file
echo "Checking for input file..."
if [ -f "nif_hdk_hits_enriched_with_quality_checkm.csv" ]; then
    echo "Input file found: nif_hdk_hits_enriched_with_quality_checkm.csv"
else
    echo "ERROR: Input file not found: nif_hdk_hits_enriched_with_quality_checkm.csv"
    exit 1
fi
echo ""

# Check for required scripts
echo "Checking for required scripts..."
if [ -f "pangenome_pipeline_consolidated.py" ]; then
    echo "Main script found: pangenome_pipeline_consolidated.py"
else
    echo "ERROR: Main script not found: pangenome_pipeline_consolidated.py"
    exit 1
fi

if [ -f "02_download_proteins.py" ]; then
    echo "Download script found: 02_download_proteins.py"
else
    echo "WARNING: Download script not found: 02_download_proteins.py"
    echo "Protein download step may fail."
fi
echo ""

# Check disk space
echo "Disk space:"
df -h . | head -2
echo ""

echo "======================================================================"
echo "STARTING PIPELINE EXECUTION"
echo "======================================================================"
echo ""

# Run the consolidated pipeline
# Pass number of threads from SLURM allocation
python3 pangenome_pipeline_consolidated.py \
    --input nif_hdk_hits_enriched_with_quality_checkm.csv \
    --threads $SLURM_CPUS_PER_TASK \
    --non-interactive \
    --steps 1-8

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "======================================================================"
    echo "PIPELINE COMPLETED SUCCESSFULLY"
    echo "======================================================================"
else
    echo "======================================================================"
    echo "PIPELINE FAILED WITH EXIT CODE: $EXIT_CODE"
    echo "======================================================================"
fi
echo "End time: $(date)"
echo ""

echo ""
echo "======================================================================"
echo "JOB COMPLETE"
echo "======================================================================"

exit $EXIT_CODE
