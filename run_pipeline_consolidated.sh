#!/bin/bash
#SBATCH --job-name=diazotroph_consolidated
#SBATCH --output=diazotroph_consolidated_%j.out
#SBATCH --error=diazotroph_consolidated_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=bigmem

# Diazotroph Classification Pipeline (Consolidated)
# Runs the full pipeline via the single Python driver (pangenome_pipeline_consolidated.py).

set -euo pipefail

echo "======================================================================"
echo "DIAZOTROPH CONSOLIDATED PIPELINE"
echo "======================================================================"
echo "Job ID: ${SLURM_JOB_ID:-N/A}"
echo "Node: ${SLURM_NODELIST:-N/A}"
echo "Partition: ${SLURM_JOB_PARTITION:-N/A}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Memory: ${SLURM_MEM_PER_NODE:-N/A}MB"
echo "Start time: $(date)"
echo "======================================================================"
echo ""

# Resolve script directory to ensure relative paths work even when submitted elsewhere
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Working directory: $(pwd)"

echo "======================================================================"
echo "Loading modules..."
echo "======================================================================"
module purge
module load StdEnv
module load mmseqs2/15-6f452

# Optional conda environment activation (set CONDA_ENV_NAME before submitting)
if [ -n "${CONDA_ENV_NAME:-}" ]; then
    echo "Activating conda environment: ${CONDA_ENV_NAME}"
    if ! command -v conda >/dev/null 2>&1; then
        # Attempt to source a common conda profile script if the shell does not know about conda yet
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

# Optional NCBI datasets module (used by the consolidated script if available)
if module avail ncbi-datasets 2>&1 | grep -q ncbi-datasets; then
    module load ncbi-datasets
    echo "Loaded ncbi-datasets module"
elif module avail datasets 2>&1 | grep -q datasets; then
    module load datasets
    echo "Loaded datasets module"
elif module avail ncbi-blast 2>&1 | grep -q ncbi-blast; then
    module load ncbi-blast
    echo "Loaded ncbi-blast module (provides NCBI tools)"
elif module avail ncbi 2>&1 | grep -q ncbi; then
    module load ncbi
    echo "Loaded ncbi module"
else
    echo "NCBI datasets module not found; HTTPS downloads will be used if needed"
    echo "Run 'module keyword ncbi' to see available NCBI-related modules on this cluster"
fi

module list

echo "Using pre-activated Python environment (venv or conda)"
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
    "biopython",
]
missing = [pkg for pkg in required if importlib.util.find_spec(pkg) is None]
if missing:
    print("ERROR: Missing Python packages: " + ", ".join(missing))
    print("Install them in your environment before submitting this job.")
    sys.exit(1)
print("All required Python packages found.")
PY

echo "======================================================================"
echo "Validating inputs..."
echo "======================================================================"
INPUT_FILE="${1:-nif_hdk_hits_enriched_with_quality_checkm.csv}"
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

if [ ! -f "pangenome_pipeline_consolidated.py" ]; then
    echo "ERROR: pangenome_pipeline_consolidated.py not found in $(pwd)"
    exit 1
fi

# Show disk space for awareness
printf "\nDisk space:\n"
df -h . | head -2

echo "\n======================================================================"
echo "STARTING CONSOLIDATED PIPELINE"
echo "======================================================================"
python3 pangenome_pipeline_consolidated.py \
    --input "$INPUT_FILE" \
    --threads "${SLURM_CPUS_PER_TASK:-32}" \
    --non-interactive \
    --steps 1-8

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "\n======================================================================"
    echo "PIPELINE COMPLETED SUCCESSFULLY"
    echo "======================================================================"
else
    echo "\n======================================================================"
    echo "PIPELINE FAILED WITH EXIT CODE: $EXIT_CODE"
    echo "======================================================================"
fi

echo "End time: $(date)"
exit $EXIT_CODE
