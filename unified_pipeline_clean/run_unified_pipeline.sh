#!/usr/bin/env bash
set -euo pipefail

# --- User-tunable settings (edit as needed) ---------------------------------
ENTREZ_EMAIL="${ENTREZ_EMAIL:-you@example.org}"   # required for NCBI downloads
UPSTREAM_DIR="nif_hdk_scan_release_clean"          # path to nifH/D/K scan pipeline
UPSTREAM_SUBSET=""                                 # e.g., 200 to smoke test
UPSTREAM_WITH_METADATA=1                            # 1 to fetch metadata/checkM
UPSTREAM_FORCE=""                                   # e.g., query|download|scan|all
UPSTREAM_ASSEMBLY_LEVELS="complete genome"          # comma-separated; empty = config defaults
HMMSEARCH_CMD="${HMMSEARCH_CMD:-hmmsearch}"        # override to a fully qualified hmmsearch path
HMMER_MODULE="${HMMER_MODULE:-}"                    # optional: module name to auto-load if hmmsearch is missing

RUN_DIR="unified_pipeline_run"                     # downstream working dir
DOWNSTREAM_OUTPUT_ROOT="results"                   # saved under ${RUN_DIR}/
DOWNSTREAM_THREADS=16                               # threads for clustering
DOWNSTREAM_MIN_GENOMES=40                           # min genomes per family
DOWNSTREAM_MMSEQS_ID=0.80                           # MMseqs2 identity threshold
DOWNSTREAM_MMSEQS_COV=0.80                          # MMseqs2 coverage threshold
DOWNSTREAM_PURITY=0.90                              # purity filter for CV
DOWNSTREAM_CV_GROUP="genus"                        # genus|family|order
DOWNSTREAM_CV_FOLDS=5
DOWNSTREAM_CV_SEED=42
DOWNSTREAM_TOP_N=100
DOWNSTREAM_SKIP_DOWNLOAD=0                          # set 1 if proteins already present

DOWNSTREAM_MODE="full"                              # full|experiment|comparative_experiment_matrix
DOWNSTREAM_EXPERIMENT_MODEL="xgboost"              # only used in experiment mode
# ----------------------------------------------------------------------------

ROOT="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${ROOT}/.." && pwd)"
UPSTREAM_ROOT="${ROOT}/${UPSTREAM_DIR}"
DOWNSTREAM_SCRIPT="${REPO_ROOT}/nif_downstream_code/pangenome_pipeline_consolidated2.py"

[[ -x "${UPSTREAM_ROOT}/run_nifhdk_repro.sh" ]] || { echo "Missing upstream driver at ${UPSTREAM_ROOT}"; exit 1; }
[[ -f "${DOWNSTREAM_SCRIPT}" ]] || { echo "Missing downstream script at ${DOWNSTREAM_SCRIPT}"; exit 1; }

python - <<'PY'
import importlib, sys
deps = {
    "yaml": "PyYAML",
    "Bio": "biopython",
    "pandas": "pandas",
    "tqdm": "tqdm",
    "numpy": "numpy",
    "sklearn": "scikit-learn",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "xgboost": "xgboost",
}
missing = []
for mod, pkg in deps.items():
    try:
        importlib.import_module(mod)
    except ImportError:
        missing.append(pkg)

if missing:
    print("Missing Python packages: " + ", ".join(sorted(set(missing))), file=sys.stderr)
    print("Install them with: python -m pip install -r requirements.txt", file=sys.stderr)
    sys.exit(1)
PY

# Try to make hmmsearch available (honor HMMSEARCH_CMD and optional module)
if [[ "${HMMSEARCH_CMD}" == */* && -x "${HMMSEARCH_CMD}" ]]; then
  PATH="$(dirname "${HMMSEARCH_CMD}"):${PATH}"
fi

if ! command -v "${HMMSEARCH_CMD}" >/dev/null 2>&1; then
  if [[ -n "${HMMER_MODULE}" ]] && command -v module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    module load "${HMMER_MODULE}" >/dev/null 2>&1 || true
  fi
fi

# Fail fast if HMMER is still missing before starting downloads
if ! command -v "${HMMSEARCH_CMD}" >/dev/null 2>&1; then
  cat >&2 <<'EOF'
ERROR: hmmsearch (HMMER) is not on your PATH.
Either:
  * export HMMSEARCH_CMD=/path/to/hmmsearch, or
  * set HMMER_MODULE (e.g., hmmer/3.4) so the script can `module load` it, or
  * manually load your cluster's HMMER module before running.
EOF
  exit 1
fi

export HMMSEARCH_CMD

export ENTREZ_EMAIL

pushd "${UPSTREAM_ROOT}" >/dev/null
UP_ARGS=()
[[ -n "${UPSTREAM_SUBSET}" ]] && UP_ARGS+=(--subset "${UPSTREAM_SUBSET}")
[[ -n "${UPSTREAM_FORCE}" ]] && UP_ARGS+=(--force "${UPSTREAM_FORCE}")
[[ -n "${UPSTREAM_ASSEMBLY_LEVELS}" ]] && UP_ARGS+=(--assembly-levels "${UPSTREAM_ASSEMBLY_LEVELS}")
[[ ${UPSTREAM_WITH_METADATA} -eq 1 ]] && UP_ARGS+=(--with-metadata)
bash run_nifhdk_repro.sh "${UP_ARGS[@]:+${UP_ARGS[@]}}"
popd >/dev/null

UPSTREAM_OUT="${UPSTREAM_ROOT}/results/nif_hdk_hits.csv"
if [[ ${UPSTREAM_WITH_METADATA} -eq 1 ]]; then
  UPSTREAM_OUT="${UPSTREAM_ROOT}/results/nif_hdk_hits.enriched_with_quality_checkm.csv"
fi
[[ -s "${UPSTREAM_OUT}" ]] || { echo "Upstream output missing: ${UPSTREAM_OUT}"; exit 1; }

mkdir -p "${ROOT}/${RUN_DIR}"
cp -f "${UPSTREAM_OUT}" "${ROOT}/${RUN_DIR}/nif_hdk_hits_enriched_with_quality_checkm.csv"

# stage assembly quality metadata for filtering
mkdir -p "${ROOT}/${RUN_DIR}/${DOWNSTREAM_OUTPUT_ROOT}/assemblies"
if [[ -f "${UPSTREAM_ROOT}/results/assemblies/assembly_quality.tsv" ]]; then
  cp -f "${UPSTREAM_ROOT}/results/assemblies/assembly_quality.tsv" \
    "${ROOT}/${RUN_DIR}/${DOWNSTREAM_OUTPUT_ROOT}/assemblies/assembly_quality.tsv"
else
  echo "ERROR: Expected assembly quality table at ${UPSTREAM_ROOT}/results/assemblies/assembly_quality.tsv" >&2
  exit 1
fi

pushd "${ROOT}/${RUN_DIR}" >/dev/null

DOWN_ARGS=(
  --mode "${DOWNSTREAM_MODE}"
  --input nif_hdk_hits_enriched_with_quality_checkm.csv
  --threads "${DOWNSTREAM_THREADS}"
  --min-genomes "${DOWNSTREAM_MIN_GENOMES}"
  --mmseqs-identity "${DOWNSTREAM_MMSEQS_ID}"
  --mmseqs-coverage "${DOWNSTREAM_MMSEQS_COV}"
  --purity-threshold "${DOWNSTREAM_PURITY}"
  --cv-group-col "${DOWNSTREAM_CV_GROUP}"
  --cv-folds "${DOWNSTREAM_CV_FOLDS}"
  --cv-seed "${DOWNSTREAM_CV_SEED}"
  --top-n-features "${DOWNSTREAM_TOP_N}"
  --output-root "${DOWNSTREAM_OUTPUT_ROOT}"
)
[[ ${DOWNSTREAM_SKIP_DOWNLOAD} -eq 1 ]] && DOWN_ARGS+=(--skip-download)
if [[ "${DOWNSTREAM_MODE}" == "comparative_experiment_matrix" ]]; then
  DOWN_ARGS+=(--experiment-model "${DOWNSTREAM_EXPERIMENT_MODEL}")
fi

python "${DOWNSTREAM_SCRIPT}" "${DOWN_ARGS[@]}"

popd >/dev/null

echo "Unified pipeline complete. Outputs are under ${RUN_DIR}/${DOWNSTREAM_OUTPUT_ROOT}."
