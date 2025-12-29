#!/usr/bin/env bash
set -euo pipefail

# --- User-tunable settings (edit as needed) ---------------------------------
ENTREZ_EMAIL="${ENTREZ_EMAIL:-you@example.org}"   # required for NCBI downloads
UPSTREAM_DIR="nif_hdk_scan_release_clean"          # path to nifH/D/K scan pipeline
UPSTREAM_SUBSET=""                                 # e.g., 200 to smoke test
UPSTREAM_WITH_METADATA=1                            # 1 to fetch metadata/checkM
UPSTREAM_FORCE=""                                   # e.g., query|download|scan|all
UPSTREAM_ASSEMBLY_LEVELS="complete genome"          # comma-separated; empty = config defaults

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

DOWNSTREAM_MODE="pipeline"                          # pipeline|comparative_experiment_matrix
DOWNSTREAM_EXPERIMENT_MODEL="xgboost"              # only used in experiment mode
# ----------------------------------------------------------------------------

ROOT="$(cd "$(dirname "$0")" && pwd)"
UPSTREAM_ROOT="${ROOT}/${UPSTREAM_DIR}"
DOWNSTREAM_SCRIPT="${ROOT}/nif_downstream_code/pangenome_pipeline_consolidated2.py"

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
