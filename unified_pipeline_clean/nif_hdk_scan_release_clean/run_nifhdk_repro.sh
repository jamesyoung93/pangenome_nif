#!/usr/bin/env bash
# Resumable driver for nifH/D/K scan + robust summarize (+ optional metadata)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"; cd "$ROOT"
STAMP_DIR="results/.stamps"; mkdir -p "$STAMP_DIR" logs results
SUMMARIZER="scripts/summarize_robust_fixcase.py"

subset=""; force=""; with_meta=0; assembly_levels=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --subset) subset="$2"; shift 2;;
    --force)  force="$2"; shift 2;;
    --assembly-levels) assembly_levels="$2"; shift 2;;
    --with-metadata|--with_meta|--metadata) with_meta=1; shift;;
    *) echo "Unknown arg: $1"; exit 2;;
  esac
done

if [[ -n "$assembly_levels" ]]; then
  export NIF_ASSEMBLY_LEVELS="$assembly_levels"
fi

have(){ test -s "$1"; }
touch_ok(){ mkdir -p "$(dirname "$1")"; : > "$1"; }

# cascading clear
clear_from(){
  case "$1" in
    query)      rm -f "$STAMP_DIR"/{query,download,combine,scan,summarize,metadata}.ok;;
    download)   rm -f "$STAMP_DIR"/{download,combine,scan,summarize,metadata}.ok;;
    combine)    rm -f "$STAMP_DIR"/{combine,scan,summarize,metadata}.ok;;
    scan)       rm -f "$STAMP_DIR"/{scan,summarize,metadata}.ok;;
    summarize)  rm -f "$STAMP_DIR"/{summarize,metadata}.ok;;
    metadata)   rm -f "$STAMP_DIR"/metadata.ok;;
    all)        rm -f "$STAMP_DIR"/*.ok 2>/dev/null || true;;
  esac
}
[[ -n "${force:-}" ]] && clear_from "$force"

# Required for NCBI
: "${ENTREZ_EMAIL:?Please export ENTREZ_EMAIL (e.g. export ENTREZ_EMAIL=you@org)}"

step_query(){
  echo "[query] fetching assemblies…"
  python pipeline.py query |& tee -a logs/pipeline.query.log
  if [[ -n "$subset" ]]; then
    awk -v n="$subset" 'NR==1 || NR<=n+1' results/assemblies/assemblies.tsv > results/assemblies/assemblies.tmp
    mv results/assemblies/assemblies.tmp results/assemblies/assemblies.tsv
    echo "[subset] keeping first $subset assemblies"
  fi
  touch_ok "$STAMP_DIR/query.ok"
}
step_download(){
  echo "[download] downloading proteomes (resumable)…"
  python pipeline.py download |& tee -a logs/pipeline.download.log
  dl_count=$(find results/downloads -type f -name '*_protein.faa*' | wc -l | awk '{print $1}')
  [[ "$dl_count" -gt 0 ]] || { echo "ERROR: no proteomes found after download"; exit 1; }
  touch_ok "$STAMP_DIR/download.ok"
}
step_combine(){
  echo "[combine] building combined_proteins.faa…"
  python pipeline.py combine |& tee -a logs/pipeline.combine.log
  have results/proteomes/combined_proteins.faa || { echo "ERROR: combine failed"; exit 1; }
  touch_ok "$STAMP_DIR/combine.ok"
}
step_scan(){
  echo "[scan] running hmmsearch scans…"
  python pipeline.py scan |& tee -a logs/pipeline.scan.log
  # harmonize directory name if pipeline wrote to hmmsearch
  if [[ ! -d results/hmmer && -d results/hmmsearch ]]; then ln -s hmmsearch results/hmmer; fi
  hit_tbls=$(find results/hmmer -type f -name '*.domtblout' 2>/dev/null | wc -l | awk '{print $1}')
  [[ "$hit_tbls" -gt 0 ]] || echo "WARNING: no *.domtblout found under results/hmmer"
  touch_ok "$STAMP_DIR/scan.ok"
}
step_summarize(){
  echo "[summarize] robust summarization…"
  domdir="results/hmmer"; [[ -d "$domdir" ]] || domdir="results/hmmsearch"
  python "$SUMMARIZER" \
    --domdir "$domdir" \
    --downloads results/downloads \
    --combined results/proteomes/combined_proteins.faa \
    --assemblies results/assemblies/assemblies.tsv \
    --out results/nif_hdk_hits.csv \
    --thr 1e-5 |& tee -a logs/pipeline.summarize.log
  touch_ok "$STAMP_DIR/summarize.ok"
}

ensure_datasets_cli(){
  if command -v datasets >/dev/null 2>&1; then
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    echo "[metadata] installing NCBI 'datasets' CLI via pip…" >&2
    PIP_DISABLE_PIP_VERSION_CHECK=1 python3 -m pip install --quiet --user ncbi-datasets-pylib || true
    # refresh PATH for user installs if needed
    PATH="$HOME/.local/bin:$PATH"
    if command -v datasets >/dev/null 2>&1; then
      echo "[metadata] 'datasets' CLI available after pip install" >&2
      return 0
    fi
  fi
  return 1
}
step_metadata(){
  [[ "$with_meta" -eq 1 ]] || return 0
  if ! ensure_datasets_cli; then
    echo "ERROR: required NCBI 'datasets' CLI not found and pip install failed." >&2
    echo "Install via 'mamba install -c conda-forge ncbi-datasets-cli' (or 'conda install ...')," >&2
    echo "load an HPC module such as 'module load ncbi-datasets', or see" >&2
    echo "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/." >&2
    exit 1
  fi
  echo "[metadata] fetching NCBI assembly metadata…"
  mkdir -p results/assemblies
  cut -f1 results/assemblies/assemblies.tsv | sed '1d' | sort -u > results/assemblies/asm_list.txt
  datasets summary genome accession --inputfile results/assemblies/asm_list.txt --as-json-lines > results/assemblies/assembly_data_report.jsonl
  python scripts/extract_ncbi_assembly_quality.py --in results/assemblies/assembly_data_report.jsonl --out results/assemblies/assembly_quality.tsv
  python scripts/enrich_and_join.py --hits results/nif_hdk_hits.csv --quality results/assemblies/assembly_quality.tsv --out results/nif_hdk_hits.enriched_with_quality_checkm.csv
  touch_ok "$STAMP_DIR/metadata.ok"
}

# Orchestrate with stamps
[[ "$force" == "query"      ]] && clear_from query
[[ "$force" == "download"   ]] && clear_from download
[[ "$force" == "combine"    ]] && clear_from combine
[[ "$force" == "scan"       ]] && clear_from scan
[[ "$force" == "summarize"  ]] && clear_from summarize
[[ "$force" == "metadata"   ]] && clear_from metadata

if ! have "$STAMP_DIR/query.ok";      then step_query;     else echo "[query] already done."; fi
if ! have "$STAMP_DIR/download.ok";   then step_download;  else echo "[download] already done."; fi
if ! have "$STAMP_DIR/combine.ok";    then step_combine;   else echo "[combine] already done."; fi
if ! have "$STAMP_DIR/scan.ok";       then step_scan;      else echo "[scan] already done."; fi
if ! have "$STAMP_DIR/summarize.ok";  then step_summarize; else echo "[summarize] already done."; fi
if ! have "$STAMP_DIR/metadata.ok";   then step_metadata;  else echo "[metadata] already done."; fi

echo "All done. Output: results/nif_hdk_hits.csv"
