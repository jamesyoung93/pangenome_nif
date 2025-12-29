#!/usr/bin/env bash
set -euo pipefail

# ---- config / usage ---------------------------------------------------------
usage() {
  cat <<'USAGE'
Usage:
  ./setup_refs.sh --source /path/to/NFixDB_or_other --commit [--build-from-fasta]

What it does:
  - Finds *.hmm under --source and copies ONLY nifH/nifD/nifK profiles into:
      refs/hmms/nifH/
      refs/hmms/nifD/
      refs/hmms/nifK/
  - Dry-run by default; use --commit to actually copy.
  - If a subunit has no HMMs and --build-from-fasta is given, it looks for FASTA
    seeds (*.fa|*.fasta) and runs hmmbuild to create HMMs in the right folder.

Examples:
  # (A) copy from your diazocyano_v2 NFixDB:
  ./setup_refs.sh --source ~/projects/diazocyano_v2/data/hmms/NFixDB

  # (B) actually copy (not dry-run) and build from FASTAs if needed:
  ./setup_refs.sh --source ~/projects/diazocyano_v2/data/hmms/NFixDB --commit --build-from-fasta

  # (C) if you keep seeds elsewhere:
  ./setup_refs.sh --source /path/to/seeds_or_hmms --commit --build-from-fasta
USAGE
}

SRC=""
DO_COMMIT=0
DO_BUILD=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) SRC="${2:-}"; shift 2;;
    --commit) DO_COMMIT=1; shift;;
    --build-from-fasta) DO_BUILD=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

if [[ -z "${SRC}" ]]; then
  echo "ERROR: --source is required"; usage; exit 1
fi
if [[ ! -d "${SRC}" ]]; then
  echo "ERROR: --source not found: ${SRC}" >&2; exit 1
fi

# where to place things (inside nif_hdk_scan):
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_H= "${ROOT}/refs/hmms/nifH"
DEST_D= "${ROOT}/refs/hmms/nifD"
DEST_K= "${ROOT}/refs/hmms/nifK"
mkdir -p "${ROOT}/refs/hmms/nifH" "${ROOT}/refs/hmms/nifD" "${ROOT}/refs/hmms/nifK"

echo "Source: ${SRC}"
echo "Dest  : ${ROOT}/refs/hmms/{nifH,nifD,nifK}"
echo "Mode  : $([[ ${DO_COMMIT} -eq 1 ]] && echo COMMIT || echo DRY-RUN)"
echo

# ---- helper to classify filename -> subunit -------------------------------
classify() {
  local f="$(basename "$1")"
  local lf="$(echo "$f" | tr '[:upper:]' '[:lower:]')"
  # strictly nif(H|D|K); ignore vnf/anf
  if   echo "$lf" | grep -Eq '(?:^|[^a-z])nifh([^a-z]|$)'; then echo "nifH"
  elif echo "$lf" | grep -Eq '(?:^|[^a-z])nifd([^a-z]|$)'; then echo "nifD"
  elif echo "$lf" | grep -Eq '(?:^|[^a-z])nifk([^a-z]|$)'; then echo "nifK"
  else echo "skip"
  fi
}

# ---- collect HMMs under source --------------------------------------------
mapfile -t HMM_LIST < <(find "${SRC}" -type f -iname '*.hmm' 2>/dev/null | sort || true)
COUNT_H=0; COUNT_D=0; COUNT_K=0

if [[ ${#HMM_LIST[@]} -gt 0 ]]; then
  echo "Discovered $((${#HMM_LIST[@]})) HMM files under source; classifying..."
  for hmm in "${HMM_LIST[@]}"; do
    sub="$(classify "$hmm")"
    case "$sub" in
      nifH) TARGET="${DEST_H}"; COUNT_H=$((COUNT_H+1));;
      nifD) TARGET="${DEST_D}"; COUNT_D=$((COUNT_D+1));;
      nifK) TARGET="${DEST_K}"; COUNT_K=$((COUNT_K+1));;
      *)    continue;;
    esac
    if [[ ${DO_COMMIT} -eq 1 ]]; then
      cp -v "$hmm" "$TARGET/"
    else
      echo "[DRY-RUN] would copy: $hmm  ->  $TARGET/"
    fi
  done
else
  echo "No *.hmm files found under --source (will try FASTA if --build-from-fasta)."
fi

echo
echo "HMMs (nifH/nifD/nifK) discovered: $COUNT_H / $COUNT_D / $COUNT_K"
echo

# ---- optionally build from FASTA if missing any subunit -------------------
if [[ ${DO_BUILD} -eq 1 ]]; then
  # does hmmbuild exist?
  if ! command -v hmmbuild >/dev/null 2>&1; then
    echo "WARNING: hmmbuild not found on PATH; cannot build from FASTA." >&2
  else
    build_for() {
      local sub="$1"
      local dest="${ROOT}/refs/hmms/${sub}"
      local built=0
      # look for FASTAs under source (and also under nif_hdk_scan/refs/fasta/sub)
      mapfile -t FA_LIST < <(
        find "${SRC}" -type f \( -iname '*.fa' -o -iname '*.fasta' -o -iname '*.faa' \) -printf '%p\n' 2>/dev/null | grep -i "/${sub}/" || true
      )
      # include local seed folder if present
      if [[ -d "${ROOT}/refs/fasta/${sub}" ]]; then
        while IFS= read -r p; do FA_LIST+=("$p"); done < <(find "${ROOT}/refs/fasta/${sub}" -type f -iname '*.fa*' -printf '%p\n' 2>/dev/null || true)
      fi
      if [[ ${#FA_LIST[@]} -eq 0 ]]; then
        echo "No ${sub} seed FASTAs found for hmmbuild."
        return
      fi
      for fa in "${FA_LIST[@]}"; do
        bn="$(basename "${fa}")"
        hmm_out="${dest}/${bn%.*}.hmm"
        if [[ ${DO_COMMIT} -eq 1 ]]; then
          echo "hmmbuild ${hmm_out} ${fa}"
          hmmbuild "${hmm_out}" "${fa}" >/dev/null
        else
          echo "[DRY-RUN] would hmmbuild ${hmm_out} from ${fa}"
        fi
        built=$((built+1))
      done
      echo "Built ${built} HMM(s) for ${sub}"
    }

    if [[ $COUNT_H -eq 0 ]]; then build_for "nifH"; fi
    if [[ $COUNT_D -eq 0 ]]; then build_for "nifD"; fi
    if [[ $COUNT_K -eq 0 ]]; then build_for "nifK"; fi
  fi
fi

echo
echo "Done. Next steps:"
echo "  1) python pipeline.py scan"
echo "  2) python pipeline.py summarize"
