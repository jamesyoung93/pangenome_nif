# Unified nifH/D/K Pipeline (Minimal Bundle)

This folder is a trimmed copy of the unified pipeline so you can seed a fresh repository without the extra scratch files from the parent project. It contains only the pieces needed to run the end-to-end workflow:

- `run_unified_pipeline.sh` – single entry point that drives the upstream nifH/D/K scan and the downstream pangenome analysis.
- `nif_hdk_scan_release_clean/` – upstream scanner (HMM search + optional metadata enrichment).
- `nif_downstream_code/` – consolidated downstream pipeline and the helper scripts it expects.

## Quickstart

```bash
# start in the directory that contains this folder
cd unified_pipeline_clean

# ensure dependencies (bash, python3, mmseqs2, hmmer, curl)
# install Python packages once (includes PyYAML + Biopython needed upstream)
python -m pip install -r requirements.txt
# set your NCBI contact email for downloads
export ENTREZ_EMAIL="you@example.org"

# run the unified workflow (defaults: full upstream scan + downstream pipeline mode)
./run_unified_pipeline.sh
```

Outputs land in `unified_pipeline_run/results/`.

## Customization

- Adjust user-tunable variables at the top of `run_unified_pipeline.sh` (e.g., `UPSTREAM_SUBSET`, `DOWNSTREAM_SKIP_DOWNLOAD`, `DOWNSTREAM_MMSEQS_ID`).
- To run only the downstream portion on pre-downloaded proteins, set `DOWNSTREAM_SKIP_DOWNLOAD=1` and place `.faa` files under `unified_pipeline_run/downloads/proteins/` before running.
- For the comparative experiment matrix, set `DOWNSTREAM_MODE="comparative_experiment_matrix"`.

## Contents

```
run_unified_pipeline.sh
nif_hdk_scan_release_clean/
  README.md
  config.yaml
  pipeline.py
  run_nifhdk_repro.sh
  setup_refs.sh
  util.py
  scripts/
    enrich_and_join.py
    extract_ncbi_assembly_quality.py
    summarize_robust_fixcase.py
nif_downstream_code/
  pangenome_pipeline_consolidated2.py
  02_download_proteins.py
  06_merge_annotations_selected.py
  07_expand_gene_families.py
  08_build_narrative_table.py
```

## Notes

- No sample results or archives are included; the upstream scan will populate `nif_hdk_scan_release_clean/results/` on first run.
- Keep this folder as the root of your new repository so relative paths in the scripts stay valid.
