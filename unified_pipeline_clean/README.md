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

# run the unified workflow (defaults: full upstream scan + downstream **full** mode)
./run_unified_pipeline.sh
```

Outputs land in `unified_pipeline_run/results/`. Modeling artifacts are copied into `unified_pipeline_run/results/modeling/`.

### Genome filtering policy

- By default, downstream steps **only** accept assemblies that are both `GCF_*` (RefSeq) **and** `Complete Genome`.
- The upstream `assembly_quality.tsv` is staged to `unified_pipeline_run/results/assemblies/` and filtered into:
  - `assembly_quality.filtered.tsv`
  - `filtered_accessions.txt`
  - `filter_summary.json`
- To relax the defaults, pass `--allow-gca` and/or `--allow-noncomplete` through downstream args (or edit `config.yaml`).

### Modes and modeling outputs

- The default **full** mode runs the complete modeling flow (classification + directionality) via `04_classify.py` and `05_analyze_feature_directionality.py`.
- Key artifacts are copied into `unified_pipeline_run/results/modeling/`:
  - `classification_summary.csv`, `feature_importance_top.csv`, `feature_importance_*.csv`, `fold_metrics_*.csv`, `roc_curves.png`
  - `feature_directionality_full.csv`, `feature_directionality_summary.csv`, `feature_directionality_*.csv`
  - `narrative_top_features.csv` and plots such as `feature_importance_plot.png`
- Use **experiment** mode to skip heavy directionality (Step 5). The pipeline will print a warning and directionality outputs will be absent.
- The **comparative_experiment_matrix** mode remains available for the multi-arm benchmarking workflow.

## Customization

- Adjust user-tunable variables at the top of `run_unified_pipeline.sh` (e.g., `UPSTREAM_SUBSET`, `DOWNSTREAM_SKIP_DOWNLOAD`, `DOWNSTREAM_MMSEQS_ID`).
- To run only the downstream portion on pre-downloaded proteins, set `DOWNSTREAM_SKIP_DOWNLOAD=1` and place `.faa` files under `unified_pipeline_run/downloads/proteins/` before running.
- For the comparative experiment matrix, set `DOWNSTREAM_MODE="comparative_experiment_matrix"`.
- If HMMER is installed outside your default `PATH`, set `HMMSEARCH_CMD` to the full `hmmsearch` path or set `HMMER_MODULE` (e.g., `hmmer/3.4`) so the wrapper can `module load` it automatically.
- Ensure `mmseqs` and `hmmsearch` are available on PATH (for example HPC modules `mmseqs2/15-6f452` and `hmmer/3.4`). Python requirements include scikit-learn, seaborn, matplotlib, and xgboost for modeling steps.

## Contents

```
run_unified_pipeline.sh
nif_downstream_code/            # uses root-level canonical scripts
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
- Post-run validation (under `unified_pipeline_run/results/`):
  - `assemblies/filtered_accessions.txt` exists and lists only `GCF_*` complete genomes
  - `modeling/feature_importance_top.csv`, `modeling/fold_metrics_*.csv`, and `modeling/feature_importance_*.csv` are present
