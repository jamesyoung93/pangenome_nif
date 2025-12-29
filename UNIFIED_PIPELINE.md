# Unified nifH/D/K + Downstream Pangenome Pipeline

This wrapper connects the `nif_hdk_scan_release_clean` scan with the downstream
pangenome analysis so the full workflow can be launched with one command.
Key knobs are defined at the top of `run_unified_pipeline.sh` so you can change
threads, MMseqs2 cutoffs, CV grouping, or the upstream subset without hunting
through multiple scripts.

## Prerequisites
- Internet access for NCBI downloads
- `bash`, `python3`, `mmseqs2`, and `hmmer` on your PATH (the wrapper now exits
  immediately if `hmmsearch` is missing so you don't waste time downloading
  proteomes)
- Python packages required by `nif_hdk_scan_release_clean` and
  `nif_downstream_code/pangenome_pipeline_consolidated2.py`
- `ENTREZ_EMAIL` set (can also be edited in the script header)

For the nifH/D/K scan, the quickest setup is:
```bash
mamba env create -f nif_hdk_scan_release_clean/env/environment.yml
mamba activate nifhdk
module load hmmer/3.4   # or install hmmer into the env
# Verify HMMER is active
which hmmsearch
which hmmbuild
```

## One-line run (default settings)
```bash
./run_unified_pipeline.sh
```
Outputs land in `unified_pipeline_run/results/`.

### Fast rerun from an HPC scratch directory
If you cloned into a scratch path (for example
`/mmfs1/scratch/jacks.local/jyoung67391/pangenome_nif`) and want to restart
cleanly:
1. Stop any running job for this pipeline.
2. Remove the old checkout to clear prior outputs: `rm -rf /mmfs1/scratch/jacks.local/jyoung67391/pangenome_nif`.
3. Recreate and clone fresh:
   ```bash
   mkdir -p /mmfs1/scratch/jacks.local/jyoung67391/pangenome_nif
   cd /mmfs1/scratch/jacks.local/jyoung67391/pangenome_nif
   git clone https://github.com/jamesyoung93/pangenome_nif.git .
   ```
4. Activate your env and reinstall requirements if needed:
   `python -m pip install -r requirements.txt` (and `module load hmmer/3.4` on this
   cluster since HMMER is provided as a module).
5. Export your NCBI contact email and rerun:
   ```bash
   export ENTREZ_EMAIL="you@example.org"
   ./run_unified_pipeline.sh
   ```


## Adjusting parameters
Open `run_unified_pipeline.sh` and edit the block under
`# --- User-tunable settings ---`. Examples:
- Restrict the upstream fetch to complete genomes (default). To broaden it, set
  `UPSTREAM_ASSEMBLY_LEVELS="complete genome,chromosome,scaffold,contig"` or any
  comma-separated list recognized by NCBI.
- Limit the upstream fetch to 200 assemblies for a smoke test:
  `UPSTREAM_SUBSET=200`
- Skip the downstream protein downloads if you already have them in place:
  `DOWNSTREAM_SKIP_DOWNLOAD=1`
- Tighten clustering: set `DOWNSTREAM_MMSEQS_ID=0.9` and `DOWNSTREAM_MMSEQS_COV=0.9`
- Switch to order-blocked CV: `DOWNSTREAM_CV_GROUP="order"`

## What the wrapper does
1. Runs `nif_hdk_scan_release_clean/run_nifhdk_repro.sh` (resumable) with your
   subset/metadata options.
2. Copies the enriched hits table into `unified_pipeline_run/` as the downstream
   starting point.
3. Launches `nif_downstream_code/pangenome_pipeline_consolidated2.py` with the
   knobs you set at the top of the wrapper.

If you prefer batch mode on an HPC, wrap `./run_unified_pipeline.sh` in your
SLURM/LSF submission script after exporting `ENTREZ_EMAIL` and activating your
conda/venv.
