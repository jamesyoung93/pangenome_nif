# nif_hdk_scan (release)
Quick steps:
1) `mamba env create -f env/environment.yml && mamba activate nifhdk`
2) `module load hmmer/3.4`  (or install hmmer in the env)
3) `export ENTREZ_EMAIL="you@yourorg"`
4) Put HMMs under `refs/hmms/nif{H,D,K}/*.hmm` (already included here if you had them).
5) Smoke test: `bash run_nifhdk_repro.sh --subset 200`
6) Full run:  `bash run_nifhdk_repro.sh`
7) Add metadata (optional): `bash run_nifhdk_repro.sh --with-metadata` (the script now auto-installs the NCBI `datasets` CLI via `pip` if missing; you can also install manually with `mamba install -c conda-forge ncbi-datasets-cli`, `conda install -c conda-forge ncbi-datasets-cli`, or `module load ncbi-datasets`)
