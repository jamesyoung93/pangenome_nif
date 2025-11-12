#!/bin/bash
#SBATCH -J diazo-gf
#SBATCH -p bigmem
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o step3_%j.out
#SBATCH -e step3_%j.err

module load mmseqs2
module load python

export TMPDIR=/mmfs1/scratch/jacks.local/$USER/mmseqs_tmp_$SLURM_JOB_ID
export MMSEQS_TMP=$TMPDIR
mkdir -p "$TMPDIR"

# clean stale concat, then run
rm -f all_proteins.faa gene_families_DB* gene_families_cluster* gene_families_clusters.tsv

python 03_create_gene_families.py 16
