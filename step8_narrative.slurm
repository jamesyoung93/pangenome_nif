#!/bin/bash
#SBATCH -J diazo-narrative
#SBATCH -p bigmem
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o step8_%j.out
#SBATCH -e step8_%j.err

module load python
# choose how many features to include: 50 or 100
export TOP_N=${TOP_N:-100}

python 08_build_narrative_table.py
