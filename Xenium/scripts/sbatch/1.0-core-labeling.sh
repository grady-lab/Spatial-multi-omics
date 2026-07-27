#!/bin/bash

#SBATCH --mail-user=lallen2@fredhutch.org
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=./slurmout/%j.out
#SBATCH --cpus-per-task=4

# >>> Conda setup >>>
source /home/lallen2/miniforge3/etc/profile.d/conda.sh
conda activate spatial-clean
# <<< Conda setup <<<

# Reads per-slide h5ad files from data/interim/adata/
# Merges metadata from data/metadata/core_annotations.csv
# Writes combined output to data/processed/adata/all_cells/
/usr/bin/time -v ./scripts/python/1.0-core-labeling.py -k A1 A2
