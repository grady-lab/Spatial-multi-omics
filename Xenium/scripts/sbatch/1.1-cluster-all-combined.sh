#!/bin/bash

#SBATCH --mail-user=lallen2@fredhutch.org
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=./slurmout/%j.out
#SBATCH --cpus-per-task=30

# >>> Conda setup >>>
source /home/lallen2/miniforge3/etc/profile.d/conda.sh
conda activate spatial-clean
# <<< Conda setup <<<

# Reads data/processed/adata/all_cells/combined_filtered_normalized_xenium.h5ad
# Writes clustered outputs to the same directory.
/usr/bin/time -v ./scripts/python/1.1-cluster-all-combined.py -k combo_all_cells -p 50 -n 16
# /usr/bin/time -v ./scripts/python/1.1-cluster-all-combined.py -k combo_all_cells -p 50 -n 20
