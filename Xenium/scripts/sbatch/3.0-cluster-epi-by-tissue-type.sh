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

# Reads data/processed/adata/epithelial/epithelial_adata.h5ad
# Writes per-tissue-type clustered outputs to the same directory.
/usr/bin/time -v ./scripts/python/3.0-cluster-epi-by-tissue-type.py -p 50 -n 30
