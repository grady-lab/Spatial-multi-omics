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

# Reads data/processed/adata/nonepithelial/nonepithelial_adata.h5ad
# Writes clustered outputs to the same directory.
/usr/bin/time -v ./scripts/python/2.0-cluster-nonepithelial.py -p 50 -n 16
# /usr/bin/time -v ./scripts/python/2.0-cluster-nonepithelial.py -p 50 -n 30
