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

# Reads epithelial adata + misannotated_fibroblasts.txt from defaults
# (data/processed/adata/epithelial/ and data/annotations/).
/usr/bin/time -v ./scripts/python/3.2-recluster-epi-normal.py -p 50 -n 16
