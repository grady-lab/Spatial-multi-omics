#!/bin/bash

#SBATCH --mail-user=lallen2@fredhutch.org
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=./slurmout/%j.out
#SBATCH --cpus-per-task=10

# >>> Conda setup >>>
source /home/lallen2/miniforge3/etc/profile.d/conda.sh
conda activate spatial-clean
# <<< Conda setup <<<

# Reads data/processed/adata/all_cells/all_cells_final_anno.h5ad
# Writes T-cell sub-clustered output to data/processed/adata/tcells/.
/usr/bin/time -v ./scripts/python/4.1-subcluster-t-cells.py -p 50 -n 30
