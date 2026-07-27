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

/usr/bin/time -v ./scripts/python/2.1-nonepi-analysis.py
