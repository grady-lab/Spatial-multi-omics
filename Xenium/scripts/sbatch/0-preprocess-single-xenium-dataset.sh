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

# Slide keys to preprocess. Each must have an entry under `slides:` in
# config.local.yaml. The script resolves the raw Xenium dir and the metadata
# CSV from config; outputs land in data/interim/{adata,sdata}/.
KEYS=("A1" "A2")

for KEY in "${KEYS[@]}"; do
    /usr/bin/time -v ./scripts/python/0-preprocess-single-xenium-dataset.py -k "$KEY"
done
