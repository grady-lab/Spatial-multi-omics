#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs
from cluster_utils import cluster_xenium

import scanpy as sc
import anndata as ad
import pandas as pd
import spatialdata as sd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def get_args():
        parser = argparse.ArgumentParser("files!")
        parser.add_argument("-a", "--adata", help="xenium filtered normalized h5ad file",
                            default=str(P.processed.adata.all_cells / P.fn.all_cells_final))
        parser.add_argument("-o", "--adata_output_dir", help="path for adata output",
                            default=str(P.processed.adata.tcells))
        parser.add_argument("-s", "--scale", help="whether to scale", default="No")
        parser.add_argument("-v", "--hvgs", help="whether to use hvgs", default="No")
        parser.add_argument("-p", "--pcs", help="number of pcs", required=True)
        parser.add_argument("-n", "--nns", help="number of ns", required=True)
        return parser.parse_args()

args = get_args()
adata_fp = str(args.adata)
scale = str(args.scale)
hvgs = str(args.hvgs)
pcs = int(str(args.pcs))
nns = int(str(args.nns))

adata_output_dir = str(args.adata_output_dir)

ensure_dirs()



adata_all = sc.read_h5ad(adata_fp) # load your single adata here



cell_ids = adata_all.obs["cell_id"].astype(str)
mask = ~cell_ids.duplicated(keep="first")
adata_all = adata_all[mask].copy()

# Verify
print(f"Cells before: {len(mask)}")
print(f"Cells after:  {adata_all.n_obs}")
print(f"Removed:      {(~mask).sum()}")



adata_t_cells = adata_all[adata_all.obs["annotation_final_fine"] == "T-Cells"].copy()



# First run (harmony)
# adata_t_cells_harmony = adata_t_cells.copy()
# cluster_xenium(
#     adata_t_cells_harmony, nns, pcs,
#     key="T-Cells",
#     output_path=f"{adata_output_dir}/tcells_{nns}_{pcs}_harmony_05_05.h5ad",
#     scale=scale, hvgs=hvgs,
#     harmony=["batch", "patient_id"], thetas=[0.5, 0.5]
# )

# # Second run (no harmony)
# adata_t_cells_no_harmony = adata_t_cells.copy()
# cluster_xenium(
#     adata_t_cells_no_harmony, nns, pcs,
#     key="T-Cells",
#     output_path=f"{adata_output_dir}/tcells_{nns}_{pcs}.h5ad",
#     scale=scale, hvgs=hvgs,
# )

# third run (harmony)
adata_t_cells_harmony = adata_t_cells.copy()
cluster_xenium(
    adata_t_cells_harmony, nns, pcs,
    key="T-Cells",
    output_path=f"{adata_output_dir}/tcells_{nns}_{pcs}_harmony_05_pt_only.h5ad",
    scale=scale, hvgs=hvgs,
    harmony=["patient_id"], thetas=[0.5]
)
