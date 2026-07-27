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
import senepy as sp


def get_args():
        parser = argparse.ArgumentParser("files!")
        parser.add_argument("-e", "--epithelial", help="xenium filtered normalized h5ad file",
                            default=str(P.processed.adata.epithelial / "epithelial_adata.h5ad"))
        parser.add_argument("-o", "--adata_output_dir", help="path for adata output",
                            default=str(P.processed.adata.epithelial))
        parser.add_argument("-s", "--scale", help="whether to scale", default="No")
        parser.add_argument("-v", "--hvgs", help="whether to use hvgs", default="No")
        parser.add_argument("-p", "--pcs", help="number of pcs", required=True)
        parser.add_argument("-n", "--nns", help="number of ns", required=True)
        return parser.parse_args()

args = get_args()
epi_fp = str(args.epithelial)
scale = str(args.scale)
hvgs = str(args.hvgs)
pcs = int(str(args.pcs))
nns = int(str(args.nns))

adata_output_dir = str(args.adata_output_dir)

ensure_dirs()


def run_senepy(adata, celltypes):
    hubs = sp.load_hubs(species='Human')
    translator = sp.translator(hub=hubs.hubs, data=adata)
    if celltypes == "epithelial":
        adata.obs["senepy_intestine_epi_0"] = sp.score_hub(adata, hubs.hubs[('intestine', 'epithelial cell', 0)], translator=translator)
    else:
        hubs.merge_hubs(hubs.metadata, new_name="Universal", calculate_thresh=True, p_thres=0.01)
        adata.obs["senepy_universal"] = sp.score_hub(adata, hubs.hubs['Universal'], translator=translator)
    return adata



def score_stemness_lau(adata):
    # Stemness scoring
    stem_score_genes = ["CLDN2", "CD44", "AXIN2", "RNF43", "TGFBI",
                        "EPHB2", "TEAD2", "CDX2", "LGR5", "OLFM4", "ASCL2"]

    # Check which genes are available
    stem_available_genes = [g for g in stem_score_genes if g in adata.var_names]
    stem_missing_genes = [g for g in stem_score_genes if g not in adata.var_names]

    print(f"\nStemness genes available: {stem_available_genes}")
    print(f"Stemness genes missing: {stem_missing_genes}")

    # Set X to normalized layer for scoring
    if "normalized" in adata.layers:
        adata.X = adata.layers["normalized"].copy()
    else:
        print("Warning: 'normalized' layer not found, using current X")

    # Calculate stemness score
    if len(stem_available_genes) > 0:
        sc.tl.score_genes(
            adata,
            gene_list=stem_available_genes,
            score_name='stemness_score',
            use_raw=False
        )
        print(f"\nStemness score calculated. Mean: {adata.obs['stemness_score'].mean():.3f}")
    else:
        print("Warning: No stemness genes available for scoring")
    return adata

def assign_senepy_high(adata, column='senepy_intestine_epi_0'):
    data = adata.obs[column]
    mean = np.mean(data)
    sd = np.std(data)
    upper_2sd = mean + 2 * sd
    adata.obs['senepy_high'] = adata.obs[column] > upper_2sd
    return adata



adata_epi = sc.read_h5ad(epi_fp)

# Remove duplicated cell_ids, keeping the first occurrence
cell_ids = adata_epi.obs["cell_id"].astype(str)
mask = ~cell_ids.duplicated(keep="first")
adata_epi = adata_epi[mask].copy()

# Verify
print(f"Cells before: {len(mask)}")
print(f"Cells after:  {adata_epi.n_obs}")
print(f"Removed:      {(~mask).sum()}")



adata_epi.obs["GDF15_positive"] = (adata_epi[:, "GDF15"].layers["counts"].toarray().flatten() > 0)
adata_epi.obs["GDF15_more_than_1_transcript"] = (adata_epi[:, "GDF15"].layers["counts"].toarray().flatten() > 1)

adata_epi = run_senepy(adata_epi, "epithelial")
adata_epi = score_stemness_lau(adata_epi)
adata_epi = assign_senepy_high(adata_epi)


adata_epi.obs["tissue_type_cell_level_normal_split"] = adata_epi.obs["tissue_type_cell_level"].astype(str)
adata_epi.obs.loc[adata_epi.obs["distant_normal"] == True, "tissue_type_cell_level_normal_split"] = "Dist_N"
print(adata_epi.obs["tissue_type_cell_level_normal_split"].value_counts())

adata_epi.obs["tissue_type_dysplasia_cell_level_normal_split"] = adata_epi.obs["tissue_type_dysplasia_cell_level"].astype(str)
adata_epi.obs.loc[adata_epi.obs["distant_normal"] == True, "tissue_type_dysplasia_cell_level_normal_split"] = "Dist_N"
print(adata_epi.obs["tissue_type_dysplasia_cell_level_normal_split"].value_counts())




for tissue_type in [t for t in adata_epi.obs["tissue_type_cell_level"].unique() if t != "NA"]:

    print(f"Processing tissue type: {tissue_type}")

    adata_subset = adata_epi[adata_epi.obs["tissue_type_cell_level"] == tissue_type].copy()

    tissue_type_clean = tissue_type.replace(" ", "_").replace("/", "_")
    output_fp = f"{adata_output_dir}/epithelial_{tissue_type_clean}_{nns}_{pcs}_tt_harmony_batch_05_pt_05_ssg.h5ad"

    cluster_xenium(
        adata_subset,
        nns,
        pcs,
        key=f"epithelial_{tissue_type_clean}",
        output_path=output_fp,
        scale=scale,
        hvgs=hvgs,
        harmony=["batch", "patient_id"],
        thetas=[0.5, 0.5]
    )

    print(f"Saved {tissue_type} to {output_fp}")


adata_dist_N = adata_epi[adata_epi.obs["tissue_type_cell_level_normal_split"] == "Dist_N"].copy()

output_fp = f"{adata_output_dir}/epithelial_distant_normal_{nns}_{pcs}.h5ad"
cluster_xenium(
    adata_dist_N, nns, pcs,
    key="distant_normal",
    output_path=output_fp,
    scale=scale,
    hvgs=hvgs,
    harmony=["batch", "patient_id"],
    thetas=[0.5, 0.5]
)


adata_N_not_sep = adata_epi[adata_epi.obs["tissue_type_cell_level_normal_split"] == "Adj_N"].copy()

output_fp = f"{adata_output_dir}/epithelial_adjacent_normal_{nns}_{pcs}.h5ad"


cluster_xenium(
    adata_N_not_sep, nns, pcs,
    key="adjacent_normal",
    output_path=output_fp,
    scale=scale,
    hvgs=hvgs,
    harmony=["batch", "patient_id"],
    thetas=[0.5, 0.5]
)
