#!/usr/bin/env python3
"""Merge immune/stromal annotations, assign cell-level tissue type labels
via polygon annotations, create normal-split columns, and write annotated
all-cells and epithelial adatas.

Usage:
    python 2.4-final-annotation.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc
import pandas as pd
import numpy as np

INTERMEDIATE_CLUSTERS = ["10"]


def assign_label(row, lookup_dict, fallback_col):
    if row["annotation_final_coarse"] != "Epithelial":
        return "NA"
    core = row["core_id"]
    if core not in lookup_dict:
        return row[fallback_col]
    entry = lookup_dict[core]
    return entry["poly_label"] if row["cell_id"] in entry["cells"] else entry["other_label"]


def main():
    ensure_dirs()

    adata = sc.read_h5ad(
        P.processed.adata.all_cells
        / "xenium_allcells_combo_all_cells_16_50_harmonized_clustered.h5ad"
    )
    print(f"Loaded {adata.n_obs} cells")

    immune_anno = pd.read_csv(P.annotations / "immune-annotations_04.csv")
    stromal_anno_fine = pd.read_csv(P.annotations / "stromal-annotations_fine_08.csv")

    immune_map = immune_anno.drop_duplicates("cell_id").set_index("cell_id")[
        "immune_annotation_fine"
    ]
    stromal_map_fine = stromal_anno_fine.drop_duplicates("cell_id").set_index(
        "cell_id"
    )["stromal_annotation_fine"]

    adata.obs["immune_anno_fine"] = adata.obs["cell_id"].map(immune_map)
    adata.obs["stromal_anno_fine"] = adata.obs["cell_id"].map(stromal_map_fine)

    adata.obs["annotation_final_fine"] = np.where(
        adata.obs["stromal_anno_fine"].notna(),
        adata.obs["stromal_anno_fine"],
        np.where(
            adata.obs["immune_anno_fine"].notna(),
            adata.obs["immune_anno_fine"],
            "Epithelial",
        ),
    )

    adata.obs["annotation_final_coarse"] = np.where(
        adata.obs["immune_anno_fine"].notna(),
        "Immune",
        np.where(
            adata.obs["stromal_anno_fine"].notna()
            & (adata.obs["stromal_anno_fine"] != "Epithelial"),
            "Stromal",
            "Epithelial",
        ),
    )

    print(f"\nannotation_final_fine:\n{adata.obs['annotation_final_fine'].value_counts().to_string()}")
    print(f"\nannotation_final_coarse:\n{adata.obs['annotation_final_coarse'].value_counts().to_string()}")

    adata = adata[
        ~adata.obs["combo_all_cells_0.6"].isin(INTERMEDIATE_CLUSTERS)
    ].copy()
    print(f"\nAfter removing intermediate clusters: {adata.n_obs} cells")

    a1_annos = pd.read_csv(P.metadata / "A1_annos_3-4-26.csv", header=10)
    a2_annos = pd.read_csv(P.metadata / "A2_annos_3-4-26.csv", header=10)
    mixed_core_groups = pd.read_csv(
        P.metadata / "a1_a2_mixed_core_anno_groups.csv"
    )

    a1_annos["core"] = a1_annos["Annotation Name"].str.split("_").str[0]
    a2_annos["core"] = a2_annos["Annotation Name"].str.split("_").str[0]
    all_annos = pd.concat([a1_annos, a2_annos], ignore_index=True)

    cols_to_add = [
        "core", "polygons_are_tissue_type", "other_is_tissue_type",
        "polygons_are_dysplasia", "other_is_dysplasia",
    ]
    all_annos = all_annos.merge(mixed_core_groups[cols_to_add], on="core", how="left")

    unmatched = all_annos[all_annos["polygons_are_tissue_type"].isna()]["core"].unique()
    if len(unmatched) > 0:
        print(f"Warning: cores with no match in mixed_core_groups: {unmatched}")

    tissue_type_dict = {}
    dysplasia_dict = {}
    for core_id, group in all_annos.groupby("core"):
        tissue_type_dict[core_id] = {
            "cells": set(group["Cell ID"].values),
            "poly_label": group["polygons_are_tissue_type"].iloc[0],
            "other_label": group["other_is_tissue_type"].iloc[0],
        }
        dysplasia_dict[core_id] = {
            "cells": set(group["Cell ID"].values),
            "poly_label": group["polygons_are_dysplasia"].iloc[0],
            "other_label": group["other_is_dysplasia"].iloc[0],
        }

    adata.obs["tissue_type_cell_level"] = adata.obs.apply(
        assign_label, axis=1, lookup_dict=tissue_type_dict, fallback_col="tissue_type"
    )
    adata.obs["tissue_type_dysplasia_cell_level"] = adata.obs.apply(
        assign_label, axis=1, lookup_dict=dysplasia_dict, fallback_col="tissue_type_dysplasia"
    )

    epi_mask = adata.obs["annotation_final_coarse"] == "Epithelial"
    nonepi_mask = ~epi_mask
    for col in ["tissue_type_cell_level", "tissue_type_dysplasia_cell_level"]:
        epi_na = (epi_mask & (adata.obs[col] == "NA")).sum()
        nonepi_ok = (nonepi_mask & (adata.obs[col] != "NA")).sum()
        print(f"{col}: epithelial with NA = {epi_na} | non-epithelial with non-NA = {nonepi_ok}")

    cores_in_annos = set(all_annos["core"].unique())
    not_in_annos = set(adata.obs["core_id"].unique()) - cores_in_annos
    fallback_check = (
        adata.obs[adata.obs["core_id"].isin(not_in_annos) & epi_mask]
        .groupby("core_id")["tissue_type_cell_level"]
        .nunique()
    )
    dirty_fallback = fallback_check[fallback_check > 1]
    if len(dirty_fallback) > 0:
        print(f"WARNING: fallback cores with multiple labels: {dirty_fallback.index.tolist()}")
    else:
        print("All fallback cores are clean.")

    mixed_cores = adata.obs[adata.obs["mixed_core"]]["core_id"].unique()
    print("\n=== Verifying mixed cores have multiple epithelial labels ===\n")
    for core in sorted(mixed_cores):
        core_epi = adata.obs[epi_mask & (adata.obs["core_id"] == core)]
        tissue_labels = core_epi["tissue_type_cell_level"].unique().tolist()
        dysplasia_labels = core_epi["tissue_type_dysplasia_cell_level"].unique().tolist()
        print(f"{core} (n_epi={len(core_epi)})")
        print(f"  tissue_type_cell_level:           {tissue_labels}")
        print(f"  tissue_type_dysplasia_cell_level: {dysplasia_labels}")

    for col in ["tissue_type_cell_level", "tissue_type_dysplasia_cell_level"]:
        adata.obs[col] = adata.obs[col].astype(str).replace({"N": "Adj_N", "TA": "AD"})

    for col in ["tissue_type_cell_level", "tissue_type_dysplasia_cell_level"]:
        split_col = col + "_normal_split"
        adata.obs[split_col] = adata.obs[col].astype(str)
        adata.obs.loc[adata.obs["distant_normal"], split_col] = "Dist_N"

    print(f"\ntissue_type_cell_level_normal_split:\n{adata.obs['tissue_type_cell_level_normal_split'].value_counts().to_string()}")
    print(f"\ntissue_type_dysplasia_cell_level_normal_split:\n{adata.obs['tissue_type_dysplasia_cell_level_normal_split'].value_counts().to_string()}")

    for col in ["immune_anno_fine", "stromal_anno_fine"]:
        adata.obs[col] = adata.obs[col].fillna("NA").astype(str)

    epithelial = adata[adata.obs["annotation_final_coarse"] == "Epithelial"].copy()

    all_out = P.processed.adata.all_cells / "all_cells_annotated.h5ad"
    epi_out = P.processed.adata.epithelial / "epithelial_adata.h5ad"
    all_out.parent.mkdir(parents=True, exist_ok=True)
    epi_out.parent.mkdir(parents=True, exist_ok=True)

    adata.write_h5ad(all_out)
    epithelial.write_h5ad(epi_out)
    print(f"\nWrote {all_out}")
    print(f"Wrote {epi_out}")


if __name__ == "__main__":
    main()
