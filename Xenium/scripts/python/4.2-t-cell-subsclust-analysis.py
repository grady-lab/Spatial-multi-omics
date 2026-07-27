#!/usr/bin/env python3
"""Annotate T-cell subclusters, create annotation_final_fine_cd8 on the
all-cells adata using cluster membership, export annotation CSV, and
write CD8 T-cell subset.

The CD8 cluster is identified automatically as the cluster with the
highest median CD8A expression at the T-Cells_0.6 resolution.

Usage:
    python 4.2-t-cell-subsclust-analysis.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import numpy as np
import scanpy as sc
import pandas as pd

CLUSTER_COL = "T-Cells_0.6"
CD8A_GENE = "CD8A"


def find_cd8_cluster(adata, cluster_col=CLUSTER_COL, gene=CD8A_GENE):
    """Return the cluster label with the highest median CD8A expression."""
    expr = adata[:, gene].layers["normalized"]
    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()
    else:
        expr = np.asarray(expr).flatten()

    clusters = adata.obs[cluster_col].astype(str)
    medians = {}
    for cl in clusters.unique():
        medians[cl] = float(np.median(expr[clusters == cl]))

    best = max(medians, key=medians.get)
    print(f"\nMedian {gene} expression per cluster:")
    for cl in sorted(medians, key=lambda x: int(x)):
        marker = " <-- CD8" if cl == best else ""
        print(f"  cluster {cl}: {medians[cl]:.4f}{marker}")
    return best


def export_xenium_explorer_csv(adata, column, output_dir):
    import os
    os.makedirs(output_dir, exist_ok=True)
    for batch in adata.obs["batch"].unique():
        df = adata.obs.loc[adata.obs["batch"] == batch, ["cell_id", column]].copy()
        df.columns = ["cell_id", "group"]
        out_path = os.path.join(output_dir, f"{column}_{batch}.csv")
        df.to_csv(out_path, index=False)
        print(f"Saved {out_path} ({len(df)} cells)")


def main():
    ensure_dirs()

    t_cell_adata = sc.read_h5ad(
        P.processed.adata.tcells / "tcells_30_50_harmony_05_pt_only.h5ad"
    )
    print(f"Loaded {t_cell_adata.n_obs} T-cells")

    cd8_cluster = find_cd8_cluster(t_cell_adata)

    anno_dict = {}
    for cl in t_cell_adata.obs[CLUSTER_COL].astype(str).unique():
        anno_dict[cl] = "CD8_T_Cells" if cl == cd8_cluster else "T_Cells"

    t_cell_adata.obs["T-Cells_annotation"] = (
        t_cell_adata.obs[CLUSTER_COL].astype(str).map(anno_dict)
    )
    print(f"\nT-Cells_annotation:\n{t_cell_adata.obs['T-Cells_annotation'].value_counts().to_string()}")

    anno_out = P.annotations / "t-cell_annotations.csv"
    t_cell_adata.obs[["cell_id", "T-Cells_annotation"]].to_csv(anno_out, index=False)
    print(f"Wrote {anno_out}")

    export_xenium_explorer_csv(
        t_cell_adata, CLUSTER_COL, str(P.results.xenium_explorer)
    )

    cd8 = t_cell_adata[t_cell_adata.obs[CLUSTER_COL] == cd8_cluster].copy()
    cd8_out = P.processed.adata.tcells / "cd8_t_cells_unclustered.h5ad"
    cd8_out.parent.mkdir(parents=True, exist_ok=True)
    cd8.write(cd8_out)
    print(f"Wrote {cd8.n_obs} CD8 T-cells to {cd8_out}")

    all_cells_path = P.processed.adata.all_cells / P.fn.all_cells_final
    adata = sc.read_h5ad(all_cells_path)
    print(f"\nLoaded {adata.n_obs} all-cells from {all_cells_path}")

    cd8_cell_ids = set(cd8.obs["cell_id"])
    adata.obs["annotation_final_fine_cd8"] = adata.obs["annotation_final_fine"].copy()
    if "CD8-T-Cell" not in adata.obs["annotation_final_fine_cd8"].cat.categories:
        adata.obs["annotation_final_fine_cd8"] = (
            adata.obs["annotation_final_fine_cd8"].cat.add_categories("CD8-T-Cell")
        )
    cd8_mask = adata.obs["cell_id"].isin(cd8_cell_ids)
    adata.obs.loc[cd8_mask, "annotation_final_fine_cd8"] = "CD8-T-Cell"

    print(f"Labeled {cd8_mask.sum()} cells as CD8-T-Cell")
    print(f"\nannotation_final_fine_cd8:\n{adata.obs['annotation_final_fine_cd8'].value_counts().to_string()}")

    adata.write_h5ad(all_cells_path)
    print(f"Updated {all_cells_path}")


if __name__ == "__main__":
    main()
