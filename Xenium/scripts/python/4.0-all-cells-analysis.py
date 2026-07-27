#!/usr/bin/env python3
"""Transfer epithelial labels, fix fibroblast misannotations, and export
Xenium Explorer CSVs.

Note: annotation_final_fine_cd8 is created in step 4.2 after T-cell
subclustering identifies CD8 T-cells by cluster membership.

Usage:
    python 4.0-all-cells-analysis.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc
import numpy as np


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

    adata = sc.read_h5ad(P.processed.adata.all_cells / "all_cells_annotated.h5ad")
    print(f"Loaded {adata.n_obs} cells")

    with open(P.annotations / "misannotated_fibroblasts.txt") as f:
        fibro_misanno = f.read().splitlines()
    print(f"Loaded {len(fibro_misanno)} misannotated fibroblast cell_ids")

    epi_adata = sc.read_h5ad(
        P.processed.adata.epithelial / "epithelial_20_50_harmony_batch_05_pt_05_ssg.h5ad"
    )

    cols_to_transfer = [
        "epithelial_0.8",
        "GDF15_positive", "GDF15_more_than_1_transcript",
        "senepy_intestine_epi_0", "stemness_score", "senepy_high",
    ]
    epi_obs_indexed = epi_adata.obs.set_index("cell_id")[cols_to_transfer]

    for col in cols_to_transfer:
        adata.obs[col] = adata.obs["cell_id"].map(epi_obs_indexed[col])

    for col in cols_to_transfer:
        n_filled = adata.obs[col].notna().sum()
        n_nan = adata.obs[col].isna().sum()
        print(f"{col}: {n_filled} filled, {n_nan} NaN")

    for col in cols_to_transfer:
        adata.obs[col] = adata.obs[col].astype(str).replace("nan", "NA")

    out = P.processed.adata.all_cells / P.fn.all_cells_final
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"\nWrote {out}")

    adata.obs["annotation_final_tissue_type_dysplasia"] = np.where(
        adata.obs["annotation_final_fine"] == "Epithelial",
        adata.obs["tissue_type_dysplasia_cell_level"].astype(str),
        adata.obs["annotation_final_fine"].astype(str),
    )
    print(f"\nannotation_final_tissue_type_dysplasia:\n{adata.obs['annotation_final_tissue_type_dysplasia'].value_counts().to_string()}")

    xenium_dir = str(P.results.xenium_explorer)
    export_xenium_explorer_csv(adata, "annotation_final_fine", xenium_dir)
    export_xenium_explorer_csv(adata, "annotation_final_tissue_type_dysplasia", xenium_dir)


if __name__ == "__main__":
    main()
