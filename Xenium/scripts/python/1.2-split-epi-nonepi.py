#!/usr/bin/env python3
"""Remove intermediate cluster and split all-cells adata into epithelial
and non-epithelial subsets.

Usage:
    python 1.2-split-epi-nonepi.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc

INTERMEDIATE_CLUSTERS = ["10"]
EPI_CLUSTERS = ["0", "1", "4", "6", "9", "11"]


def main():
    ensure_dirs()

    adata = sc.read_h5ad(
        P.processed.adata.all_cells
        / "xenium_allcells_combo_all_cells_16_50_harmonized_clustered.h5ad"
    )
    print(f"Loaded {adata.n_obs} cells")

    # ── Remove intermediate cluster ──────────────────────────────────────────
    adata = adata[
        ~adata.obs["combo_all_cells_0.6"].isin(INTERMEDIATE_CLUSTERS)
    ].copy()
    print(f"After removing intermediate clusters: {adata.n_obs} cells")

    # ── Split ────────────────────────────────────────────────────────────────
    nonepi = adata[
        ~adata.obs["combo_all_cells_0.6"].isin(EPI_CLUSTERS)
    ].copy()
    print(f"Non-epithelial: {nonepi.n_obs} cells")

    # ── Write ────────────────────────────────────────────────────────────────
    out = P.processed.adata.nonepithelial / "nonepithelial_adata.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    nonepi.write(out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
