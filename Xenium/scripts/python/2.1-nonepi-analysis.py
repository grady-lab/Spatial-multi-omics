#!/usr/bin/env python3
"""Split clustered non-epithelial adata into immune and stromal subsets.

Usage:
    python 2.1-nonepi-analysis.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc

IMMUNE_CLUSTERS = ["0", "3"]
STROMAL_CLUSTERS = ["1", "2"]


def main():
    ensure_dirs()

    nonepi = sc.read_h5ad(
        P.processed.adata.nonepithelial / "nonepi_16_50_harmonized_clustered.h5ad"
    )
    print(f"Loaded {nonepi.n_obs} non-epithelial cells")

    # ── Split ────────────────────────────────────────────────────────────────
    immune = nonepi[nonepi.obs["nonepi_0.2"].isin(IMMUNE_CLUSTERS)]
    stromal = nonepi[nonepi.obs["nonepi_0.2"].isin(STROMAL_CLUSTERS)]
    print(f"Immune: {immune.n_obs} cells | Stromal: {stromal.n_obs} cells")

    # ── Write ────────────────────────────────────────────────────────────────
    immune_out = P.processed.adata.immune / "immune_cells.h5ad"
    stromal_out = P.processed.adata.stromal / "stromal_cells.h5ad"
    immune_out.parent.mkdir(parents=True, exist_ok=True)
    stromal_out.parent.mkdir(parents=True, exist_ok=True)

    immune.write_h5ad(immune_out)
    stromal.write_h5ad(stromal_out)
    print(f"Wrote {immune_out}")
    print(f"Wrote {stromal_out}")


if __name__ == "__main__":
    main()
