#!/usr/bin/env python3
"""Annotate immune and stromal clusters and export annotation CSVs.

Usage:
    python 2.3-immune-and-stromal-analysis.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc

IMMUNE_ANNO_DICT = {
    "T-Cells": ["0"],
    "B-Cells": ["3"],
    "Mast Cells": ["4"],
    "Plasma Cells": ["2"],
    "Myeloid Cells": ["1"],
    "Unknown": ["5"],
}

STROMAL_ANNO_DICT = {
    "Fibroblasts": ["2", "5", "4", "6", "0", "8", "9"],
    "Endothelial": ["3"],
    "Pericytes": ["7"],
    "Smooth Muscle": ["1"],
    "Schwann Cells": ["10"],
    "Unknown": ["11"],
}


def main():
    ensure_dirs()

    immune = sc.read_h5ad(
        P.processed.adata.immune / "immune_clustered_harmonized.h5ad"
    )
    stromal = sc.read_h5ad(
        P.processed.adata.stromal / "stromal_clustered_harmonized.h5ad"
    )
    print(f"Immune: {immune.n_obs} cells | Stromal: {stromal.n_obs} cells")

    cluster_to_anno = {
        c: a for a, clusters in IMMUNE_ANNO_DICT.items() for c in clusters
    }
    immune.obs["immune_annotation_fine"] = immune.obs["immune_0.4"].map(
        cluster_to_anno
    )
    print(f"\nimmune_annotation_fine:\n{immune.obs['immune_annotation_fine'].value_counts().to_string()}")

    immune_out = P.annotations / "immune-annotations_04.csv"
    immune.obs[["cell_id", "immune_annotation_fine"]].to_csv(immune_out)
    print(f"Wrote {immune_out}")

    cluster_to_anno = {
        c: a for a, clusters in STROMAL_ANNO_DICT.items() for c in clusters
    }
    stromal.obs["stromal_annotation_fine"] = stromal.obs["stromal_0.8"].map(
        cluster_to_anno
    )
    print(f"\nstromal_annotation_fine:\n{stromal.obs['stromal_annotation_fine'].value_counts().to_string()}")

    stromal_out = P.annotations / "stromal-annotations_fine_08.csv"
    stromal.obs[["cell_id", "stromal_annotation_fine"]].to_csv(stromal_out)
    print(f"Wrote {stromal_out}")


if __name__ == "__main__":
    main()
