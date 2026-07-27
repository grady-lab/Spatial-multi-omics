"""Shared clustering utilities used across the pipeline scripts."""
from __future__ import annotations

from typing import Optional, Sequence

import scanpy as sc


DEFAULT_RESES: Sequence[float] = (
    0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
)


def cluster_xenium(
    adata,
    n_neighbors: int,
    n_pcs: int,
    key: str,
    output_path,
    scale: str = "No",
    hvgs: str = "No",
    reses: Sequence[float] = DEFAULT_RESES,
    harmony: Optional[Sequence[str]] = None,
    thetas: Optional[Sequence[float]] = None,
):
    """Set X to the normalized layer, optionally Harmony-integrate, then PCA → UMAP → Louvain at multiple resolutions, and write the result.

    Louvain keys are stored in `adata.obs` as `{key}_{resolution}`.
    """
    adata.X = adata.layers["normalized"].copy()
    if hvgs != "No":
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        print("done with hvgs")
    if scale != "No":
        sc.pp.scale(adata)
        print("done with scaling")
    sc.tl.pca(adata, svd_solver="arpack")
    print("done with pca")

    # Restore X so any caller-side use of .X sees normalized values, not the
    # scaled matrix that sc.pp.scale leaves behind when scale != "No".
    adata.X = adata.layers["normalized"].copy()

    if harmony is not None:
        sc.external.pp.harmony_integrate(
            adata,
            key=harmony,
            basis="X_pca",
            adjusted_basis="X_pca_harmony",
            theta=thetas,
        )
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca_harmony")
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca")
    print("done with finding neighbors")

    sc.tl.umap(adata)
    print("done with umap")

    for r in reses:
        sc.tl.louvain(adata, resolution=r, key_added=f"{key}_{r}", random_state=42)
        print(f"done with resolution {r}")

    adata.write_h5ad(output_path)
    return adata
