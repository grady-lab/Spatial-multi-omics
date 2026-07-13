#!/usr/bin/env python3
"""Run stLearn spatial pseudotime independently within tissue islands.

This standalone script is distilled from
Data/Python/h5ad_files/LoopforPseudotime.ipynb. It is intentionally not
launched by the R workflow.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors


EPITHELIAL_CLUSTERS = {"0", "1", "2", "3", "4", "6", "8"}
EPITHELIAL_ANNOTATIONS = {"Normal", "Low-grade"}


def tissue_islands(adata, k=8, cutoff=2.0):
    coords = adata.obs[["array_row", "array_col"]].to_numpy(dtype=float)
    n_obs = coords.shape[0]
    neighbors = NearestNeighbors(n_neighbors=min(k + 1, n_obs)).fit(coords)
    distances, indices = neighbors.kneighbors(coords)
    rows = []
    cols = []
    for i in range(n_obs):
        for j, distance in zip(indices[i, 1:], distances[i, 1:]):
            if distance <= cutoff:
                rows.extend((i, int(j)))
                cols.extend((int(j), i))
    adjacency = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_obs, n_obs))
    _, labels = connected_components(adjacency, directed=False)
    return labels


def adaptive_epsilon(adata, quantile=0.9, multiplier=1.5):
    coords = adata.obs[["imagerow", "imagecol"]].to_numpy(dtype=float)
    nearest = NearestNeighbors(n_neighbors=2).fit(coords)
    distances, _ = nearest.kneighbors(coords)
    return float(np.quantile(distances[:, 1], quantile) * multiplier)


def root_spot(adata, cluster_col, root_cluster="1"):
    labels = adata.obs[cluster_col].astype(str)
    mask = (labels == root_cluster).to_numpy()
    if not mask.any():
        raise ValueError(f"Root cluster {root_cluster!r} is absent from this tissue island")
    coordinates = adata.obsm["X_umap"][mask]
    center = np.median(coordinates, axis=0)
    within_cluster = np.where(mask)[0]
    return int(within_cluster[np.argmin(np.linalg.norm(coordinates - center, axis=1))])


def run_one_island(adata, cluster_col, root_cluster):
    st.pp.filter_genes(adata, min_cells=3)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)
    adata.raw = adata
    st.pp.scale(adata)

    n_pcs = min(50, adata.n_obs - 1, adata.n_vars - 1)
    n_neighbors = min(25, adata.n_obs - 1)
    adata.obs[cluster_col] = adata.obs[cluster_col].astype(str).astype("category")
    st.em.run_pca(adata, n_comps=n_pcs, random_state=0)
    st.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=0)
    sc.tl.umap(adata, random_state=0)

    root = root_spot(adata, cluster_col=cluster_col, root_cluster=root_cluster)
    adata.uns["iroot"] = root
    epsilon = adaptive_epsilon(adata)

    if "sub_cluster_labels" in adata.obs:
        adata.obs.drop(columns=["sub_cluster_labels"], inplace=True)
    for key in ("global_graph", "paga", "split_node", "centroid_dict", "threshold_spots"):
        adata.uns.pop(key, None)

    st.spatial.clustering.localization(
        adata, use_label=cluster_col, eps=epsilon, min_samples=1
    )
    st.spatial.trajectory.pseudotime(
        adata,
        use_label=cluster_col,
        eps=epsilon,
        threshold=0.001,
        threshold_spots=0,
        max_nodes=15,
        run_knn=True,
        n_neighbors=8,
        use_rep="X_pca",
        reverse=False,
    )
    pseudotime = adata.obs["dpt_pseudotime"].astype(float)
    finite = np.isfinite(pseudotime)
    if finite.any():
        low, high = pseudotime[finite].min(), pseudotime[finite].max()
        adata.obs.loc[finite, "dpt_pseudotime"] = (
            (pseudotime[finite] - low) / (high - low) if high > low else 0.0
        )
    return adata.obs.copy(), epsilon, root


def process_file(
    path,
    output_dir,
    annotation_col,
    cluster_col,
    root_cluster,
    min_island_spots,
):
    sample_id = path.stem
    adata = st.convert_scanpy(sc.read_h5ad(path))
    labels = adata.obs[cluster_col].astype(str)
    annotation = adata.obs[annotation_col].astype(str).str.strip()
    keep = labels.isin(EPITHELIAL_CLUSTERS) & annotation.isin(EPITHELIAL_ANNOTATIONS)
    adata = adata[keep].copy()
    adata.obs["tissue_island"] = tissue_islands(adata).astype(str)

    outputs = []
    log = []
    for island, count in adata.obs["tissue_island"].value_counts().items():
        if count < min_island_spots:
            log.append({"sample_id": sample_id, "island": island, "status": "too_few_spots"})
            continue
        subset = adata[adata.obs["tissue_island"] == island].copy()
        try:
            result, epsilon, root = run_one_island(subset, cluster_col, root_cluster)
            result["sample_id"] = sample_id
            result["tissue_island"] = island
            result["localization_epsilon"] = epsilon
            result["root_spot_index"] = root
            outputs.append(result)
            log.append({"sample_id": sample_id, "island": island, "status": "ok"})
        except Exception as error:  # keep other samples running and make exclusions auditable
            log.append({"sample_id": sample_id, "island": island, "status": repr(error)})

    if outputs:
        output = pd.concat(outputs, axis=0)
        output.index.name = "barcode"
        output.to_csv(output_dir / f"{sample_id}.pseudotime.csv")
    return log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--annotation-col", default="Adenoma")
    parser.add_argument("--cluster-col", default="combined_cluster")
    parser.add_argument("--root-cluster", default="1")
    parser.add_argument("--min-island-spots", type=int, default=30)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log = []
    for path in sorted(args.input_dir.glob("*.h5ad")):
        log.extend(
            process_file(
                path,
                args.output_dir,
                args.annotation_col,
                args.cluster_col,
                args.root_cluster,
                args.min_island_spots,
            )
        )
    pd.DataFrame(log).to_csv(args.output_dir / "pseudotime_run_log.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
