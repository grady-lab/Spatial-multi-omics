#!/usr/bin/env python3
"""Run the Cell2location workflow used in the study.

This is a path-parameterized release version of the authoritative analysis
script at ``Data/Python/Deconv/run_cell2location.py``. Analysis stages and
defaults match that script; only file paths and output names are made portable.
The workflow is run independently and is not launched by the R notebooks.
"""

import argparse
import gc
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
import torch
from cell2location import run_colocation
from cell2location.utils.filtering import filter_genes


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train the study Cell2location models and export posterior results."
    )
    parser.add_argument("--spatial", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cell-type-key", default="my_annotation")
    parser.add_argument("--reference-batch-key", default="sample_id")
    parser.add_argument("--spatial-batch-key", default="sample_id")
    parser.add_argument("--reference-epochs", type=int, default=500)
    parser.add_argument("--spatial-epochs", type=int, default=5000)
    parser.add_argument("--train-batch-size", type=int, default=2500)
    parser.add_argument("--posterior-batch-size", type=int, default=5000)
    parser.add_argument("--posterior-samples", type=int, default=1000)
    parser.add_argument("--mean-cells-per-spot", type=float, default=10)
    parser.add_argument("--alpha", type=float, default=20)
    parser.add_argument("--num-workers", type=int, default=12)
    parser.add_argument("--colocation-min-factors", type=int, default=6)
    parser.add_argument("--colocation-max-factors", type=int, default=14)
    parser.add_argument("--colocation-restarts", type=int, default=3)
    parser.add_argument("--cpu", action="store_true")
    return parser.parse_args()


def filter_reference(reference):
    names = reference.var_names
    remove = (
        names.str.startswith(("MT-", "mt-"))
        | names.str.startswith(("RPS", "RPL", "rps", "rpl"))
        | names.str.startswith(("HBA", "HBB", "hba", "hbb"))
    )
    reference = reference[:, ~remove].copy()
    selected = filter_genes(
        reference,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    return reference[:, selected].copy()


def train_or_load_reference_signatures(reference, spatial, args, accelerator):
    signature_path = args.output_dir / "reference_signatures.csv"
    if signature_path.exists():
        print(f"Loading existing reference signatures: {signature_path}")
        signatures = pd.read_csv(signature_path, index_col=0)
        return signatures.loc[spatial.var_names]

    print("Training the reference regression model")
    cell2location.models.RegressionModel.setup_anndata(
        adata=reference,
        labels_key=args.cell_type_key,
        batch_key=args.reference_batch_key,
    )
    reference_model = cell2location.models.RegressionModel(reference)
    reference_model.train(
        max_epochs=args.reference_epochs,
        batch_size=args.train_batch_size,
        train_size=1,
        accelerator=accelerator,
        logger=True,
        enable_checkpointing=False,
    )
    reference = reference_model.export_posterior(
        reference,
        sample_kwargs={
            "num_samples": args.posterior_samples,
            "batch_size": args.train_batch_size,
        },
    )
    if "means_per_cluster_mu_fg" not in reference.varm:
        raise KeyError("Reference signatures are absent from reference.varm")
    factor_names = reference.uns["mod"]["factor_names"]
    signatures = reference.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{name}" for name in factor_names]
    ]
    signatures.columns = factor_names
    signatures.to_csv(signature_path)

    del reference_model
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return signatures


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    accelerator = "cpu" if args.cpu else "gpu"
    if accelerator == "gpu" and not torch.cuda.is_available():
        raise RuntimeError("GPU requested but unavailable; pass --cpu to run on CPU")
    scvi.settings.dl_num_workers = args.num_workers

    print("Torch:", torch.__version__)
    print("CUDA available:", torch.cuda.is_available())
    print("Loading raw-count AnnData inputs")
    spatial = sc.read_h5ad(args.spatial)
    reference = filter_reference(sc.read_h5ad(args.reference))

    for adata in (spatial, reference):
        if sp.issparse(adata.X) and not sp.isspmatrix_csr(adata.X):
            adata.X = adata.X.tocsr()
    intersect = spatial.var_names.intersection(reference.var_names)
    spatial = spatial[:, intersect].copy()
    reference = reference[:, intersect].copy()
    print(f"Intersected gene count: {len(intersect)}")

    signatures = train_or_load_reference_signatures(
        reference, spatial, args, accelerator
    )

    print("Training the spatial mapping model")
    cell2location.models.Cell2location.setup_anndata(
        adata=spatial,
        batch_key=args.spatial_batch_key,
    )
    spatial_model = cell2location.models.Cell2location(
        spatial,
        cell_state_df=signatures,
        N_cells_per_location=args.mean_cells_per_spot,
        detection_alpha=args.alpha,
    )
    spatial_model.train(
        max_epochs=args.spatial_epochs,
        batch_size=args.train_batch_size,
        train_size=1,
        accelerator=accelerator,
    )
    spatial = spatial_model.export_posterior(
        spatial,
        sample_kwargs={
            "num_samples": args.posterior_samples,
            "batch_size": args.posterior_batch_size,
        },
    )

    posterior_path = args.output_dir / "cell2location_posterior.h5ad"
    model_path = args.output_dir / "cell2location_model"
    spatial.write_h5ad(posterior_path)
    spatial_model.save(model_path, overwrite=True)

    q05 = spatial.obsm["q05_cell_abundance_w_sf"].copy()
    q05.columns = [
        column.replace("q05cell_abundance_wsf_", "") for column in q05.columns
    ]
    abundance = pd.DataFrame(
        q05,
        index=spatial.obs_names,
        columns=spatial.uns["mod"]["factor_names"],
    )
    abundance.to_csv(args.output_dir / "cell2location_q05.csv")

    print("Running NMF cell-type colocalization")
    _, spatial = run_colocation(
        spatial,
        model_name="CoLocatedGroupsSklearnNMF",
        train_args={
            "n_fact": np.arange(
                args.colocation_min_factors,
                args.colocation_max_factors + 1,
            ),
            "sample_name_col": args.spatial_batch_key,
            "n_restarts": args.colocation_restarts,
        },
        model_kwargs={
            "alpha": 0.001,
            "init": "random",
            "nmf_kwd_args": {"tol": 0.000001},
        },
        export_args={"path": str(args.output_dir / "NMF_colocation")},
    )
    spatial.write_h5ad(posterior_path)

    print("Computing cell-type-specific expected expression")
    expected = spatial_model.module.model.compute_expected_per_cell_type(
        spatial_model.samples["post_sample_q05"],
        spatial_model.adata_manager,
    )
    for index, cell_type in enumerate(spatial_model.factor_names_):
        spatial.layers[cell_type] = expected["mu"][index]
    spatial.write_h5ad(posterior_path)

    print("Cell2location workflow complete")
    print("Abundances:", args.output_dir / "cell2location_q05.csv")
    print("Posterior with expected-expression layers:", posterior_path)
    print("Saved model:", model_path)
    print("NMF colocalization:", args.output_dir / "NMF_colocation")


if __name__ == "__main__":
    main()
