#!/usr/bin/env python3
"""Concatenate per-slide AnnData objects, apply core-id renames, and merge
metadata from the core-annotations spreadsheet.

Usage:
    python 1.0-core-labeling.py -k A1 A2
    python 1.0-core-labeling.py -k A1 A2 -o /path/to/output.h5ad
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc
import anndata as ad
import pandas as pd


CORE_ID_RENAMES = {
    "Y2-CA-2": "Y2-TA-3",
    "Q1-CA-2": "Q1-TA-3",
    "O1-CA-1": "O1-TA-2",
    "Q1-CA-1": "Q1-TA-4",
    "L1-CA-2": "L1-TA-4",
}


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-k", "--keys", nargs="+", default=["A1", "A2"],
                        help="slide keys to concatenate (default: A1 A2)")
    parser.add_argument("-i", "--input_dir", type=Path, default=None,
                        help=f"directory containing per-slide h5ad files (default: {P.interim.adata})")
    parser.add_argument("-c", "--core_anno", type=Path, default=None,
                        help=f"core annotations CSV (default: {P.metadata}/core_annotations_final.csv)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help=f"output h5ad path (default: {P.processed.adata.all_cells}/<combined_adata>)")
    args = parser.parse_args()
    if args.input_dir is None:
        args.input_dir = P.interim.adata
    if args.core_anno is None:
        args.core_anno = P.metadata / "core_annotations_final.csv"
    if args.output is None:
        args.output = P.processed.adata.all_cells / P.fn.combined_adata
    return args


def main():
    args = get_args()
    ensure_dirs()

    adatas = []
    for key in args.keys:
        fp = args.input_dir / P.fn.filtered_normalized_adata.format(key=key)
        print(f"Reading {fp}")
        adatas.append(sc.read_h5ad(fp))

    adata = ad.concat(adatas, axis=0, join="outer",
                      label="batch", keys=args.keys)
    print(f"Concatenated: {adata.n_obs} cells from {len(args.keys)} slides")

    adata.obs["core_id"] = adata.obs["core_id"].replace(CORE_ID_RENAMES)

    core_anno = pd.read_csv(args.core_anno)
    valid_cores = set(core_anno["core_id"])
    keep_mask = adata.obs["core_id"].isin(valid_cores)
    n_removed = (~keep_mask).sum()
    n_cores_removed = adata.obs.loc[~keep_mask, "core_id"].nunique()
    adata = adata[keep_mask].copy()
    print(f"Removed {n_removed} cells from {n_cores_removed} excluded cores")

    adata.obs["patient_id"] = adata.obs["patient_id"].str.replace(r"\d+$", "", regex=True)

    merge_cols = ["core_id", "tissue_type", "tissue_type_dysplasia",
                  "mixed_core", "mixed_core_tissue_type", "mixed_core_dysplasia",
                  "mix_type", "distant_normal", "age", "gender", "colon_location"]
    anno = core_anno[merge_cols].copy()

    for col in ["mixed_core", "mixed_core_tissue_type", "mixed_core_dysplasia", "distant_normal"]:
        anno[col] = anno[col].fillna("").astype(str).str.upper() == "TRUE"

    obs = adata.obs.drop(
        columns=[c for c in merge_cols if c in adata.obs.columns and c != "core_id"],
        errors="ignore",
    )
    obs = obs.merge(anno, on="core_id", how="left")
    obs.index = adata.obs.index
    adata.obs = obs

    print(f"\nFinal: {adata.n_obs} cells, {adata.obs['core_id'].nunique()} cores, "
          f"{adata.obs['patient_id'].nunique()} patients")
    print(f"\ntissue_type:\n{adata.obs['tissue_type'].value_counts().to_string()}")
    print(f"\ntissue_type_dysplasia:\n{adata.obs['tissue_type_dysplasia'].value_counts().to_string()}")
    print(f"\nmixed_core: {adata.obs['mixed_core'].sum()} cells in "
          f"{adata.obs[adata.obs['mixed_core']]['core_id'].nunique()} cores")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
