#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs
from cluster_utils import cluster_xenium

import scanpy as sc
import anndata as ad


RESES_SHORT = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0)


def get_args():
        parser = argparse.ArgumentParser("files!")
        parser.add_argument("-a", "--xenium_fp", help="xenium filtered normalized h5ad file",
                            default=str(P.processed.adata.all_cells / P.fn.combined_adata))
        parser.add_argument("-k", "--key", help="key", required=True)
        parser.add_argument("-o", "--adata_output_dir", help="path for adata output",
                            default=str(P.processed.adata.all_cells))
        parser.add_argument("-s", "--scale", help="whether to scale", default="No")
        parser.add_argument("-v", "--hvgs", help="whether to use hvgs", default="No")
        parser.add_argument("-p", "--pcs", help="number of pcs", required=True)
        parser.add_argument("-n", "--nns", help="number of ns", required=True)
        return parser.parse_args()

args = get_args()
xenium_fp= str(args.xenium_fp)
key = str(args.key)
scale = str(args.scale)
hvgs = str(args.hvgs)
pcs = int(str(args.pcs))
nns = int(str(args.nns))
adata_output_dir = str(args.adata_output_dir)

ensure_dirs()

xenium_adata = sc.read_h5ad(xenium_fp)

cluster_xenium(xenium_adata.copy(), nns, pcs, key=key, output_path=f"{adata_output_dir}/xenium_allcells_{key}_{nns}_{pcs}_harmonized_clustered.h5ad", scale=scale, hvgs=hvgs, reses=RESES_SHORT, harmony=["batch","patient_id"], thetas=[0.5,0.5])
