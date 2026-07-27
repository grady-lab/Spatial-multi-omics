#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs
from cluster_utils import cluster_xenium

import spatialdata as sd
import scanpy as sc


def get_args():
        parser = argparse.ArgumentParser("files!")
        parser.add_argument("-a", "--nonepithelial", help="xenium filtered normalized h5ad file",
                            default=str(P.processed.adata.nonepithelial / "nonepithelial_adata.h5ad"))
        parser.add_argument("-o", "--adata_output_dir", help="path for adata output",
                            default=str(P.processed.adata.nonepithelial))
        parser.add_argument("-s", "--scale", help="whether to scale", default="No")
        parser.add_argument("-v", "--hvgs", help="whether to use hvgs", default="No")
        parser.add_argument("-p", "--pcs", help="number of pcs", required=True)
        parser.add_argument("-n", "--nns", help="number of ns", required=True)
        return parser.parse_args()

args = get_args()
nonepi_fp = str(args.nonepithelial)
scale = str(args.scale)
hvgs = str(args.hvgs)
pcs = int(str(args.pcs))
nns = int(str(args.nns))

adata_output_dir = str(args.adata_output_dir)

ensure_dirs()


adata_nonepi = sc.read_h5ad(nonepi_fp)
cluster_xenium(adata_nonepi.copy(), nns, pcs, key="nonepi",output_path=f"{adata_output_dir}/nonepi_{nns}_{pcs}_harmonized_clustered.h5ad", scale=scale, hvgs=hvgs, harmony=["batch","patient_id"],thetas=[0.5,0.5])
