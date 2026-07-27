#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import spatialdata as sd
import sopa
from spatialdata_io import xenium
import scanpy as sc
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
import os


def get_args():
    parser = argparse.ArgumentParser("Preprocess a single Xenium dataset")
    parser.add_argument("-k", "--key", required=True,
                        help="slide key (e.g. A1) — must exist in config.local.yaml `slides:` unless -i is given")
    parser.add_argument("-i", "--xenium_dir", type=Path, default=None,
                        help="raw Xenium output dir (default: config.local.yaml slides[<key>])")
    parser.add_argument("-c", "--core_cell_fp", type=Path, default=None,
                        help=f"core/cell CSV (default: {P.metadata}/slide_<key>_cells_and_cores.csv)")
    parser.add_argument("-s", "--sdata_output_dir", type=Path, default=None,
                        help=f"sdata output dir (default: {P.interim.sdata})")
    parser.add_argument("-a", "--adata_output_dir", type=Path, default=None,
                        help=f"adata output dir (default: {P.interim.adata}; also receives the cell-id list)")
    args = parser.parse_args()
    if args.xenium_dir is None:
        if args.key not in P.slides:
            parser.error(f"slide key '{args.key}' not found in config.local.yaml; pass -i explicitly")
        args.xenium_dir = Path(P.slides[args.key])
    if args.core_cell_fp is None:
        args.core_cell_fp = P.metadata / f"slide_{args.key}_cells_and_cores.csv"
    if args.sdata_output_dir is None:
        args.sdata_output_dir = P.interim.sdata
    if args.adata_output_dir is None:
        args.adata_output_dir = P.interim.adata
    return args

# ╭────────────────────────────────────────────────────────╮
# │                         CONFIG                         │
# ╰────────────────────────────────────────────────────────╯

args = get_args()
xenium_dir = args.xenium_dir
key = args.key
core_cell_fp = args.core_cell_fp
sdata_output_dir = args.sdata_output_dir
adata_output_dir = args.adata_output_dir
ensure_dirs()


# ╭────────────────────────────────────────────────────────╮
# │ 1: Create Xenium SpatialData objects, filter, and save │
# ╰────────────────────────────────────────────────────────╯
    
def read_core_dict(filename, cores_to_remove=None):
    """Read CSV and return dictionary with core_ids as keys and lists of cell_ids as values"""
    df = pd.read_csv(filename)
    # Group by core_id and convert cell_ids to lists
    core_dict = df.groupby('core_id')['cell_id'].apply(list).to_dict()
    if cores_to_remove:
        for core in cores_to_remove:
            core_dict.pop(core, None)  
    return core_dict

def assign_core_ids_to_cells(sdata, core_dict):
    """
    Assign core IDs to cells and remove cells without core assignments
    """
    
    # Create reverse mapping: cell_id -> core_id
    cell_to_core = {}
    for core_id, cell_list in core_dict.items():
        for cell_id in cell_list:
            cell_to_core[cell_id] = core_id
    
    sdata["table"].obs['core_id'] = sdata["table"].obs['cell_id'].map(cell_to_core)
    
    # Count cells before filtering
    cells_before = sdata["table"].n_obs
    
    mask = sdata["table"].obs['core_id'].notna()
    adata_filtered = sdata["table"][mask].copy()
    
    # Add patient_id and tissue_type columns after filtering
    adata_filtered.obs['patient_id'] = adata_filtered.obs['core_id'].str.split('-').str[0]
    adata_filtered.obs['tissue_type'] = adata_filtered.obs['core_id'].str.split('-').str[1]
    
    # Update the spatialdata object with filtered data
    sdata["table"] = adata_filtered
    
    # Print summary
    cells_removed = cells_before - adata_filtered.n_obs
    print(f"Cells before core assignment: {cells_before}")
    print(f"Cells after core assignment: {adata_filtered.n_obs}")
    print(f"Cells removed (no core assignment): {cells_removed}")
    
    return sdata

def aggregate_dapi(sdata, plot_output_dir=None):
    
    sdata.images["morphology_focus"] = sdata.images["morphology_focus"].sel(c=["DAPI"])
    
    modes = ['average']
    for mode in modes:
        intensities = sopa.aggregation.aggregate_channels(
            sdata,
            image_key="morphology_focus",
            shapes_key='cell_boundaries',
            mode=mode
        )
        sdata.tables['table'].obsm[f'intensities_{mode}'] = intensities
        sdata.tables['table'].obs[f'DAPI_{mode}'] = intensities[:, 0]
    
    # Plot histograms of DAPI average intensity
    if plot_output_dir is not None:
        os.makedirs(plot_output_dir, exist_ok=True)
        
        dapi_values = sdata.tables['table'].obs['DAPI_average']
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(dapi_values, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('DAPI Average Intensity', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('DAPI Intensity Distribution (50 bins)', fontsize=14)
        ax.grid(alpha=0.3)
        plot_path_50 = os.path.join(plot_output_dir, 'dapi_histogram_50bins.png')
        plt.savefig(plot_path_50, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved histogram (50 bins) to {plot_path_50}")
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(dapi_values, bins=100, edgecolor='black', alpha=0.7)
        ax.set_xlabel('DAPI Average Intensity', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('DAPI Intensity Distribution (100 bins)', fontsize=14)
        ax.grid(alpha=0.3)
        plot_path_100 = os.path.join(plot_output_dir, 'dapi_histogram_100bins.png')
        plt.savefig(plot_path_100, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved histogram (100 bins) to {plot_path_100}")
    
    
    return sdata


CORE_ID_RENAMES = {
    "Y2-CA-2": "Y2-TA-3",
    "Q1-CA-2": "Q1-TA-3",
    "O1-CA-1": "O1-TA-2",
    "Q1-CA-1": "Q1-TA-4",
    "L1-CA-2": "L1-TA-4",
}

sdata = xenium(xenium_dir)
core_dict = read_core_dict(core_cell_fp)
sdata = aggregate_dapi(
    sdata,
    plot_output_dir="./plots",
)
sdata = assign_core_ids_to_cells(sdata, core_dict)

table_obs = sdata["table"].obs
original = table_obs["core_id"].copy()
table_obs["core_id"] = table_obs["core_id"].replace(CORE_ID_RENAMES)
table_obs["tissue_type"] = table_obs["core_id"].str.split("-").str[1]
n_renamed = (original != table_obs["core_id"]).sum()
if n_renamed > 0:
    print(f"Renamed {n_renamed} cells across {len(CORE_ID_RENAMES)} core_id corrections")
sdata.write(f"{sdata_output_dir}/xenium_{key}.zarr", overwrite=True)



# ╭────────────────────────────────────────────────────────╮
# │                    2: Filter Xenium                    │
# ╰────────────────────────────────────────────────────────╯

def filter_sdata(sdata, sdata_save_path=None, adata_save_path=None, 
                 cell_ids_save_path=None):
    """
    Filter a single SpatialData object.
    """
    import os
    
    adata = sdata["table"].copy()  
    initial_cells = adata.n_obs
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # QC filtering
    sc.pp.filter_cells(adata, min_counts=15)
    sc.pp.filter_cells(adata, min_genes=3)
    sc.pp.filter_genes(adata, min_cells=5)
    
    cells_after_qc = adata.n_obs  # Save count after QC filtering

    # DAPI filtering
    dapi_p1 = np.percentile(adata.obs['DAPI_average'], 1)
    filter_mask = (adata.obs['DAPI_average'] > dapi_p1)
    adata_filtered = adata[filter_mask].copy()
    
    # Update the spatialdata object with filtered data
    sdata["table"] = adata_filtered
    
    print(f"Cells before QC filtering: {initial_cells}")
    print(f"Cells after QC filtering: {cells_after_qc}")
    print(f"Cells after DAPI filtering: {adata_filtered.n_obs}")
    print(f"Total removed: {initial_cells - adata_filtered.n_obs}")
    print(f"Genes: {adata_filtered.n_vars}")
    
    # Save filtered SpatialData object
    if sdata_save_path is not None:
        print(f"\nSaving filtered SpatialData to {sdata_save_path}...")
        sdata.write(sdata_save_path, overwrite=True)
        print("✓ Saved")
    
    # Save filtered AnnData
    if adata_save_path is not None:
        print(f"\nSaving filtered AnnData to {adata_save_path}...")
        adata_filtered.write_h5ad(adata_save_path)
        print("✓ Saved")
    
    # Save list of all cell IDs to txt file
    if cell_ids_save_path is not None:
        print(f"\nSaving cell IDs to {cell_ids_save_path}...")
        with open(cell_ids_save_path, 'w') as f:
            for cell_id in adata_filtered.obs_names:
                f.write(f"{cell_id}\n")
        print(f"✓ Saved {len(adata_filtered.obs_names)} cell IDs")
    
    return sdata, adata_filtered

sdata, adata_filtered = filter_sdata(
    sdata,
    sdata_save_path=sdata_output_dir / P.fn.filtered_sdata.format(key=key),
    adata_save_path=adata_output_dir / P.fn.filtered_adata.format(key=key),
    cell_ids_save_path=adata_output_dir / f"filtered_xenium_{key}_cell_ids.txt",
)




# ╭────────────────────────────────────────────────────────╮
# │                  3: Normalize Xenium                   │
# ╰────────────────────────────────────────────────────────╯

def normalize_xenium(adata, norm_target=100, log1p=True):
    # Save raw counts to a layer before normalizing
    adata.layers['counts'] = adata.X.copy()
    
    # Normalize the main X matrix
    sc.pp.normalize_total(adata, target_sum=norm_target)
    if log1p:
        sc.pp.log1p(adata)
    
    # Save normalized counts to a layer
    adata.layers['normalized'] = adata.X.copy()
    
    return adata

adata_filtered_normalized = normalize_xenium(adata_filtered)

del sdata.tables["table"]

# Assign the new one
sdata["table"] = adata_filtered_normalized

sdata.write(sdata_output_dir / f"filtered_normalized_xenium_{key}.zarr", overwrite=True)
adata_filtered_normalized.write_h5ad(adata_output_dir / P.fn.filtered_normalized_adata.format(key=key))


