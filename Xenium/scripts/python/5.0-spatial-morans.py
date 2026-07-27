#!/usr/bin/env python3
"""
Spatial Dense-Cell Annotation Pipeline
=======================================

Generalized pipeline that, for each user-specified feature (gene or .obs column),
performs:
  1. Epithelial region restriction (boundary segmentation)
  2. Local Moran's I spatial autocorrelation within epithelial regions
  3. "Dense" cell labeling (High-High + Low-High LISA quadrants)

Also supports pre-computed gene signature scores (e.g. TNFα/NF-κB pathway)
that are scored on epithelial cells before spatial analysis.

The final output is the input AnnData written to file with new .obs columns:
  - `{feature}_dense`           : bool, True if cell is in a LISA dense quadrant
  - `{feature}_morans_quad`     : str, LISA quadrant label
  - `{feature}_morans_pval`     : float, local Moran's I p-value
  - `{feature}_morans_Ii`       : float, local Moran's I statistic

Usage
-----
Edit the CONFIG section below, then run:
    python spatial_dense_annotation.py
"""

import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Union, Dict

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import spatialdata as sd

from scipy.spatial import cKDTree
from shapely.geometry import Point, Polygon as ShapelyPolygon
from shapely.ops import unary_union
from esda.moran import Moran_Local
from libpysal.weights import DistanceBand

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

warnings.simplefilter("ignore")

ensure_dirs()



# -- Paths -------------------------------------------------------------------
ADATA_PATH  = str(P.processed.adata.all_cells / P.fn.all_cells_final)
SDATA_DIR   = str(P.interim.sdata)
OUTPUT_PATH = str(P.processed.adata.all_cells / P.fn.all_cells_morans)


FEATURES: List[str] = [
    "GDF15", "OLFM4",
    "stemness_score", "senepy_intestine_epi_0",
]

# -- Gene signature scores ---------------------------------------------------

GENE_SCORES: Dict[str, List[str]] = {
    # "TNFA_NFKB_score": [
    #     'ABCA1','AREG','ATF3','ATP2B1','B4GALT1','B4GALT5','BCL2A1','BCL3','BCL6','BHLHE40',
    #     'BIRC2','BIRC3','BMP2','BTG1','BTG2','BTG3','CCL2','CCL20','CCL4','CCL5',
    #     'CCND1','CCNL1','CCRL2','CD44','CD69','CD80','CD83','CDKN1A','CEBPB','CEBPD',
    #     'CFLAR','CLCF1','CSF1','CSF2','CXCL1','CXCL10','CXCL11','CXCL2','CXCL3','CXCL6',
    #     'ACKR3','CCN1','RIGI','DENND5A','DNAJB4','DRAM1','DUSP1','DUSP2','DUSP4','DUSP5',
    #     'EDN1','EFNA1','EGR1','EGR2','EGR3','EHD1','EIF1','ETS2','F2RL1','F3',
    #     'FJX1','FOS','FOSB','FOSL1','FOSL2','FUT4','G0S2','GADD45A','GADD45B','GCH1',
    #     'GEM','GFPT2','GPR183','HBEGF','HES1','ICAM1','ICOSLG','ID2','IER2','IER3',
    #     'IER5','IFIH1','IFIT2','IFNGR2','IL12B','IL15RA','IL18','IL1A','IL1B','IL23A',
    #     'IL6','IL6ST','IL7R','INHBA','IRF1','IRS2','JAG1','JUN','JUNB','KDM6B',
    #     'KLF10','KLF2','KLF4','KLF6','KLF9','KYNU','LAMB3','LDLR','LIF','LITAF',
    #     'MAFF','MAP2K3','MAP3K8','MARCKS','MCL1','MSC','MXD1','MYC','NAMPT','NFAT5',
    #     'NFE2L2','NFIL3','NFKB1','NFKB2','NFKBIA','NFKBIE','NINJ1','NR4A1','NR4A2','NR4A3',
    #     'OLR1','PANX1','PDE4B','PDLIM5','PER1','PFKFB3','PHLDA1','PHLDA2','PLAU','PLAUR',
    #     'PLEK','PLK2','PMEPA1','PNRC1','PLPP3','PPP1R15A','PTGER4','PTGS2','PTPRE','PTX3',
    #     'RCAN1','REL','RELA','RELB','RHOB','RIPK2','RNF19B','SAT1','SDC4','SERPINB2',
    #     'SERPINB8','SERPINE1','SGK1','SIK1','SLC16A6','SLC2A3','SLC2A6','SMAD3','SNN','SOCS3',
    #     'SOD2','SPHK1','SPSB1','SQSTM1','STAT5A','TANK','TAP1','TGIF1','TIPARP','TLR2',
    #     'TNC','TNF','TNFAIP2','TNFAIP3','TNFAIP6','TNFAIP8','TNFRSF9','TNFSF9','TNIP1','TNIP2',
    #     'TRAF1','TRIB1','TRIP10','TSC22D1','TUBB2A','VEGFA','YRDC','ZBTB10','ZC3H12A','ZFP36',
    # ],
}

EPITHELIAL_COL       = "annotation_final_fine_cd8"   # obs column identifying epithelial cells
EPITHELIAL_LABEL     = "Epithelial"                  # value in that column
MIN_EPI_CELLS        = 100       # minimum epithelial cells per core to process
NEIGHBOR_RADIUS      = 50        # radius for epithelial neighbor density filtering
MIN_EPI_NEIGHBORS    = 2         # minimum neighbors to keep a cell
BUFFER_DISTANCE      = 45        # buffer around cells for mask construction
MIN_MASK_AREA        = 200       # minimum mask polygon area
MIN_CELLS_PER_MASK   = 100       # minimum cells per mask region
SIMPLIFY_TOLERANCE   = 10        # boundary simplification tolerance

MORANS_DISTANCE      = 50        # distance threshold for spatial weights
MORANS_PERMUTATIONS  = 999       # permutations for significance testing
MORANS_SEED          = 31        # random seed 
MIN_CELLS_FOR_MORANS = 30        # minimum epithelial cells in mask for Moran's I

DENSE_QUADRANTS      = ["High-High", "Low-High"]  # LISA quadrants treated as "dense"


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PIPELINE FUNCTIONS                                                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def load_spatialdata(sdata_dir: str) -> dict:
    """Load all filtered/normalized SpatialData .zarr files into a dictionary."""
    spatialdata_dict = {}
    entries = [f for f in os.listdir(sdata_dir) if f.startswith("fil")]
    for filename in sorted(entries):
        sdata = sd.read_zarr(os.path.join(sdata_dir, filename))
        key = filename.replace("filtered_normalized_xenium_", "").replace(".zarr", "")
        spatialdata_dict[key] = sdata
    print(f"Loaded {len(spatialdata_dict)} SpatialData sections")
    return spatialdata_dict



def filter_valid_cores(
    spatialdata_dict: dict,
    all_cells_adata: ad.AnnData,
    epithelial_col: str,
    epithelial_label: str,
    min_epi_cells: int,
) -> dict:
    """Filter spatialdata_dict to (section, core) pairs with sufficient epithelial cells."""
    epi_obs = all_cells_adata.obs[all_cells_adata.obs[epithelial_col] == epithelial_label]
    epi_counts = epi_obs.groupby(["batch", "core_id"]).size()
    valid_pairs = epi_counts[epi_counts >= min_epi_cells].reset_index()
    valid_section_cores = set(zip(valid_pairs["batch"], valid_pairs["core_id"]))
    valid_cell_ids = set(all_cells_adata.obs["cell_id"].unique())

    filtered = {}
    for section_key, sdata in spatialdata_dict.items():
        table = sdata.tables["table"]
        valid_cores_for_section = {
            core for (batch, core) in valid_section_cores if batch == section_key
        }
        mask = (
            table.obs["cell_id"].isin(valid_cell_ids)
            & table.obs["core_id"].isin(valid_cores_for_section)
        )
        if mask.sum() > 0:
            sdata.tables["table"] = table[mask].copy()
            filtered[section_key] = sdata
            n_cores = sdata.tables["table"].obs["core_id"].nunique()
            print(f"  {section_key}: {mask.sum()} cells, {n_cores} cores")

    print(f"Valid (section, core) pairs: {len(valid_section_cores)}")
    return filtered


def transfer_obs_columns(
    spatialdata_dict: dict, reference_adata: ad.AnnData, columns: list
) -> None:
    """Transfer annotation columns from reference AnnData to SpatialData tables."""
    existing_cols = [c for c in columns if c in reference_adata.obs.columns]
    ref_obs = reference_adata.obs.set_index("cell_id")[existing_cols]
    for section_key, sdata in spatialdata_dict.items():
        table = sdata.tables["table"]
        for col in existing_cols:
            table.obs[col] = table.obs["cell_id"].map(ref_obs[col])


def compute_gene_scores(
    all_cells_adata: ad.AnnData,
    gene_scores: Dict[str, List[str]],
    epithelial_col: str,
    epithelial_label: str,
) -> None:
    """
    Compute gene signature scores on epithelial cells and store in all_cells_adata.obs.

    For each score_name → gene_list pair:
      1. Subset to epithelial cells
      2. Filter gene_list to genes present in .var_names
      3. Run sc.tl.score_genes on the epithelial subset
      4. Map the score back to all_cells_adata.obs (NaN for non-epithelial cells)

    Modifies all_cells_adata.obs in-place.
    """
    epi_mask = all_cells_adata.obs[epithelial_col] == epithelial_label
    epi_adata = all_cells_adata[epi_mask].copy()
    print(f"   Epithelial subset: {len(epi_adata)} cells")

    for score_name, gene_list in gene_scores.items():
        
        valid_genes = [g for g in gene_list if g in epi_adata.var_names]
        n_total = len(gene_list)
        n_found = len(valid_genes)
        print(f"   Scoring '{score_name}': {n_found}/{n_total} genes found in panel")

        if n_found == 0:
            print(f"   WARNING: No genes found for '{score_name}'. Skipping.")
            all_cells_adata.obs[score_name] = np.nan
            continue

        
        sc.tl.score_genes(epi_adata, gene_list=valid_genes, score_name=score_name)

        
        all_cells_adata.obs[score_name] = np.nan
        score_map = epi_adata.obs.set_index(
            epi_adata.obs["cell_id"] if "cell_id" in epi_adata.obs.columns
            else epi_adata.obs.index
        )[score_name]
        all_cells_adata.obs.loc[epi_mask, score_name] = (
            all_cells_adata.obs.loc[epi_mask, "cell_id"].map(score_map).values
        )

        scored = all_cells_adata.obs[score_name].notna().sum()
        print(f"   → '{score_name}' scored for {scored} epithelial cells "
              f"(mean={all_cells_adata.obs[score_name].mean():.4f})")


def create_epithelial_boundaries(
    filtered_spatialdata_dict: dict,
    core_ids: list,
    epithelial_col: str = "annotation_final_fine_cd8",
    epithelial_label: str = "Epithelial",
    neighbor_radius: float = 50,
    min_epi_neighbors: int = 3,
    buffer_distance: float = 100,
    min_mask_area: float = 10000,
    min_cells_per_mask: int = 50,
    simplify_tolerance: float = 10,
    use_cell_boundaries: bool = True,
) -> dict:
    """
    Create epithelial region masks via density-based filtering and buffering.

    Returns
    -------
    dict : {section_key: {core_id: {'masks': [Polygon, ...], ...}}}
    """
    results = {}

    for section_key, sdata in filtered_spatialdata_dict.items():
        results[section_key] = {}
        table = sdata.tables["table"]
        core_mask = table.obs["core_id"].isin(core_ids)
        cores_in_section = table.obs.loc[core_mask, "core_id"].unique()

        for core_id in cores_in_section:
            core_cell_mask = table.obs["core_id"] == core_id
            epi_mask = core_cell_mask & (table.obs[epithelial_col] == epithelial_label)

            if epi_mask.sum() < min_cells_per_mask:
                continue

            epi_cells = table.obs[epi_mask]
            epi_cell_ids = epi_cells["cell_id"].values

            shapes_key = "cell_boundaries" if use_cell_boundaries else "cell_circles"
            cell_shapes = sdata.shapes[shapes_key]
            epi_spatial = cell_shapes[cell_shapes.index.isin(epi_cell_ids)]

            if len(epi_spatial) == 0:
                continue

            
            centroids = np.array(
                [[g.centroid.x, g.centroid.y] for g in epi_spatial.geometry]
            )
            centroid_ids = epi_spatial.index.values
            tree = cKDTree(centroids)
            neighbor_counts = np.array(
                [len(nbrs) - 1 for nbrs in tree.query_ball_point(centroids, neighbor_radius)]
            )
            dense_mask = neighbor_counts >= min_epi_neighbors
            filtered_cell_ids = centroid_ids[dense_mask]
            filtered_spatial = epi_spatial[epi_spatial.index.isin(filtered_cell_ids)]

            if len(filtered_spatial) < min_cells_per_mask:
                continue

           
            buffered_cells = []
            for geom in filtered_spatial.geometry:
                try:
                    b = (
                        geom.buffer(buffer_distance)
                        if use_cell_boundaries
                        else Point(geom.centroid.x, geom.centroid.y).buffer(buffer_distance)
                    )
                    if b.is_valid and not b.is_empty:
                        buffered_cells.append(b)
                except Exception:
                    continue

            if not buffered_cells:
                continue

            
            merged = unary_union(buffered_cells)
            mask_polygons = (
                [merged]
                if merged.geom_type == "Polygon"
                else list(merged.geoms) if merged.geom_type == "MultiPolygon"
                else []
            )

            
            valid_masks, mask_stats = [], []
            for mask in mask_polygons:
                area = mask.area
                cells_in = sum(
                    1
                    for g in filtered_spatial.geometry
                    if (
                        mask.intersects(g)
                        if use_cell_boundaries
                        else mask.contains(Point(g.centroid.x, g.centroid.y))
                    )
                )
                if area >= min_mask_area and cells_in >= min_cells_per_mask:
                    if simplify_tolerance > 0:
                        mask = mask.simplify(simplify_tolerance, preserve_topology=True)
                    valid_masks.append(mask)
                    mask_stats.append({"area": area, "n_cells": cells_in})

            if not valid_masks:
                continue

            results[section_key][core_id] = {
                "masks": valid_masks,
                "n_masks": len(valid_masks),
                "mask_stats": mask_stats,
                "epi_cell_ids": epi_cell_ids,
                "filtered_cell_ids": filtered_cell_ids,
            }
            total_cells = sum(s["n_cells"] for s in mask_stats)
            print(f"  {section_key}/{core_id}: {len(valid_masks)} mask(s), {total_cells} cells")

    return results


def resolve_feature_values(
    adata_subset: ad.AnnData,
    feature: str,
    all_cells_adata: ad.AnnData,
) -> Optional[np.ndarray]:
    """
    Resolve a feature to a numeric array for the cells in adata_subset.

    If `feature` is a gene in var_names → extract from .X
    If `feature` is an .obs column     → extract from obs (must be numeric)

    Returns None if the feature cannot be resolved to numeric values.
    """
    if feature in adata_subset.var_names:
        vals = adata_subset[:, feature].X
        vals = vals.toarray().flatten() if hasattr(vals, "toarray") else np.asarray(vals).flatten()
        return vals.astype(np.float64)

    
    for source in [adata_subset, all_cells_adata]:
        if feature in source.obs.columns:
            
            if source is all_cells_adata and "cell_id" in adata_subset.obs.columns:
                cell_ids = adata_subset.obs["cell_id"].values
                ref = all_cells_adata.obs.set_index("cell_id")[feature]
                vals = pd.Series(cell_ids).map(ref).values
            else:
                vals = source.obs[feature].values

            
            vals = pd.to_numeric(vals, errors="coerce")
            if np.all(np.isnan(vals)):
                return None
            return vals.astype(np.float64)

    return None


def run_morans_i_epithelial(
    spatialdata_dict: dict,
    boundaries: dict,
    feature: str,
    all_cells_adata: ad.AnnData,
    epithelial_col: str = "annotation_final_fine",
    epithelial_label: str = "Epithelial",
    distance_threshold: float = 50,
    n_permutations: int = 999,
    seed: int = 31,
    min_cells: int = 30,
) -> pd.DataFrame:
    """
    Local Moran's I for a given feature on epithelial cells within boundary masks.

    Parameters
    ----------
    feature : str
        Gene name (from .var_names / .X) or .obs column name.

    Returns
    -------
    DataFrame with per-cell Moran's I statistics and LISA quadrant labels.
    """
    QUADRANT_LABELS = {1: "High-High", 2: "Low-High", 3: "Low-Low", 4: "High-Low"}
    all_results = []

    for sample_name, sdata in spatialdata_dict.items():
        if sample_name not in boundaries:
            continue

        adata = sdata.tables["table"]

        for core in adata.obs["core_id"].unique():
            if core not in boundaries[sample_name]:
                continue

            epi_masks = boundaries[sample_name][core]["masks"]
            if not epi_masks:
                continue

            
            adata_core = adata[adata.obs["core_id"] == core]
            epi_cell_mask = adata_core.obs[epithelial_col] == epithelial_label
            adata_epi = adata_core[epi_cell_mask].copy()

            if len(adata_epi) == 0:
                continue

            epi_coords = adata_epi.obsm["spatial"]
            epi_cell_ids = (
                adata_epi.obs["cell_id"].values
                if "cell_id" in adata_epi.obs.columns
                else adata_epi.obs.index.values
            )

            
            in_mask = np.zeros(len(epi_coords), dtype=bool)
            for i, (x, y) in enumerate(epi_coords):
                pt = Point(x, y)
                if any(m.contains(pt) for m in epi_masks):
                    in_mask[i] = True

            if in_mask.sum() < min_cells:
                continue

            adata_epi_zone = adata_epi[in_mask].copy()
            coords_zone = epi_coords[in_mask]
            ids_zone = epi_cell_ids[in_mask]

           
            feature_vals = resolve_feature_values(adata_epi_zone, feature, all_cells_adata)
            if feature_vals is None:
                print(f"  WARNING: Could not resolve feature '{feature}' for "
                      f"{sample_name}/{core}. Skipping.")
                continue

            
            feature_vals = np.nan_to_num(feature_vals, nan=0.0)

           
            if np.std(feature_vals) == 0:
                continue

            
            w = DistanceBand.from_array(coords_zone, threshold=distance_threshold, binary=True)
            if w.mean_neighbors < 1:
                continue

            np.random.seed(seed)
            moran = Moran_Local(feature_vals, w, permutations=n_permutations)

            quadrant_names = [
                QUADRANT_LABELS.get(q, "NS") if p < 0.05 else "Not significant"
                for q, p in zip(moran.q, moran.p_sim)
            ]

            core_df = pd.DataFrame({
                "sample": sample_name,
                "core": core,
                "cell_id": ids_zone,
                "x": coords_zone[:, 0],
                "y": coords_zone[:, 1],
                "feature_value": feature_vals,
                "local_morans_i": moran.Is,
                "z_score": moran.z_sim,
                "p_value": moran.p_sim,
                "significant": moran.p_sim < 0.05,
                "quadrant": quadrant_names,
                "quadrant_num": moran.q,
            })
            all_results.append(core_df)

            sig_pct = (moran.p_sim < 0.05).mean() * 100
            hh = sum(1 for q in quadrant_names if q == "High-High")
            print(f"  {sample_name}/{core}: {len(core_df)} cells, "
                  f"{sig_pct:.0f}% significant, {hh} High-High")

    if not all_results:
        print(f"  WARNING: No results for feature '{feature}'")
        return pd.DataFrame()

    return pd.concat(all_results, ignore_index=True)


def label_dense_cells(
    morans_df: pd.DataFrame,
    dense_quadrants: list,
) -> pd.DataFrame:
    """
    Label cells in dense regions based on LISA quadrant membership.

    Returns the input DataFrame with a new boolean column 'is_dense'.
    """
    result = morans_df.copy()
    result["is_dense"] = result["quadrant"].isin(dense_quadrants)

    n_dense = result["is_dense"].sum()
    print(f"  Dense cells: {n_dense}/{len(result)} ({n_dense / len(result) * 100:.1f}%)")
    return result


def annotate_cells_with_morans(
    all_cells_adata: ad.AnnData,
    morans_labeled: pd.DataFrame,
    feature: str,
    col_prefix: str = None,
) -> None:
    """
    Add per-cell Moran's I annotations to all_cells_adata.obs for a given feature.

    New columns added:
        {prefix}_dense       : bool — cell is in a LISA dense quadrant (High-High or Low-High)
        {prefix}_morans_quad : str  — LISA quadrant label
        {prefix}_morans_pval : float
        {prefix}_morans_Ii   : float
    """
    if col_prefix is None:
        col_prefix = feature

    
    safe_prefix = col_prefix.replace(" ", "_").replace("-", "_")

    
    all_cells_adata.obs[f"{safe_prefix}_dense"] = False
    all_cells_adata.obs[f"{safe_prefix}_morans_quad"] = pd.NA
    all_cells_adata.obs[f"{safe_prefix}_morans_pval"] = np.nan
    all_cells_adata.obs[f"{safe_prefix}_morans_Ii"] = np.nan

    if len(morans_labeled) == 0:
        return

    morans_indexed = morans_labeled.set_index("cell_id")
    cell_ids = all_cells_adata.obs["cell_id"]

    
    all_cells_adata.obs[f"{safe_prefix}_morans_quad"] = cell_ids.map(morans_indexed["quadrant"]).values
    all_cells_adata.obs[f"{safe_prefix}_morans_pval"] = cell_ids.map(morans_indexed["p_value"]).values
    all_cells_adata.obs[f"{safe_prefix}_morans_Ii"] = cell_ids.map(morans_indexed["local_morans_i"]).values
    all_cells_adata.obs[f"{safe_prefix}_dense"] = cell_ids.map(morans_indexed["is_dense"]).fillna(False).values

    
    n_dense = all_cells_adata.obs[f"{safe_prefix}_dense"].sum()
    n_analyzed = cell_ids.isin(morans_indexed.index).sum()
    print(f"  → {safe_prefix}: {n_dense} dense cells out of {n_analyzed} analyzed")


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  MAIN PIPELINE                                                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def main():
    print("=" * 70)
    print("1. Loading AnnData")
    print("=" * 70)
    all_cells_adata = sc.read_h5ad(ADATA_PATH)
    print(f"   {all_cells_adata.shape[0]} cells × {all_cells_adata.shape[1]} genes")

    before = len(all_cells_adata)
    all_cells_adata = all_cells_adata[
        ~all_cells_adata.obs["cell_id"].duplicated(keep=False)
    ].copy()
    print(f"   {before} → {len(all_cells_adata)} cells "
          f"({before - len(all_cells_adata)} duplicates removed)")


    if GENE_SCORES:
        print("\n" + "=" * 70)
        print("2. Computing gene signature scores (epithelial cells only)")
        print("=" * 70)
        compute_gene_scores(
            all_cells_adata, GENE_SCORES,
            epithelial_col=EPITHELIAL_COL,
            epithelial_label=EPITHELIAL_LABEL,
        )

   
    print("\n" + "=" * 70)
    print("3. Validating features")
    print("=" * 70)
    valid_features = []
    for feat in FEATURES:
        
        if not feat or not feat.strip():
            continue
        if feat in all_cells_adata.var_names:
            print(f"   ✓ '{feat}' found in .var_names (gene)")
            valid_features.append(feat)
        elif feat in all_cells_adata.obs.columns:
            vals = pd.to_numeric(all_cells_adata.obs[feat], errors="coerce")
            non_null = vals.notna().sum()
            print(f"   ✓ '{feat}' found in .obs ({non_null} non-null numeric values)")
            valid_features.append(feat)
        else:
            print(f"   ✗ '{feat}' NOT FOUND in .var_names or .obs — skipping")

    if not valid_features:
        raise ValueError("No valid features found. Check FEATURES list.")

    print("\n" + "=" * 70)
    print("4. Loading SpatialData")
    print("=" * 70)
    spatialdata_dict = load_spatialdata(SDATA_DIR)

    print("\n" + "=" * 70)
    print("5. Filtering valid cores")
    print("=" * 70)
    filtered_sdata = filter_valid_cores(
        spatialdata_dict, all_cells_adata,
        EPITHELIAL_COL, EPITHELIAL_LABEL, MIN_EPI_CELLS,
    )

    print("\n" + "=" * 70)
    print("6. Transferring annotations to SpatialData")
    print("=" * 70)
    transfer_cols = list(set(
        [EPITHELIAL_COL, "annotation_final_fine", "core_id", "patient_id"]
        + [f for f in valid_features if f in all_cells_adata.obs.columns]
    ))
    transfer_obs_columns(filtered_sdata, all_cells_adata, transfer_cols)
    print(f"   Transferred {len(transfer_cols)} columns")

    print("\n" + "=" * 70)
    print("7. Building epithelial boundaries")
    print("=" * 70)
    all_cores = sorted({
        core
        for sdata in filtered_sdata.values()
        for core in sdata.tables["table"].obs["core_id"].unique()
    })
    print(f"   Processing {len(all_cores)} unique cores")

    boundaries_all = create_epithelial_boundaries(
        filtered_sdata,
        core_ids=all_cores,
        epithelial_col=EPITHELIAL_COL,
        epithelial_label=EPITHELIAL_LABEL,
        neighbor_radius=NEIGHBOR_RADIUS,
        min_epi_neighbors=MIN_EPI_NEIGHBORS,
        buffer_distance=BUFFER_DISTANCE,
        min_mask_area=MIN_MASK_AREA,
        min_cells_per_mask=MIN_CELLS_PER_MASK,
        simplify_tolerance=SIMPLIFY_TOLERANCE,
        use_cell_boundaries=True,
    )
    n_boundaries = sum(len(cores) for cores in boundaries_all.values())
    print(f"\n   Built boundaries for {n_boundaries} cores")

    for feat_idx, feat in enumerate(valid_features, 1):
        print("\n" + "=" * 70)
        print(f"8.{feat_idx} Processing feature: {feat}")
        print("=" * 70)

        print(f"\n  a. Running Local Moran's I for '{feat}'")
        morans_df = run_morans_i_epithelial(
            filtered_sdata, boundaries_all,
            feature=feat,
            all_cells_adata=all_cells_adata,
            epithelial_col="annotation_final_fine",
            epithelial_label=EPITHELIAL_LABEL,
            distance_threshold=MORANS_DISTANCE,
            n_permutations=MORANS_PERMUTATIONS,
            seed=MORANS_SEED,
            min_cells=MIN_CELLS_FOR_MORANS,
        )

        if morans_df.empty:
            print(f"  No Moran's I results for '{feat}'. Skipping.")
            continue

        print(f"  Total cells analyzed: {len(morans_df)}")

        print(f"\n  b. Labeling dense cells for '{feat}'")
        morans_labeled = label_dense_cells(morans_df, DENSE_QUADRANTS)

        print(f"\n  c. Annotating cells for '{feat}'")
        annotate_cells_with_morans(
            all_cells_adata, morans_labeled, feature=feat,
        )

    print("\n" + "=" * 70)
    print("9. Saving annotated AnnData")
    print("=" * 70)

    all_cells_adata.obs_names_make_unique()

    suffixes = ["_dense", "_morans_quad", "_morans_pval", "_morans_Ii"]
    new_cols = [
        c for c in all_cells_adata.obs.columns
        if any(c.endswith(s) for s in suffixes)
    ]
    print(f"   New annotation columns: {new_cols}")

    score_cols = [c for c in all_cells_adata.obs.columns if c in GENE_SCORES]
    if score_cols:
        print(f"   Gene score columns: {score_cols}")

    output_path = OUTPUT_PATH
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    all_cells_adata.write_h5ad(output_path)
    print(f"   Written to: {output_path}")
    print(f"   Shape: {all_cells_adata.shape}")
    print("\nDone.")


if __name__ == "__main__":
    main()
