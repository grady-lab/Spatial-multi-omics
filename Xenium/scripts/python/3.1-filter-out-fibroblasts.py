#!/usr/bin/env python3
"""Filter misannotated fibroblasts from per-tissue-type epithelial adatas
and write the fibroblast cell-id list for downstream correction.

Usage:
    python 3.1-filter-out-fibroblasts.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import P, ensure_dirs

import scanpy as sc
import anndata as ad


def validate_adata(adata, name):
    """Check for minority labels and duplicate cell_ids."""
    epi_obs = adata.obs[adata.obs["annotation_final_coarse"] == "Epithelial"]

    for col in ["tissue_type_cell_level", "tissue_type_dysplasia_cell_level"]:
        print(f"\n=== {name} / {col} ===")
        counts = epi_obs.groupby(["core_id", col]).size().reset_index(name="n_cells")
        core_label_counts = counts.groupby("core_id").filter(lambda x: len(x) > 1)
        small_labels = core_label_counts[core_label_counts["n_cells"] < 10]
        if len(small_labels) > 0:
            print(f"WARNING: cores with a label < 10 cells:")
            print(small_labels.to_string(index=False))
        else:
            print("All labels in mixed cores have >= 10 cells.")

    cell_ids = adata.obs["cell_id"].astype(str)
    n_dupes = len(cell_ids) - cell_ids.nunique()
    if n_dupes > 0:
        dupe_ids = cell_ids[cell_ids.duplicated(keep=False)].unique().tolist()
        print(f"WARNING: {n_dupes} duplicate cell_ids found: {dupe_ids}")
    else:
        print(f"No duplicate cell_ids ({len(cell_ids)} cells, all unique).")


def main():
    ensure_dirs()
    epi_dir = P.processed.adata.epithelial

    adjacent_normal = sc.read_h5ad(epi_dir / "epithelial_adjacent_normal_30_50.h5ad")
    separate_normal = sc.read_h5ad(epi_dir / "epithelial_distant_normal_30_50.h5ad")
    ta = sc.read_h5ad(epi_dir / "epithelial_AD_30_50_tt_harmony_batch_05_pt_05_ssg.h5ad")
    ca = sc.read_h5ad(epi_dir / "epithelial_CA_30_50_tt_harmony_batch_05_pt_05_ssg.h5ad")

    print(f"Adjacent normal: {adjacent_normal.n_obs} | Distant normal: {separate_normal.n_obs}")
    print(f"AD: {ta.n_obs} | CA: {ca.n_obs}")

    validate_adata(adjacent_normal, "adjacent_normal")
    validate_adata(separate_normal, "distant_normal")
    validate_adata(ta, "AD")
    validate_adata(ca, "CA")

    sep_filt = separate_normal[separate_normal.obs["distant_normal_0.8"] != "6"].copy()
    adj_filt = adjacent_normal[adjacent_normal.obs["adjacent_normal_0.8"] != "6"].copy()

    sep_fibro = separate_normal[separate_normal.obs["distant_normal_0.8"] == "6"].copy()
    adj_fibro = adjacent_normal[adjacent_normal.obs["adjacent_normal_0.8"] == "6"].copy()

    print(f"\nDistant normal: {separate_normal.n_obs} -> {sep_filt.n_obs} (removed {sep_fibro.n_obs} fibroblasts)")
    print(f"Adjacent normal: {adjacent_normal.n_obs} -> {adj_filt.n_obs} (removed {adj_fibro.n_obs} fibroblasts)")
    assert sep_filt.n_obs + sep_fibro.n_obs == separate_normal.n_obs
    assert adj_filt.n_obs + adj_fibro.n_obs == adjacent_normal.n_obs

    sep_filt.write_h5ad(epi_dir / "epithelial_distant_normal_filt.h5ad")
    adj_filt.write_h5ad(epi_dir / "epithelial_adjacent_normal_filt.h5ad")
    print(f"Wrote filtered adatas to {epi_dir}")

    combined_fibro = ad.concat([sep_fibro, adj_fibro], join="outer")
    fibro_out = P.annotations / "misannotated_fibroblasts.txt"
    with open(fibro_out, "w") as f:
        for cell_id in combined_fibro.obs["cell_id"]:
            f.write(f"{cell_id}\n")
    print(f"Wrote {len(combined_fibro)} fibroblast cell_ids to {fibro_out}")


if __name__ == "__main__":
    main()
