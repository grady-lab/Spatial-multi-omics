# Validation performed

The release draft was checked on 2026-07-09 as follows:

- All R source files parse under R 4.4.0.
- All R chunks extracted from the three R Markdown notebooks parse under R 4.4.0.
- Both Python scripts pass Python 3 syntax compilation.
- Cell2location stages and defaults were checked against the authoritative `Data/Python/Deconv/run_cell2location.py`, including reference/spatial training, posterior export, NMF colocalization, and cell-type-specific expected-expression layers.
- Synthetic-data tests pass for the before/after GDF15 comparison, tissue connected components, CD8 neighborhood density, and luminal-surface distance.
- No analyst-specific absolute paths occur in the release directory.
- The sample manifest resolves to 22 included sections (19 RNA+ADT and 3 RNA-only).
- No source file in the release is larger than 1 MB; raw and derived scientific data types are ignored by Git.

Full numerical reproduction was not run because the raw data, reference files, and large derived objects are external inputs. HTML rendering was not tested in this compute session because Pandoc is unavailable; notebook code extraction and R parsing succeeded.

On 2026-07-20, synthetic-data parity tests confirmed that the public `compare_score_clusters_paired()` and `run_comp_tests()` implementations match the authoritative helper script in both paired spot-level and unpaired sample-level modes. The categorical-proportion parity test used the locally available vegan 2.6.4; the release environment remains pinned to vegan 2.7.2 in `DESCRIPTION`.

The Marteau reference curation was also checked against the realized Cell2location input H5AD. Starting from the archived `marteau_TA_updated.qs`, the public rules reproduce 52,871 cells, 24,593 genes, 66 samples from `Chen_2021_Cell`, all 28 annotation names, and every per-annotation cell count. The 9.7 GB H5AD itself was inspected without loading its dense expression matrix; a new full H5AD export was not written during validation.

The diffusion-map implementation was checked against `Scripts/DiffusionMap.R`: the same epithelial clusters, per-sample RNA-layer split, normalization/PCA/FastMNN sequence, MNN dimensions 1–20, `sigma = 1`, diffusion eigenvectors 1 and 3, and pseudotime-based DC1 sign orientation are explicit in the public workflow.

The CytoTRACE2 input preparation follows the active source chunk: only clusters C0, C1, C2, C3, C4, and C6 are retained, the joined `Spatial` assay is copied to `RNA`, `DietSeurat()` retains the RNA assay and WNN UMAP, and `cytotrace2()` receives raw counts with `species = "human"` and `is_seurat = TRUE`. Returned scores are matched back by barcode so excluded clusters remain missing rather than receiving recycled values.

The plotting and spatial-association code was audited against the four named authoritative utilities on 2026-07-20. `DotPlot_by_class_horizontal()` retains the source feature ordering, class facets, scaled-expression gradient, and percent-expression point sizes. The three-way Cell2location plot uses the source 99th-percentile cap, per-feature min-max scaling, subtractive RGB blend, colors, and legends; the two inconsistent variable names in the exploratory source script were corrected to the declared `CD8_T_cells` feature and the generated spatial plot. The released local Moran wrapper preserves metadata-first feature lookup, per-image kNN weights, `Not sig` labels, and target-specific metadata names. For spatial lag correlations, only the final definitions in `Spatial smoothed spatial correlation.R` were retained: symmetric row-standardized kNN weights, duplicate-coordinate jitter, `x_vs_Wy`/`Wx_vs_Wy` modes, within-sample residualization, and permutation tests. Synthetic numerical parity with the authoritative final definitions passed for unadjusted `x_vs_Wy`, GDF15-adjusted `x_vs_Wy` with both variables residualized, and unadjusted `Wx_vs_Wy`, including per-sample statistics and group-label permutation p-values. Source files and extracted notebook chunks parse under R 4.4.0. The Seurat/rgeoda wrapper and final plots were not executed end-to-end because the large derived Seurat object is not included in the repository and Seurat did not load successfully in the active compute environment.
