# Methods-to-code map

| Manuscript method | Primary implementation |
|---|---|
| Spot RNA and ADT quality control | `R/preprocessing.R`: `load_visium_samples()`, `qc_rna()`, `prepare_adt()` |
| RNA normalization, FastMNN, ADT DSB/Harmony, WNN | `R/preprocessing.R`: `build_multimodal_reference()` |
| RNA-only label and embedding transfer | `R/preprocessing.R`: `map_rna_only_spots()` |
| Cluster markers and annotation export | `Rmd/Visium_analyses.Rmd`, cluster annotation section |
| stLearn spatial pseudotime and tissue islands | Independent Python workflow: `python/run_pseudotime.py` |
| FastMNN-based epithelial diffusion map | `R/scoring.R`: `add_diffusion_map()`; `Rmd/Visium_analyses.Rmd`, `import-stlearn-output` section |
| Seurat and ccAFv2 cell-cycle scoring | `R/scoring.R`: `add_cell_cycle_scores()` |
| CytoTRACE 2 | `R/scoring.R`: `prepare_cytotrace2_input()`, `run_cytotrace2()`; clusters C0, C1, C2, C3, C4, and C6 only |
| irGSEA/RRA enrichment | `R/scoring.R`: `run_irgsea()` |
| SenePy and cohort-wide outliers | `R/Cal_senepy.R`: `Calculate_senepy()`; `R/scoring.R`: `classify_senepy_outliers()` |
| Differential expression | `Rmd/Visium_analyses.Rmd`, `differential-expression-gdf15` section |
| Bivariate local Moran's I | `R/spatial_statistics.R`: `bivar_local_moran_seurat()` |
| Public scRNA-seq reference curation and Cell2location | `R/scrna_reference.R`: `curate_marteau_reference()`, `validate_marteau_reference()`, `export_marteau_reference_h5ad()`; `Rmd/Visium_analyses.Rmd`, `build-scrna-reference` and `import-cell2location-output` sections (including reported aggregate abundances); independent Python workflow: `python/run_cell2location.py` |
| GDF15-adjusted spatial lag correlation | `R/spatial_statistics.R`: `make_listw_knn()`, `lag_spearman_perm_mode()`, `within_sample_lagcorr()`, `compare_before_after()` |
| CD8 spatial ring-density analysis | `R/spatial_statistics.R`: `prepare_cd8_neighborhood_data()`, `add_tissue_components()`, `cd8_density_by_source()`; `Rmd/Visium_analyses.Rmd`, `cd8-neighborhood` section |
| Adenoma radial distance from the tract surface | `R/spatial_statistics.R`: `build_staffli_from_spaceranger()`, `separate_tissue_islands_semla()`, `distance_to_surface()`, `summarize_distance_bins()`; `Rmd/Visium_analyses.Rmd`, `tract-surface-distance` section |
| Continuous and categorical-proportion statistics | `R/spatial_statistics.R`: `compare_score_clusters_paired()`, `run_comp_tests()`; `Rmd/Visium_analyses.Rmd`, `sample-level-statistics` section |
| Independent public scRNA-seq validation | `Rmd/Visium_analyses.Rmd`, `single-cell-validation` section |
| Class-faceted marker dot plot and three-way Cell2location spatial blend | `R/plotting.R`: `DotPlot_by_class_horizontal()`, `three_way_spatial_plot()`; `Rmd/Clean_plotting.Rmd` |
| Main and supplementary figures | `Rmd/Clean_plotting.Rmd`, `Rmd/Supple_plotting.Rmd` |
