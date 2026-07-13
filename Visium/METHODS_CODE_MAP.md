# Methods-to-code map

| Manuscript method | Primary implementation |
|---|---|
| Spot RNA and ADT quality control | `R/preprocessing.R`: `load_visium_samples()`, `qc_rna()`, `prepare_adt()` |
| RNA normalization, FastMNN, ADT DSB/Harmony, WNN | `R/preprocessing.R`: `build_multimodal_reference()` |
| RNA-only label and embedding transfer | `R/preprocessing.R`: `map_rna_only_spots()` |
| Cluster markers and annotation export | `Rmd/My_Visium.Rmd`, cluster annotation section |
| stLearn spatial pseudotime and tissue islands | Independent Python workflow: `python/run_pseudotime.py` |
| Diffusion map | `R/scoring.R`: `add_diffusion_map()` |
| Seurat and ccAFv2 cell-cycle scoring | `R/scoring.R`: `add_cell_cycle_scores()` |
| CytoTRACE 2 | `R/scoring.R`: `run_cytotrace2()` |
| irGSEA/RRA enrichment | `R/scoring.R`: `run_irgsea()` |
| SenePy and cohort-wide outliers | `R/Cal_senepy.R`: `Calculate_senepy()`; `R/scoring.R`: `classify_senepy_outliers()` |
| Differential expression | `Rmd/My_Visium.Rmd`, differential-expression section |
| Bivariate local Moran's I | `R/spatial_statistics.R`: `bivariate_local_moran()` |
| Cell2location | Independent Python workflow: `python/run_cell2location.py` |
| GDF15-adjusted spatial lag correlation | `R/spatial_statistics.R`: `spatial_lag_correlations()`, `compare_before_after()` |
| CD8 neighborhood analysis | `R/spatial_statistics.R`: `prepare_cd8_neighborhood_data()`, `cd8_density_by_source()` |
| Luminal-surface distance | `R/spatial_statistics.R`: `distance_to_surface()` |
| Continuous and compositional statistics | `R/spatial_statistics.R`: `compare_sample_summaries()`, `composition_tests()` |
| Public scRNA-seq validation | `Rmd/My_Visium.Rmd`, validation section |
| Main and supplementary figures | `Rmd/Clean_plotting.Rmd`, `Rmd/Supple_plotting.Rmd` |
