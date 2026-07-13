# Data input contract

Data files are not committed to this repository. Paths below are relative to the repository root.

## Raw Visium inputs

For every row marked `include=TRUE` in `config/sample_manifest.csv`, provide the corresponding Space Ranger `outs/` directory at `data/raw/<sample_id>/outs/`. The directory must contain `filtered_feature_bc_matrix.h5` and the `spatial/` subdirectory.

The current manifest contains 22 analyzed sections: 3 RNA-only sections and 19 sections with RNA plus ADT. One additional section failed Space Ranger QC before object construction and is consequently absent from the manifest. `Pilot_small_2` is retained in the manifest as an explicitly excluded section because it failed downstream Seurat QC.

## ADT feature mapping

The confirmed feature-to-reagent names used by the composite CD8-positive spot gate are recorded in `config/adt_feature_map.csv`. In particular, `ADT-PTPRC` corresponds to CD45RA and `ADT-PTPRC.2` corresponds to CD45RO.

## Reference annotations and signatures

Provide these files under `data/reference/`:

- `pathology_annotations.csv`: columns `sample_id`, `barcode`, and `Adenoma`, where `Adenoma` is one of `Stroma`, `Normal`, or `Low-grade`.
- `gene_sets.rds`: a named list containing at least `SenMayo`, `IEX`, and `Stem_Signature`.
- `secreted_genes.tsv`: one gene symbol per row, with a `Gene` header.
- `scrna_reference.h5ad`: raw-count curated adenoma single-cell reference for Cell2location, with `my_annotation` and `sample_id` in `obs`.
- `validation_discovery.qs` and `validation_validation.qs`: public adenoma scRNA-seq Seurat objects, if the validation analysis is rerun.
- `luminal_surface_annotations.csv`: columns `sample_id`, `barcode`, and `surface`, with manually annotated surface anchors coded `TRUE`.

## Derived files

The pipeline writes or reads the following large files under `data/derived/`:

- `visium_integrated.qs`: final Seurat object after multimodal integration and RNA-only mapping.
- `visium_scored.qs`: final object after gene-set, SenePy, cell-cycle, pseudotime, and deconvolution metadata are added.
- `epithelial_diffusion_map.qs`: epithelial subset with diffusion-map coordinates.
- `cell2location_input.h5ad`: raw-count Visium input for Cell2location, with `sample_id` in `obs`.
- `cell2location/cell2location_q05.csv`: 5th-percentile posterior cell-type abundance estimates imported from the independent Python workflow.
- `cell2location/cell2location_posterior.h5ad`: posterior AnnData including cell-type-specific expected-expression layers.
- `cell2location/NMF_colocation/`: NMF cell-type-colocalization models and exports for 6–14 factors.
- `pseudotime/*.csv`: per-section stLearn pseudotime output.
- `visium_semla.qs`: the scored Visium object with a Staffli object constructed from the Space Ranger images and coordinates; used only for semla tissue-island separation and luminal distance.

The public data record should provide these processed objects when raw-data reruns are impractical.

The released Cell2location script is a path-parameterized version of the authoritative analysis script `Data/Python/Deconv/run_cell2location.py`. The stLearn script was distilled from `Data/Python/h5ad_files/LoopforPseudotime.ipynb`. See `python/README.md`; neither Python workflow is launched from R.
