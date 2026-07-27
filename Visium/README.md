# Code accompanying the spatial adenoma study

This repository contains the analyses and figure-generation code used for the manuscript. Raw sequencing data, large processed objects, and patient-level source files are not stored in Git. See `data/README.md` for the required inputs and expected filenames.

## Repository layout

```text
Rmd/
  Visium_analyses.Rmd        Core Visium analysis, scores, and statistical tests
  Clean_plotting.Rmd         Main-figure code
  Supple_plotting.Rmd        Supplementary-figure code
R/
  preprocessing.R            Visium RNA/ADT loading, QC, integration, and mapping
  scrna_reference.R          Public scRNA-seq reference curation and H5AD export
  Cal_senepy.R               Lightweight SenePy wrapper used from R
  scoring.R                  irGSEA, cell-cycle, CytoTRACE, and diffusion-map helpers
  spatial_statistics.R       Moran, symmetric-kNN spatial-lag, CD8-neighborhood, distance, and composition tests
  plotting.R                 Shared palettes, class-faceted dot plot, and three-way spatial blend
python/
  README.md                  Standalone Python/JupyterLab execution instructions
  run_pseudotime.py          stLearn pseudotime, including disconnected tissue islands
  run_cell2location.py       Cell2location reference and spatial models
config/
  sample_manifest.csv        De-identified sample aliases and input locations
  adt_feature_map.csv        Confirmed ADT feature-to-reagent mapping
data/                        Input contract; large data are supplied separately
results/                     Generated tables and figures
METHODS_CODE_MAP.md          Mapping from the Methods section to executable code
```

## Rendered figure notebooks

The self-contained HTML notebooks for the Visium figures are available below:

- [Main figures (`Clean_plotting.html`)](https://htmlpreview.github.io/?https://github.com/grady-lab/Spatial-multi-omics/blob/main/Visium/Rmd/Clean_plotting.html)
- [Supplementary figures (`Supple_plotting.html`)](https://htmlpreview.github.io/?https://github.com/grady-lab/Spatial-multi-omics/blob/main/Visium/Rmd/Supple_plotting.html)

## Reproducing the analysis

1. Install R 4.4.0 and the packages listed in `DESCRIPTION`. The manuscript versions of key packages are recorded there.
2. Create the isolated Python environments with `conda env create -f environment-senepy.yml`.
3. Put Space Ranger outputs and the externally distributed processed objects in the locations described in `data/README.md`.
4. Edit `config/sample_manifest.csv` only if the local data paths differ.
5. If rebuilding the public single-cell reference, render `Rmd/Visium_analyses.Rmd` with `build_scrna_reference: true`; otherwise provide the archived `data/reference/scrna_reference.h5ad` and `data/derived/marteau_cell2location_reference.qs`.
6. Run stLearn and Cell2location independently as described in `python/README.md`, or supply their archived output files.
7. Run `Rmd/Visium_analyses.Rmd` from the repository root. SenePy is selected through reticulate from the `adenoma-senepy` environment; stLearn and Cell2location results are imported from `data/derived/`. To reconstruct `visium_semla.qs` from `combined_final` and the Space Ranger files, render with the parameter `rebuild_staffli: true`.
8. Render `Rmd/Clean_plotting.Rmd` and `Rmd/Supple_plotting.Rmd` after the derived objects have been generated.

The notebooks use project-relative paths.

## Data availability

Raw sequencing data and large processed objects are distributed separately from this repository. Their required layout and filenames are documented in `data/README.md`; availability information is provided with the associated manuscript.

## Citation and license

Please cite the associated manuscript when using this code.

The code is released under the GPLv3 license, as recorded in `DESCRIPTION`.
