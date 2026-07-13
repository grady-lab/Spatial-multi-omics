# Code accompanying the spatial adenoma study

This repository contains the minimal analysis and figure-generation code used for the manuscript. Raw sequencing data, large processed objects, and patient-level source files are not stored in Git. See `data/README.md` for the required inputs and expected filenames.

## Repository layout

```text
Rmd/
  My_Visium.Rmd              Core Visium analysis, scores, and statistical tests
  Clean_plotting.Rmd         Main-figure code
  Supple_plotting.Rmd        Supplementary-figure code
R/
  preprocessing.R            Visium RNA/ADT loading, QC, integration, and mapping
  Cal_senepy.R               Lightweight SenePy wrapper used from R
  scoring.R                  irGSEA, cell-cycle, CytoTRACE, and diffusion-map helpers
  spatial_statistics.R       Moran, spatial-lag, CD8-neighborhood, distance, and composition tests
  plotting.R                 Shared palettes and plotting helpers
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
VALIDATION_NOTES.md          Items that require author confirmation before archival
```

## Reproducing the analysis

1. Install R 4.4.0 and the packages listed in `DESCRIPTION`. The manuscript versions of key packages are recorded there.
2. Create the isolated Python environments with `conda env create -f environment-senepy.yml`.
3. Put Space Ranger outputs and the externally distributed processed objects in the locations described in `data/README.md`.
4. Edit `config/sample_manifest.csv` only if the local data paths differ.
5. Run stLearn and Cell2location independently as described in `python/README.md`, or supply their archived output files.
6. Run `Rmd/My_Visium.Rmd` from the repository root. SenePy is selected through reticulate from the `adenoma-senepy` environment; stLearn and Cell2location results are imported from `data/derived/`.
7. Render `Rmd/Clean_plotting.Rmd` and `Rmd/Supple_plotting.Rmd` after the derived objects have been generated.

The notebooks use project-relative paths and never require the original analyst's home directory.

## Data availability

The accession/URL for raw and processed data must be inserted here before the repository is made public. Large objects are intentionally ignored by Git; do not force-add them.

## Citation and license

Add the final manuscript citation, repository DOI, and the authors' selected license before release.
