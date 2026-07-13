# Author decisions and remaining release items

The following analysis decisions were confirmed by the author and are implemented in the release code:

1. **RNA QC:** retain spots with total UMI count `nCount_Spatial > 300` and mitochondrial fraction `<40%`.
2. **Excluded sections:** one section failed Space Ranger QC and was never loaded; `Pilot_small_2` was loaded but excluded after failing Seurat QC. The 22-section analysis manifest therefore contains 19 RNA+ADT and 3 RNA-only sections.
3. **WNN dimensions:** use 30 RNA dimensions and 15 ADT dimensions.
4. **SenePy outliers:** use the cohort-wide mean plus two cohort-wide standard deviations.
5. **Stemness quantiles:** calculate the 25th and 75th percentiles cohort-wide across dysplastic C2/C4 spots.
6. **Senescence differential-expression regions:** define dysplastic and normal epithelium using the pathologist annotations, not cluster substitutions.
7. **GDF15 adjustment:** within each sample, residualize both senescence and stemness on spot-level GDF15, then correlate adjusted senescence with the spatial lag of adjusted neighboring stemness.
8. **CD8 source:** use normalized ADTs only. Cell2location-derived CD8 abundance is not used to define CD8-positive spots or CD8 neighborhoods.
9. **Cell2location and stLearn execution:** run these workflows independently in Python/JupyterLab. The R notebook imports their outputs and does not launch either workflow. The released Cell2location script is a path-parameterized version of the authoritative `Data/Python/Deconv/run_cell2location.py`; the stLearn script is based on `Data/Python/h5ad_files/LoopforPseudotime.ipynb`.
11. **SenePy execution:** use the lightweight R/reticulate wrapper derived from `Scripts/Utils/Cal_senepy.R`.
13. **ADT reagent mapping:** `ADT-PTPRC` corresponds to CD45RA and `ADT-PTPRC.2` corresponds to CD45RO. This mapping is recorded in `config/adt_feature_map.csv` and used by the composite CD8-positive spot gate.

Items the author will complete before making the repository public:

10. Add the final license, manuscript citation, repository DOI, raw/processed data accession, and maintainer contact.
12. Audit the compact plotting notebooks against the final figure captions and panel lettering.
