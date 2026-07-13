# Independent Python workflows

stLearn pseudotime and Cell2location deconvolution were run independently in Python/JupyterLab. They are not launched from `Rmd/My_Visium.Rmd`; the R workflow only imports their output files.

The compact scripts in this directory preserve the original analysis workflows while replacing analyst-specific paths with command-line arguments:

- `Data/Python/h5ad_files/LoopforPseudotime.ipynb` -> `run_pseudotime.py`
- `Data/Python/Deconv/run_cell2location.py` -> `run_cell2location.py`

The scripts may be run in a terminal or from a JupyterLab terminal. Paths below assume the repository root as the working directory.

## stLearn pseudotime

```bash
conda run -n adenoma-stlearn python python/run_pseudotime.py \
  --input-dir data/derived/pseudotime_input \
  --output-dir data/derived/pseudotime
```

Each input H5AD must contain raw expression, `Adenoma`, `combined_cluster`, `array_row`, `array_col`, `imagerow`, and `imagecol`, together with Visium spatial metadata. The script writes one `*.pseudotime.csv` file per section and an auditable run log.

## Cell2location

```bash
conda run -n adenoma-cell2location python python/run_cell2location.py \
  --spatial data/derived/cell2location_input.h5ad \
  --reference data/reference/scrna_reference.h5ad \
  --output-dir data/derived/cell2location \
  --reference-epochs 500 \
  --spatial-epochs 5000 \
  --train-batch-size 2500 \
  --mean-cells-per-spot 10 \
  --alpha 20
```

These defaults reproduce the authoritative script: the reference regression model is trained for 500 epochs, the spatial model for 5,000 epochs with a batch size of 2,500, and posterior export uses 1,000 samples. The spatial model uses 10 expected cells per spot and detection `alpha = 20`. Its fifth-percentile posterior abundance estimates are written to `data/derived/cell2location/cell2location_q05.csv`.

The workflow also saves the trained spatial model, runs NMF colocalization for 6 through 14 factors with three restarts, and adds cell-type-specific expected-expression matrices as layers of `cell2location_posterior.h5ad`. If `reference_signatures.csv` already exists in the output directory, it is reused, matching the original script's checkpoint behavior.

Cell2location-derived abundances are used only for the deconvolution analyses. CD8-positive spots and CD8 neighborhoods are defined exclusively from normalized ADT measurements in R.
