# Reproducibility testing

This guide tests the public workflow from the archived integrated Visium object through the statistical tables and manuscript figures. Run all commands from the `Visium` repository directory:

```bash
cd /path/to/Spatial-multi-omics/Visium
```

The workflow writes derived objects under `data/derived/` and tables and figures under `results/`. Preserve any existing outputs before testing if they must not be replaced.

## 1. Provide the required inputs

Place or symlink the following archived inputs at the project-relative paths below. See `data/README.md` for their complete definitions.

```text
data/derived/visium_integrated.qs
data/derived/marteau_cell2location_reference.qs
data/derived/pseudotime/*.csv
data/derived/cell2location/cell2location_q05.csv
data/derived/visium_semla.qs

data/reference/gene_sets.rds
data/reference/secreted_genes.tsv
data/reference/luminal_surface_annotations.csv
data/reference/validation_discovery.qs
data/reference/validation_validation.qs
```

The two `.qs` object paths may be symbolic links to archived `.qs`, `.qs2`, or `.rds` files. The workflow resolves each link and uses the reader appropriate for its target. Install the optional `qs2` package when a link targets a `.qs2` archive.

Install R 4.4.0 and the packages listed in `DESCRIPTION`. CytoTRACE2 should be version 1.1.0. The SenePy conda environment is defined in `environment-senepy.yml`.

## 2. Run the analysis notebook

This run starts from `visium_integrated.qs`. It does not rebuild the raw Space Ranger preprocessing, public single-cell reference, or semla object. CytoTRACE2 uses 24 cores, matching the source analysis.

```bash
Rscript -e 'rmarkdown::render(
  "Rmd/Visium_analyses.Rmd",
  params = list(
    execute = TRUE,
    run_preprocessing = FALSE,
    build_scrna_reference = FALSE,
    rebuild_staffli = FALSE,
    cores = 24
  ),
  envir = new.env()
)'
```

The analysis should generate or update:

```text
data/derived/visium_scored.qs
data/derived/epithelial_diffusion_map.qs
data/derived/single_cell_validation_full.qs
data/derived/single_cell_validation.qs
results/tables/*.csv
```

If Pandoc is unavailable, the notebook can still be executed with `knitr` and written to Markdown:

```bash
Rscript -e 'e <- new.env(parent = globalenv()); e$params <- list(execute = TRUE, run_preprocessing = FALSE, build_scrna_reference = FALSE, rebuild_staffli = FALSE, cores = 24); knitr::knit("Rmd/Visium_analyses.Rmd", output = "results/Visium_analyses.md", envir = e)'
```

## 3. Generate the figures

```bash
Rscript -e 'rmarkdown::render(
  "Rmd/Clean_plotting.Rmd",
  params = list(execute = TRUE),
  envir = new.env()
)'

Rscript -e 'rmarkdown::render(
  "Rmd/Supple_plotting.Rmd",
  params = list(execute = TRUE),
  envir = new.env()
)'
```

Main and supplementary PDFs are written under:

```text
results/figures/main/
results/figures/supplementary/
```

## 4. Verify the CytoTRACE2 subset

CytoTRACE2 should be calculated only for clusters C0, C1, C2, C3, C4, and C6. Spots in all other clusters should have missing CytoTRACE2 scores.

```bash
Rscript -e '
x <- qs::qread("data/derived/visium_scored.qs")
included <- as.character(x$combined_cluster) %in% c("0", "1", "2", "3", "4", "6")
stopifnot(all(is.na(x$CytoTRACE2_Score[!included])))
print(table(x$combined_cluster, Scored = !is.na(x$CytoTRACE2_Score)))
'
```

The command stops with an error if an excluded cluster received a score.

## 5. Compare CytoTRACE2 with the archived result

Replace the first path below with the location of the original `Feb4_cytotrace2_result.qs` object.

```r
old <- qs::qread("/path/to/Feb4_cytotrace2_result.qs")
new <- qs::qread("data/derived/visium_scored.qs")

common <- intersect(colnames(old), colnames(new))
old_score <- old@meta.data[common, "CytoTRACE2_Score"]
new_score <- new@meta.data[common, "CytoTRACE2_Score"]

cor(old_score, new_score, use = "complete.obs")
summary(new_score - old_score)
max(abs(new_score - old_score), na.rm = TRUE)
```

Verify that the barcode count and cluster composition also agree:

```r
length(common)
table(old$combined_cluster[common])
table(new$combined_cluster[common])
```

## 6. Compare statistical outputs

Compare each newly generated CSV under `results/tables/` with the corresponding archived table. For a table with a stable row order:

```r
old_table <- readr::read_csv("/path/to/original_table.csv", show_col_types = FALSE)
new_table <- readr::read_csv("results/tables/new_table.csv", show_col_types = FALSE)

all.equal(
  as.data.frame(old_table),
  as.data.frame(new_table),
  tolerance = 1e-8,
  check.attributes = FALSE
)
```

If the row order differs, sort both tables by their identifying columns before calling `all.equal()`.

Pay particular attention to:

- the three GDF15 differential-expression comparisons;
- sample-level GDF15 validation in Supplementary Figure 7A;
- unadjusted and GDF15-adjusted spatial-lag correlations;
- all five bivariate local Moran categories and both their Small-versus-Large and paired Normal-versus-Low-grade composition tests;
- the two CD8 neighborhood-density comparisons;
- cluster-proportion and continuous-score comparisons;
- adenoma tract-surface distance summaries.

## 7. Inspect the figures

Numerical agreement of the saved tables is the primary test. Inspect the regenerated PDFs for agreement in:

- samples and clusters included in each panel;
- factor and sample ordering;
- paired versus unpaired observations;
- spatial images and tissue regions;
- color scales and legends;
- qualitative trends and statistical annotations.

Stochastic or multithreaded methods can produce small numerical differences even with the same seed. Use the package versions in `DESCRIPTION`, compare values with a tolerance, and evaluate correlations and effect directions rather than requiring every generated file to be byte-for-byte identical.
