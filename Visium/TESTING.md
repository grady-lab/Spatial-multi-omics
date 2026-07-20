# Validation performed

The release draft was checked on 2026-07-09 as follows:

- All R source files parse under R 4.4.0.
- All R chunks extracted from the three R Markdown notebooks parse under R 4.4.0.
- Both Python scripts pass Python 3 syntax compilation.
- Cell2location stages and defaults were checked against the authoritative `Data/Python/Deconv/run_cell2location.py`, including reference/spatial training, posterior export, NMF colocalization, and cell-type-specific expected-expression layers.
- Synthetic-data tests pass for spatial lag correlations, including an explicit hand-calculated check that both senescence and stemness are residualized on GDF15 before adjusted stemness is spatially lagged. Tests also pass for the before/after GDF15 comparison, tissue connected components, CD8 neighborhood density, and luminal-surface distance.
- No analyst-specific absolute paths occur in the release directory.
- The sample manifest resolves to 22 included sections (19 RNA+ADT and 3 RNA-only).
- No source file in the release is larger than 1 MB; raw and derived scientific data types are ignored by Git.

Full numerical reproduction was not run because the raw data, reference files, and large derived objects are external inputs. HTML rendering was not tested in this compute session because Pandoc is unavailable; notebook code extraction and R parsing succeeded.

On 2026-07-20, synthetic-data parity tests confirmed that the public `compare_score_clusters_paired()` and `run_comp_tests()` implementations match the authoritative helper script in both paired spot-level and unpaired sample-level modes. The categorical-proportion parity test used the locally available vegan 2.6.4; the release environment remains pinned to vegan 2.7.2 in `DESCRIPTION`.
