# RNHI Calculator

Shiny application that computes residual normal hematopoiesis integrity (RNHI) for AML samples from bulk expression or lineage composition data, annotated scRNA-seq data and Visium spot data.

## Definition

RNHI uses the relative composition of nine lineages of residual normal hematopoiesis: MyeProg, Mono, Gran, Mega, EarlyEry, LateEry, B, NK and DC.

1. Closure: the nine lineage abundances are rescaled to sum to 1.
2. Zero replacement: zeros are replaced with 0.65 × the detection limit of the lineage (its minimum positive proportion), and non-zero parts are rescaled (following the multiplicative replacement rule of `zCompositions::multRepl`; no lineage is removed).
3. CLR transformation: `clr(x) = log(x) - mean(log(x))`.
4. Healthy center `c`: mean CLR vector of the healthy reference donors (default: 45 BoneMarrowMap donors; zero replacement within the donors).
5. `RNHI(x) = -||clr(x) - c||_2`, the negative Aitchison distance (RNHI ≤ 0; closer to 0 = closer to healthy bone marrow).

Detection limits (both options are available in every input mode):

- Fixed (default for bulk data): limits of the bulk AML compendium (3,377 transcriptomes; `data/bulk_compendium_composition.csv`), so each sample is scored independently of the other uploaded samples. Transcriptomes without any detected lineage do not affect the limits.
- Estimated (default for scRNA-seq and Visium data): minimum positive proportion of each lineage among the uploaded query samples, so scores depend on the uploaded set.
- Healthy reference samples always use limits estimated within the reference.

Bulk expression matrices are converted to relative abundances over the reference genes, genes with zero counts are removed, and each sample is deconvolved with the MuSiC estimator (weighted non-negative least squares, default settings; `music_wnnls` in `rnhi_core.R`, which reproduces `MuSiC::music.iter`) against a joint reference of 12 normal and 7 AML cell states (`data/joint_reference_basis.rds`, 7,000 genes). At least 10% of the reference genes must be present, and at least 20% of these must be expressed in each sample. Log-transformed values are rejected.

## Input

| Input | Required columns |
|---|---|
| Composition, wide | `sample_id` and the nine lineage columns |
| Composition, long | `sample_id`, `cell_type`, `abundance` (repeated pairs are summed) |
| Bulk expression | gene symbols in the first column, one column per sample, linear-scale values |
| scRNA-seq metadata | `sample_id` or `orig.ident`, `celltype`; optional `umap_1`, `umap_2` |
| Visium | `spot_id` or `barcode`, `celltype` (metadata), or `spot_id`, `cell_type`, `abundance`; optional coordinates `x`, `y`, `imagecol`, `imagerow` or `pxl_col_in_fullres`, `pxl_row_in_fullres` |

Column names and lineage labels are matched without regard to case. Cell types other than the nine lineages must be removed before upload.

Optional `healthy_reference` column: `1` = healthy reference sample, `0` = query sample. A separate healthy reference table (same format as the query data) can also be uploaded.

## Output

`RNHI`, `aitchison_distance`, `detected_lineages`, `lineage_total` (sum of the nine abundances before closure; not used in RNHI), `qc_flag` (single detected lineage), `prop_*` (closed lineage proportions) and `clr_diff_*` (CLR difference from the healthy center; the squared differences sum to the squared Aitchison distance).

## Run locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "nnls", "scales"))
shiny::runApp("shinyapp")
```

## Files

| Path | Content |
|---|---|
| `app.R`, `rnhi_core.R` | Application and R implementation |
| `rnhi_math.py` | Python implementation |
| `data/bulk_compendium_composition.csv` | Nine-lineage compositions of the 45 healthy reference donors and of the AML transcriptomes of the bulk compendium (3,377 transcriptomes, 11 public datasets) in which at least one lineage is detected (n = 3,310) |
| `data/joint_reference_basis.rds` | Reference for bulk deconvolution |
| `data/simulated_*.csv` | Simulated example inputs (see `data/DATASETS.md`) |
| `tests/` | Regression tests (`Rscript tests/test_rnhi.R`, `python3 tests/test_rnhi_python.py`) and reference values |

## Deployment

```r
rsconnect::deployApp("shinyapp", appName = "RNHI")
```

For research use only.
