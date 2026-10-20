# AML RNHI Shiny app

This app calculates the frozen AML RNHI A4 score for bulk RNA/pseudobulk, annotated scRNA-seq, and annotated Visium spot data. RNHI means **Residual Normal Hematopoiesis Index**. A4 is the nine-lineage measurable normal-hematopoiesis subcomposition used by the project: MyeProg, Mono, Gran, Mega, EarlyEry, LateEry, B, NK and DC. The final package is self-contained and includes the R calculation core in `rnhi_core.R` so it also works when the project path contains Chinese characters.

## Run locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "readr", "dplyr", "tidyr"))
# Optional but recommended for the primary bulk-expression path.
# The app has a deterministic locked-basis fallback if MuSiC is unavailable.
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install("MuSiC", ask=FALSE, update=FALSE)
shiny::runApp("/path/to/RNHI_ShinyApp")
```

## Input format

Wide composition input has one row per sample (or Visium spot) and these columns:
`sample_id` (required; `spot_id` is also accepted), `MyeProg, Mono, Gran, Mega, EarlyEry, LateEry, B, NK, DC`.

Long composition input has `sample_id, cell_type, abundance`; repeated rows are summed. For scRNA, the compact input can instead be a Seurat `meta.data` export with `celltype` and `sample_id` (or `orig.ident`): the app counts cells by sample and lineage automatically. For Visium, upload one row per annotated spot/cell with `celltype`, `spot_id` (or `sample_id`) and optional `x,y` or `imagecol,imagerow`; the app counts cells by spot and lineage automatically. The nine frozen lineage names must be used. A healthy reference can be uploaded in either composition format; otherwise the bundled 45 healthy-donor reference is used.

Bulk expression matrices are recognized automatically when the first column is a gene identifier and the remaining columns are samples. They are deconvolved with the bundled locked MuSiC reference (`data/music_basis.rds`) and then scored with the same frozen A4 RNHI core. Pre-computed nine-lineage bulk composition is also accepted. MuSiC is used when available; a deterministic locked-basis projection keeps the public example runnable on hosts where Bioconductor MuSiC is not installed, and the result message labels which path was used.
If any MuSiC sample fails quality/convergence checks, the app stops with the failed sample IDs and reasons rather than returning a partial result.

The compact metadata formats accepted by the two single-cell modes are:

| scRNA meta.data column | Example | Meaning |
|---|---|---|
| `orig.ident` (or `sample_id`) | `Van_01` | Sample used for the RNHI score |
| `celltype` | `Mono` | One of the nine frozen lineages |
| `umap_1`, `umap_2` | `1.42`, `-0.73` | Optional coordinates retained in the UMAP |

| Visium metadata column | Example | Meaning |
|---|---|---|
| `spot_id` (or `sample_id`) | `AAAC...-1` | Spot used for the RNHI score |
| `celltype` | `Gran` | Annotated lineage counted within the spot |
| `x`, `y` (or `imagecol`, `imagerow`) | `152`, `284` | Optional spatial coordinates retained in the map |

For these metadata inputs, every row is an annotated cell/spot record. The app counts rows by ID × `celltype`; you do not need to pre-compute an `abundance` column.
When image coordinates (`imagecol/imagerow` or `pxl_col_in_fullres/pxl_row_in_fullres`) are supplied, the map reverses the y-axis to preserve the usual top-left image orientation. Ordinary `x,y` coordinates are shown in their original Cartesian direction.

The same examples are shown as tables in the app guide. For copy/paste workflows, the minimum wide table is:

| sample_id | MyeProg | Mono | Gran | Mega | EarlyEry | LateEry | B | NK | DC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AML_001 | 120 | 80 | 450 | 15 | 210 | 95 | 40 | 30 | 20 |

The minimum long table is:

| sample_id | cell_type | abundance |
|---|---|---:|
| AML_001 | MyeProg | 120 |
| AML_001 | Mono | 80 |
| AML_001 | Gran | 450 |
| AML_001 | Mega | 15 |
| AML_001 | EarlyEry | 210 |
| AML_001 | LateEry | 95 |
| AML_001 | B | 40 |
| AML_001 | NK | 30 |
| AML_001 | DC | 20 |

The implementation follows the project audit: closure within the nine lineages, multiplicative zero replacement (`frac=0.65`) with healthy and query detection limits computed separately, equal-donor healthy CLR center, and `RNHI = -Aitchison distance`. Values closer to zero are more healthy-like; more negative values indicate greater compositional imbalance. `selected_lineage_total` is shown as a separate quality/context variable and does not enter RNHI. The project excludes HSC_MPP/LMPP/T from A4 because they do not meet the frozen detection criterion.

## Verified local data workflow

The real local file `composition_A4.csv` is selectable from the app with the **example** button in Bulk RNA mode. It contains 3,355 records (45 healthy reference donors and 3,310 AML query samples). The scRNA and Visium modes use small bundled examples so that their UMAP and spatial views can be tested immediately; the Visium example contains 144 spots for a readable tissue-style map. The **Refresh / reset** button returns to the input guide. The browser test uploads the same real file, clicks **Calculate RNHI**, verifies the result table, density plot and ordered composition heatmap, and writes a screenshot on demand to `tests/browser_real_data_results.png` (the screenshot is not required for deployment). The regression script compares the 3,310 AML scores with the formal `a4_rnhi_scores.csv` output.

Run the regression check with `Rscript tests/local_real_data_test.R`. The real-expression smoke test is `Rscript tests/local_music_test.R` (first eight BeatAML2 samples through MuSiC). For browser checks, start the app on port 3840 and run `Rscript tests/browser_modalities_test.R` and `Rscript tests/browser_metadata_test.R`; `Rscript tests/browser_upload_test.R` additionally uploads the full local A4 table (requires the installed `chromote` package and Chrome).

## Deployment

The folder is self-contained and can be deployed with `rsconnect::deployApp()` to shinyapps.io or Posit Connect. Do not publish identifiable clinical data.

## Bundled demo

`example_data.csv` and `example_results.csv` are included for the small demo. The full local A4 input is `composition_A4.csv`; the app's **example** button selects the mode-appropriate bundled file without a file chooser. Each input guide displays structured Markdown-style tables with the exact column names and accepted values. The Bulk guide displays both the nine-lineage composition example and the gene-by-sample MuSiC example.
