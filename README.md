# 2026_RNHI

Code accompanying the RNHI study (residual normal haematopoiesis integrity score, and its
sample-/spot-level counterpart `sRNHI`): the figure code (`Figure_*`) and the R Shiny app that
scores new samples (`shinyapp/`).

The figure code is organised one directory per figure and one sub-directory per panel; each
`Panel_*/code/` folder holds the analysis and/or plotting script(s) that generate that panel.

## Layout

| Path | Contents |
| --- | --- |
| `shinyapp/` | Self-contained R Shiny application that computes RNHI for bulk expression / lineage composition, annotated scRNA-seq and annotated Visium spot data. See [`shinyapp/README.md`](shinyapp/README.md). |
| `common/figstyle.py` | Shared Matplotlib style, A4 canvas and font resolution (`RNHI_FONT`, default `Arial`, falling back to Liberation Sans / DejaVu Sans). |
| `common/panel_io.py` | Panel I/O helpers: reads statistics tables from `data/statistics/`, writes figures to `output/`. |
| `Figure_1` … `Figure_6` | Main figures, one `Panel_*/code/` folder per panel. |
| `Figure_S1` … `Figure_S4` | Supplementary figures, same layout. |
| `run_all.sh` | Runs every configured panel runner (`code/run_panel.py`) and figure builder (`00_build_*.py`). |
| `utilities/build_code_index.py` | Writes `output/code_index.csv`, mapping figure / panel / code path. |
| `utilities/crop_panels.py` | Panel cropping helper. |

## Running

```bash
bash run_all.sh              # all figures
bash run_all.sh Figure_3     # a single figure
```

A single panel can also be run from its own directory:

```bash
cd Figure_3/Panel_B_KM_RNHI_Quartile && python3 code/run_panel.py
```

Figure-level builders (`Figure_3/00_build_Figure_3.py`, `Figure_6/00_build_Figure_6.py`)
assemble the individual panels drawn by the `panel.py` modules of that figure.

### Shiny app

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "nnls", "scales"))
shiny::runApp("shinyapp")
```

The app ships with its R core (`shinyapp/rnhi_core.R`), a Python implementation of the same maths
(`shinyapp/rnhi_math.py`), simulated demo inputs and the joint deconvolution reference, so it runs
without the figure inputs described below.

## Requirements

- Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`; `pymupdf` / `pdfedit` for figure
  assembly; `anndata`, `squidpy`, `cellcharter` for the spatial panels.
- R with `ggplot2`, `data.table`, `circlize` and `qs`.

## Data

For the figure code only code is distributed: input tables are expected under `data/statistics/`
(relative to this directory, see `common/panel_io.py`) and are not included; generated figures are
written to `output/` and to the per-figure `00_output/` folders.

The Shiny app is self-contained: its simulated demo inputs and the joint reference basis for bulk
deconvolution live in `shinyapp/data/`.

## License

Released under the MIT license; see [`LICENSE`](LICENSE).
