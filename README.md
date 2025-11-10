# Analysis of the features landscape from Human **Intrisically Disordered Regions (IDR)**

This repository contains the R code used to perform the bioinformatics analysis for the manuscript **"Mapping interactions between disordered regions reveals promiscuity in biomolecular condensate formation."**
It is used to compute physicochemical sequence features of human intrinsically disordered regions (IDRs), project them to a low-dimensional space (UMAP), and reproduce the analyses and plots reported in the paper.
Finally, we correlate the features used in the UMAP with the experimentally determined saturation concentration (C_sat) data for some candidates IDR presented in the paper.

The code is organized to facilitate the reproduction of figures from the manuscript and to allow for modification of the workflow to suit your own IDR data.

> **Repository:** [https://github.com/benjamin-elusers/human_idr](https://github.com/benjamin-elusers/human_idr)
> **License:** GNU General Public License v3.0 (GPL-3.0)

---

## Table of contents

* [Project structure](#project-structure)
* [What’s here](#whats-here)
* [System requirements](#system-requirements)
* [Installation](#installation)
* [Data acquisition & build](#data-acquisition--build)
* [Expected output](#expected-output)
* [How to use on your data](#how-to-use-on-your-data)
* [Reproducing manuscript results](#reproducing-manuscript-results)
* [FAQ & troubleshooting](#faq--troubleshooting)
* [License](#license)
* [How to cite](#how-to-cite)
* [Nature Research submission checklist mapping](#nature-research-submission-checklist-mapping)

---

## Project structure

```
├── build_datasets.r        # builds the processed feature table under data/
├── analysis_umap.r         # end‑to‑end UMAP analysis
├── analysis_csat.r         # C_sat analysis/plotting
├── src/                    # helper functions for I/O, features, plotting
├── data/                   # input/outputs (processed dataset lives here)
├── plots/                  # auto‑generated figures
├── LICENSE                 # GPL‑3.0
└── README.md               # this file
```

## What’s here

```
* **`analysis_umap.r`**: Main script to generate UMAP plots and feature correlations.
* **`analysis_csat.r`**: Main script to correlate IDR features with experimental C_sat values.
* **`build_datasets.r`**: Main script to download source data, calculate features, and build the analysis datasets.

* **`human_idr.Rproj`**: RStudio project file.

* **`src/`**: Contains all helper scripts and functions.
  * **`setup_main.r`**: Loads all packages, sets global options, and sources other helper scripts.
  * **`helpers_aa_features.r`**: Functions to calculate amino acid features (e.g., stickiness, charge, hydrophobicity).
  * **`helpers_idr.r`**: Functions to download and process data from MobiDB.
  * **`helpers_phasesep.r`**: Functions to download and process data from PhasePro/PhasePDB.
  * **`helpers_analysis.r`**: Helper functions for UMAP and correlation plotting.
  * **`helpers_yeastomics.r`**: General utility functions (e.g., string manipulation).

* **`data/`**: Contains input data and generated data.
 
  * **`ATAR_candidates_table.xlsx`**: The input table of IDRs and sequences used in the manuscript (the "demo" data).
  * **`*.rds *.Rdata files`**: Cached data downloaded or saved environment.
  * **`(Generated Files)`**: *.tsv and *.rdata files created by the analysis scripts.

* **`plots/`**: (Empty by default) This directory is where all output figures are saved.

> Tip: most reusable code lives under `src/`, while the top‑level scripts (`build_datasets.r`, `analysis_*.r`) run the complete workflows.

---

## System requirements

### Operating systems

* macOS 13–14 (Intel or Apple Silicon via native R)
* Ubuntu 22.04 LTS
* Windows 11 (via R 4.2+)

> The workflow is OS-agnostic and does **not** require a GPU. A laptop/desktop with ≥8 GB RAM is recommended for the full human IDR dataset.

### Software

* **R ≥ 4.2** (developed and tested with R 4.3–4.4)
* **RStudio** (optional but recommended)

### R packages

Install on first use; typical install time is 2–5 minutes on a standard desktop:

* _I/O_
  * `here`
  * `log`
  * `multidplyr`
  * `furrr
  * `rio`
  * `RJSONIO`
  * `readr`
  * `progressr`
  * `pbmcapply`
  * `xfun`

* _Data manipulation_
  * `tidyverse` (readr, dplyr, tidyr, ggplot2, purrr),
  * `magrittr`
  * `hablar`
  * `data.table`
  * `stringr`
  * `Biostrings`
  * `GenomicRanges`
  * `umap`
  * `Peptides`

* _Graphics_
  * `scales`
  * `RColorBrewer`
  * `ggcorrplot`
  * `ggalt`
  * `ggiraph`
  * `ggforce`
  * `ggrepel`
  * `ggpubr`
  * `ggeasy`
  * `patchwork`
  * `hrbrthemes`


> If you encounter a missing package at runtime, install it with `install.packages("<name>")` and re-run.

### Non‑standard hardware

None.

---

## Installation

1. **Clone the repo**

```bash
git clone https://github.com/benjamin-elusers/human_idr.git
cd human_idr
```

2. **Install R packages** (one time)

xfun allows you to load multiple packages and install the missing ones.

```r
install.packages("xfun")
pkg <- c(
  'here', 'tidyverse', 'dplyr', 'furrr', 'progressr', 'pbmcapply', 'multidplyr', 'Biostrings', 
  'GenomicRanges', 'RJSONIO', 'Peptides', 'ggcorrplot', 'umap', 'ggalt', 'ggplot2',
  'ggiraph', 'ggforce', 'ggrepel', 'see', 'ggpubr', 'hrbrthemes', 
  'ggeasy', 'patchwork', 'rio', 'magrittr', 'hablar', 'multidplyr', 'readr'
)
xfun::pkg_load2(pkg)
```

> **Typical install time:** ~5 minutes on a normal desktop with a good internet connection.

---

## Data acquisition & build

The analysis expects a single, processed **IDR feature table** (CSV/TSV) produced by `build_datasets.r`.

### Option A — Use the prebuilt dataset (recommended)

1. Download the processed dataset from the link provided by the authors: **[INSERT_DATASET_LINK_HERE]**.
2. Place the file under `data/` (e.g., `data/human_idr_features.csv`).

### Option B — Build from raw inputs

1. Open `build_datasets.r` and review the file paths/arguments at the top of the script. Many environments support `-h/--help`:

   ```bash
   Rscript build_datasets.r --help
   ```
2. Run the build:

   ```bash
   Rscript build_datasets.r
   ```

   The script downloads/assembles required inputs as needed and writes the processed feature table to `data/` (the exact path is printed at the end).

> **Expected runtime on a normal desktop:** building the processed feature table typically completes within minutes to tens of minutes depending on network speed and input sizes.

## Expected output

* UMAP scatter plot(s) written to `plots/` as PDF/PNG.
* A compact feature table
* A processed feature table (CSV) written under `data/` (e.g., `data/human_idr_features.csv`).
* UMAP scatter plot(s) written to `plots/` as PDF/PNG after running `a## How to use on your data

1. **Run/obtain the processed feature table** using **Data acquisition & build** above.

2. **Run the analyses**

* UMAP (from RStudio):

  ```r
  source("analysis_umap.r")
  ```

  or from a shell:

  ```bash
  Rscript analysis_umap.r
  ```

* C_sat (optional):

  ```r
  source("analysis_csat.r")
  ```

> UMAP hyperparameters (neighbors, `min_dist`) and input file paths are configurable at the top of `analysis_umap.r`.

> **Runtime guidance:** Full‑dataset UMAP + plotting typically completes on a normal desktop without special hardware; memory usage is modest (<1 GB) once the feature table is in memory.

---

## Reproducing manuscript results

To reproduce the UMAP panels and related analyses from the manuscript:

1. Ensure you have the full feature table used in the paper under `data/` (file name and path configurable at the top of `analysis_umap.r`).
2. Run:

   ```bash
   Rscript analysis_umap.r
   ```
3. The plots in `plots/` correspond to the UMAP figures described in the paper. Exact axis layouts and color palettes are set in the script for consistency.

For C_sat plots, run `analysis_csat.r` pointing to the appropriate input file.

> The paper’s **Methods** section details the feature definitions and UMAP rationale; this repository provides the concrete implementation for those steps.


---

## FAQ & troubleshooting

* **A package is missing.** Install it in R with `install.packages("<name>")` and re-run.
* **UMAP looks different between runs.** Set a seed at the top of `analysis_umap.r` (e.g., `set.seed(42)`) to make the layout deterministic.
* **Large tables load slowly.** Prefer CSV, avoid Excel; use `data.table::fread()`; and ensure there are no extraneous text columns.
* **Windows plotting fonts look odd.** Install additional fonts or write to PDF and view in a vector-capable viewer.

---

## License

This project is licensed under the **GNU General Public License v3.0** (GPL‑3.0). See `LICENSE` for details.

---

## How to cite

If you use this code, please cite the manuscript. (Insert DOI once available.) Example:

> Authors, *Title of the manuscript*, Journal (Year). DOI

---
