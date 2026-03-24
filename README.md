# LINEMAP

R package for simulating lineage histories with guide-RNA barcode dynamics, building `ape` phylogenies from simulated edge tables, and plotting / subtree sampling.

## Requirements

- R >= 4.0.0
- CRAN package **ape** (>= 5.0), installed automatically with the package.

## Install from GitHub

With [remotes](https://cran.r-project.org/package=remotes):

```r
# if the Git repository root is this package folder (contains DESCRIPTION):
remotes::install_github("methodistsmab/LINEMAP")

# if the package lives in a subfolder of the repository:
remotes::install_github("methodistsmab/LINEMAP", subdir = "LINEMAP")
```

Replace `YOUR_GITHUB_USER` and `YOUR_REPO` with your account and repository name.

## Quick start

```r
library(LINEMAP)
Q <- load_example_Q()
res <- simulate_lineage_tree(Q = Q, seed = 1, plot = FALSE)
```

A longer worked example ships as `inst/scripts/run_simulation.R` (open it from the installed package path or from this source tree). The narrative readme is under `inst/doc/LINEMAP_Readme.docx`.

## Build and check locally

```sh
R CMD build LINEMAP
R CMD check --no-manual LINEMAP_0.1.0.tar.gz
```

## Author metadata

Edit `DESCRIPTION` to set `Authors@R` and the maintainer email before publishing.
