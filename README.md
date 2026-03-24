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

# Load bundled example transition matrix
Q <- load_example_Q()

# Run simulation and build phylogenetic trees
res <- simulate_lineage_tree(
  Q = Q,
  gRNA_num = 200,
  mutation_rate = seq(0.05, 0.25, 0.05),
  recycle_blocks = 40,
  t_max = 500,
  cell_num = 1,
  which = "sg",
  use_full_edges = TRUE,
  keep_alive_only = TRUE,
  set_all_edge_time_to_one = TRUE,
  rename_cells = TRUE,
  plot = FALSE,
  seed = 1
)

# Extract trees
phy <- res$tree$phy
phy1 <- res$tree$phy1

# Sample a smaller subtree for visualization
sub <- sample_subtrees(
  phy = phy,
  phy1 = phy1,
  n_keep = 30,
  seed = 1
)

# Plot the sampled time-scaled subtree
plot_time_scaled(sub$phy_sub1)
```

## Output of the sample case
![Alt text](https://github.com/methodistsmab/LINEMAP/blob/main/inst/images/output.png)


A longer worked example ships as `inst/scripts/run_simulation.R` (open it from the installed package path or from this source tree). The narrative readme is under `inst/doc/LINEMAP_Readme.docx`.

## Build and check locally

```sh
R CMD build LINEMAP
R CMD check --no-manual LINEMAP_0.1.0.tar.gz
```

## Author metadata

Edit `DESCRIPTION` to set `Authors@R` and the maintainer email before publishing.
