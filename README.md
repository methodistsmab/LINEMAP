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
```

## Quick start

### Simulation data
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

# Extract ground-truth trees from the simulation
phy <- res$tree$phy    # time-scaled ground-truth tree
phy1 <- res$tree$phy1  # division-scaled ground-truth tree

# Sample a smaller matched subtree for visualization and benchmarking
sub <- sample_subtrees(
  phy = phy,
  phy1 = phy1,
  n_keep = 30,
  seed = 1
)

# Plot the sampled division-scaled ground-truth subtree
plot_time_scaled(sub$phy_sub1)

```

![Alt text](https://github.com/methodistsmab/LINEMAP/blob/main/inst/images/output.png|width=800)

### Step 0: extract barcode squence from simulated data for LINEMAP reconstruction
```
# Prepare sampled cells and barcode sequences for reconstruction.
sample <- sub$keep_tips
sequence <- sequence_from_history(res$sim$history$list_hg[sample])

```

### Step 1: build missing-aware distance dictionaries.
```
tabs_list_hgRNA <- build_pair_loglik_tables(
  P_or_list = res$prep$Q.list,
  division = 16
)

tabs_list_sgRNA <- build_pair_loglik_tables(
  P_or_list = res$prep$Q.list1,
  division = 16
)
```

### Step 2: compute the probability distance matrix.
```
distance.probability.matrix <- compute_prob_distance(
  sequence = sequence,
  tabs_list = tabs_list_hgRNA,
  gRNA.num = 200,
  barcode.label = 1,
  add_C0 = TRUE,
  missing_label = "MISSING"
)
```

### Step 3: build a predicted LINEMAP tree from the distance matrix.
```
tree.LINEMAP <- build_tree_from_distance(
  distance_matrix = distance.probability.matrix,
  method = "NJ",
  root_outgroup = "C0",
  drop_outgroup = TRUE
)
```
### Step 4: compare the predicted tree with the ground-truth time-scaled tree.
```
topology.comparison <- compare_tree_topology(
  tree1 = res$tree$phy,
  tree2 = tree.LINEMAP,
  plot_tangle = TRUE,
  main_left = "Ground truth",
  main_right = "Predicted"
)
```
###  Output of topology.comparison$similarity_scores
```
> topology.comparison$similarity_scores
            metric     value
1   cophenetic_cor 0.9832343
2     bakers_gamma 0.9890536
3 1 - entanglement 1.0000000
4    RF_similarity 1.0000000
5   Nye_similarity 1.0000000
6   JRF_similarity 1.0000000
```
### Output of topology.comparison$distance_scores
```
> topology.comparison$distance_scores
                 metric      value
1    1 - cophenetic_cor 0.01676568
2      1 - bakers_gamma 0.01094640
3          entanglement 0.00000000
4           RF_distance 0.00000000
5          Nye_distance 0.00000000
6          JRF_distance 0.00000000
7 TreeDistance_distance 0.00000000
8  TripletDistance_norm 0.00000000
```
![Alt text](https://github.com/methodistsmab/LINEMAP/blob/main/inst/images/tree_comparison.jpg)

### Step 5: compare pairwise MRCA heights and cell-depth relationships.
```
time.lineage.comparison <- compare_time_lineage(
  true_tree_time = res$tree$phy,
  true_tree_div = res$tree$phy1,
  pred_tree = tree.LINEMAP,
  dist_prob = distance.probability.matrix,
  true_history = res$tree$history,
  cor_method = "spearman"
)
```
### Output of time.lineage.comparison
```
> time.lineage.comparison
$cells_used
 [1] "C37"  "C105" "C129" "C187" "C270" "C277" "C299" "C307" "C330" "C382"
[11] "C471" "C485" "C494" "C509" "C591" "C597" "C601" "C677" "C679" "C725"
[21] "C729" "C775" "C801" "C802" "C836" "C841" "C852" "C874" "C878" "C907"

$pred_tree_type
[1] "phylo"

$correlation_method
[1] "spearman"

$cor_probDist_vs_trueMRCA_time
[1] 0.9190899

$cor_probDist_vs_trueMRCA_div
[1] 0.9236095

$cor_predMRCA_vs_trueMRCA_time
[1] 0.9990148

$cor_predMRCA_vs_trueMRCA_div
[1] 0.9964034

$cor_treeHeight_vs_trueDivisions
[1] 0.3006322

$cor_C0ProbDist_vs_trueDivisions
[1] 0.4524147

$pred_root2tip_depth
     C37     C105     C129     C187     C270     C277     C299     C307 
661.3133 729.8261 762.2762 718.4168 700.8992 732.8036 714.9486 740.7326 
    C330     C382     C471     C485     C494     C509     C591     C597 
749.0337 744.5566 688.8827 756.2378 712.7056 696.3559 736.6800 718.3704 
    C601     C677     C679     C725     C729     C775     C801     C802 
802.6411 711.5546 792.3080 728.6928 780.1602 740.2327 699.7758 726.0321 
    C836     C841     C852     C874     C878     C907 
794.4798 737.2712 785.9338 695.7874 748.5610 725.5439 

$C0_to_cell_prob_dist
     C37     C105     C129     C187     C270     C277     C299     C307 
1515.221 1623.757 1663.267 1554.415 1576.904 1549.802 1582.838 1622.366 
    C330     C382     C471     C485     C494     C509     C591     C597 
1629.778 1578.719 1599.701 1687.504 1594.556 1542.821 1571.431 1587.826 
    C601     C677     C679     C725     C729     C775     C801     C802 
1644.106 1599.293 1697.437 1599.400 1605.467 1602.284 1571.178 1612.546 
    C836     C841     C852     C874     C878     C907 
1660.512 1610.366 1730.166 1569.442 1579.336 1598.827 
```
# Wrapper code:
### Steps 1-3 can also be run in one call with the reconstruction wrapper:
```
tree.wrapper <- build_tree_from_sequence(
  sequence = sequence,
  P_or_list = res$prep$Q.list,
  division = 16,
  gRNA.num = 200,
  method = "NJ",
  add_C0 = TRUE,
  root_outgroup = "C0",
  drop_outgroup = TRUE,
  missing_label = "MISSING"
)

tree.wrapper$tree
```
## Steps 4-5 can also be returned together with the comparison wrapper:
```
combined.comparison <- compare_reconstructed_tree(
  pred_tree = tree.LINEMAP,
  true_tree = res$tree$phy,
  true_tree_time = res$tree$phy,
  true_tree_div = res$tree$phy1,
  dist_prob = distance.probability.matrix,
  true_history = res$tree$history,
  plot_tangle = FALSE,
  cor_method = "spearman"
)

combined.comparison
```

A longer worked example ships as `inst/scripts/run_simulation.R` 

## Build and check locally

```sh
R CMD build LINEMAP
R CMD check --no-manual LINEMAP_0.1.0.tar.gz
```

## Author metadata

Edit `DESCRIPTION` to set `Authors@R` and the maintainer email before publishing.
