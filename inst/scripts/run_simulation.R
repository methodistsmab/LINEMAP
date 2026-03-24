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
