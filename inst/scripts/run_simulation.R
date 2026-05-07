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

# -------------------------------------------------------------------------
# Reviewer example: reconstruct a LINEMAP tree from simulated barcode data
# -------------------------------------------------------------------------

# Prepare sampled cells and barcode sequences for reconstruction.
sample <- sub$keep_tips
sequence <- sequence_from_history(res$sim$history$list_hg[sample])

# Step 1: build missing-aware distance dictionaries.
tabs_list_hgRNA <- build_pair_loglik_tables(
  P_or_list = res$prep$Q.list,
  division = 16
)

tabs_list_sgRNA <- build_pair_loglik_tables(
  P_or_list = res$prep$Q.list1,
  division = 16
)

# Step 2: compute the probability distance matrix.
distance.probability.matrix <- compute_prob_distance(
  sequence = sequence,
  tabs_list = tabs_list_hgRNA,
  gRNA.num = 200,
  barcode.label = 1,
  add_C0 = TRUE,
  missing_label = "MISSING"
)

# Step 3: build a predicted LINEMAP tree from the distance matrix.
tree.LINEMAP <- build_tree_from_distance(
  distance_matrix = distance.probability.matrix,
  method = "NJ",
  root_outgroup = "C0",
  drop_outgroup = TRUE
)

# Steps 1-3 can also be run in one call with the reconstruction wrapper:
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

# Step 4: compare the predicted tree with the ground-truth time-scaled tree.
topology.comparison <- compare_tree_topology(
  tree1 = res$tree$phy,
  tree2 = tree.LINEMAP,
  plot_tangle = TRUE,
  main_left = "Ground truth",
  main_right = "Predicted"
)

topology.comparison$similarity_scores
topology.comparison$distance_scores

# Step 5: compare pairwise MRCA heights and cell-depth relationships.
time.lineage.comparison <- compare_time_lineage(
  true_tree_time = res$tree$phy,
  true_tree_div = res$tree$phy1,
  pred_tree = tree.LINEMAP,
  dist_prob = distance.probability.matrix,
  true_history = res$tree$history,
  cor_method = "spearman"
)

time.lineage.comparison

# Steps 4-5 can also be returned together with the comparison wrapper:
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
