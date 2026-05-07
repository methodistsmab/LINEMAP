#' Simulate lineage history and construct phylogenetic trees
#'
#' High-level wrapper for preparing transition matrices, simulating lineage
#' histories, and building phylogenetic trees.
#'
#' @param Q Transition matrix or simulation object.
#' @param gRNA_num Integer. Number of gRNAs.
#' @param mutation_rate Numeric vector. Mutation rates used for simulation.
#' @param recycle_blocks Integer. Number of recycled blocks in preprocessing.
#' @param t_max Integer. Maximum simulation time.
#' @param cell_num Integer. Initial number of cells.
#' @param which Character. guide RNA type, e.g. `"sg"` 'hg'.
#' @param use_full_edges Logical. Whether to use full edges.
#' @param keep_alive_only Logical. Whether to retain only alive cells.
#' @param set_all_edge_time_to_one Logical. Whether to set all edge lengths to 1.
#' @param rename_cells Logical. Whether to rename cell labels.
#' @param plot Logical. Whether to plot the output tree.
#' @param plot_which Character. Which tree to plot, e.g. `"phy"` tree with exact time  or `"phy1"` tree with cell division number.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list with three components:
#' \describe{
#'   \item{prep}{Output from \code{prepare_Q_lists()}}
#'   \item{sim}{Output from \code{simulate_lineage()}}
#'   \item{tree}{Output from \code{build_tree_from_sim()}}
#' }
#'
#' @export
simulate_lineage_tree <- function(Q,
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
                                  plot_which = "phy",
                                  seed = NULL) {
  if (missing(Q)) {
    stop("`Q` must be provided.")
  }
  if (!is.numeric(gRNA_num) || length(gRNA_num) != 1 || gRNA_num <= 0) {
    stop("`gRNA_num` must be a positive numeric scalar.")
  }
  if (!is.numeric(t_max) || length(t_max) != 1 || t_max <= 0) {
    stop("`t_max` must be a positive numeric scalar.")
  }
  if (!is.numeric(cell_num) || length(cell_num) != 1 || cell_num <= 0) {
    stop("`cell_num` must be a positive numeric scalar.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  prep <- prepare_Q_lists(
    Q = Q,
    gRNA.num = gRNA_num,
    mutation.rate = mutation_rate,
    recycle_blocks = recycle_blocks
  )
  
  sim <- simulate_lineage(
    prep$Q.list,
    prep$Q.list1,
    prep$N,
    t_max = t_max,
    cell.num = cell_num,
    gRNA.num = gRNA_num
  )
  
  tree_pack <- build_tree_from_sim(
    sim = sim,
    gRNA.num = gRNA_num,
    which = which,
    use_full_edges = use_full_edges,
    keep_alive_only = keep_alive_only,
    set_all_edge_time_to_one = set_all_edge_time_to_one,
    rename_cells = rename_cells,
    plot = plot,
    plot_which = plot_which
  )
  
  list(
    prep = prep,
    sim = sim,
    tree = tree_pack
  )
}