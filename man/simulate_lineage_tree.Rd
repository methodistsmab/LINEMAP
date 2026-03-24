\name{simulate_lineage_tree}
\alias{simulate_lineage_tree}
\title{Simulate lineage and construct phylogenetic trees}
\usage{
simulate_lineage_tree(Q, gRNA_num = 200,
  mutation_rate = seq(0.05, 0.25, 0.05), recycle_blocks = 40,
  t_max = 500, cell_num = 1, which = "sg", use_full_edges = TRUE,
  keep_alive_only = TRUE, set_all_edge_time_to_one = TRUE,
  rename_cells = TRUE, plot = FALSE, plot_which = "phy", seed = NULL)
}
\arguments{
  \item{Q}{Transition matrix.}
  \item{gRNA_num}{Number of gRNAs.}
  \item{mutation_rate}{Per-pattern mutation rates.}
  \item{recycle_blocks}{Recycling factor for gRNA list.}
  \item{t_max}{Maximum simulation time.}
  \item{cell_num}{Initial cell count.}
  \item{which}{\code{"sg"} or \code{"hg"}.}
  \item{use_full_edges, keep_alive_only, set_all_edge_time_to_one, rename_cells}{Passed to tree builder.}
  \item{plot, plot_which}{Plotting options.}
  \item{seed}{Optional RNG seed.}
}
\value{
List with \code{prep}, \code{sim}, and \code{tree}.
}
\description{
High-level wrapper: prepare Q lists, simulate, and build trees.
}
