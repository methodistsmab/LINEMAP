\name{build_tree_from_sim}
\alias{build_tree_from_sim}
\title{Build phylo objects from simulation output}
\usage{
build_tree_from_sim(sim, gRNA.num, which = c("hg", "sg"),
  use_full_edges = TRUE, keep_alive_only = TRUE,
  set_all_edge_time_to_one = FALSE, rename_cells = TRUE, plot = FALSE,
  plot_which = c("phy", "phy1", "both", "none"))
}
\arguments{
  \item{sim}{Output from \code{\link{simulate_lineage}}.}
  \item{gRNA.num}{Number of gRNAs.}
  \item{which}{\code{"hg"} or \code{"sg"}.}
  \item{use_full_edges}{Use full-state edge tables.}
  \item{keep_alive_only}{Prune dead tips.}
  \item{set_all_edge_time_to_one}{Also build unit-edge tree \code{phy1}.}
  \item{rename_cells}{Rename tips to C1, C2, ...}
  \item{plot}{Plot with \code{\link{plot_time_scaled}}.}
  \item{plot_which}{Which tree to plot.}
}
\value{
List with \code{phy}, optional \code{phy1}, \code{alive.network}, \code{history}, \code{diag}.
}
\description{
Converts simulated edge tables into \code{ape::phylo} trees.
}
