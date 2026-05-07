#' Sample matched subtrees from two phylogenies
#'
#' Randomly samples a subset of tip labels from two matched trees and returns
#' the corresponding subtrees.
#'
#' @param phy A `phylo` object, typically the time-scaled ground-truth tree.
#' @param phy1 A second `phylo` object with matching tip labels, typically the
#'   division-scaled ground-truth tree.
#' @param n_keep Integer. Number of tips to retain.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{keep_tips}{Sampled tip labels.}
#'   \item{phy_sub}{Subtree from `phy`.}
#'   \item{phy_sub1}{Subtree from `phy1`.}
#' }
#'
#' @export
sample_subtrees <- function(phy, phy1, n_keep = 30, seed = NULL) {
  if (missing(phy) || missing(phy1)) {
    stop("Both `phy` and `phy1` must be provided.")
  }
  if (!inherits(phy, "phylo")) {
    stop("`phy` must be a phylo object.")
  }
  if (!inherits(phy1, "phylo")) {
    stop("`phy1` must be a phylo object.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required to sample subtrees.")
  }
  if (!is.numeric(n_keep) || length(n_keep) != 1L || is.na(n_keep) || n_keep <= 0) {
    stop("`n_keep` must be a positive numeric scalar.")
  }
  if (!all(phy$tip.label %in% phy1$tip.label)) {
    stop("All tip labels in `phy` must also be present in `phy1`.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  tips <- phy$tip.label
  keep_tips <- sample(tips, min(as.integer(n_keep), length(tips)))

  list(
    keep_tips = keep_tips,
    phy_sub = ape::keep.tip(phy, keep_tips),
    phy_sub1 = ape::keep.tip(phy1, keep_tips)
  )
}
