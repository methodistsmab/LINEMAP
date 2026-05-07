#' Plot a time-scaled phylogenetic tree
#'
#' Basic plotting helper for time-scaled phylogenies.
#'
#' @param phy A phylo object.
#' @param expand Numeric scalar. Multiplier used to add right-side plotting
#'   space beyond the maximum root-to-tip height.
#' @param ... Additional arguments passed to \code{plot.phylo()}.
#'
#' @return Invisibly returns the input tree.
#' @export

plot_time_scaled <- function(phy, expand = 1.08, ...) {
  # Compute maximum root-to-tip height.
  h <- max(ape::node.depth.edgelength(phy))
  
  op <- graphics::par(
    mar = c(5, 4, 2, 2),
    xpd = NA
  )
  
  graphics::plot(
    phy,
    show.tip.label = TRUE,
    cex = 0.7,
    x.lim = c(0, h * expand),
    use.edge.length = TRUE,
    ...
  )
  
  ticks <- pretty(c(0, h))
  graphics::axis(1, at = ticks, labels = round(ticks, 2))
  
  graphics::par(op)
  invisible(phy)
}
