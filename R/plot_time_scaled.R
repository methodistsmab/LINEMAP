#' Plot a time-scaled phylogenetic tree
#'
#' Basic plotting helper for time-scaled phylogenies.
#'
#' @param phy A phylo object.
#' @param expand Horizontal padding factor beyond maximum tip depth.
#' @param ... Reserved for future use (passed to plotting; currently unused).
#'
#' @return Invisibly returns `phy`.
#' @export
plot_time_scaled <- function(phy, expand = 1.08, ...) {
  h <- max(node.depth.edgelength(phy))

  op <- par(
    mar = c(5, 4, 2, 2),
    xpd = NA
  )

  plot(
    phy,
    show.tip.label = TRUE,
    cex = 0.7,
    x.lim = c(0, h * expand),
    use.edge.length = TRUE
  )

  ticks <- pretty(c(0, h))
  axis(1, at = ticks, labels = round(ticks, 2))

  par(op)
  invisible(phy)
}
