#' Plot a time-scaled phylogenetic tree
#'
#' Basic plotting helper for time-scaled phylogenies.
#'
#' @param phy A phylo object.
#' @param ... Additional arguments passed to \code{plot.phylo()}.
#'
#' @return Invisibly returns the input tree.
#' @export

plot_time_scaled <- function(phy, expand = 1.08) {
  # library(ape)
  
  # 计算树的最大时间高度
  h <- max(node.depth.edgelength(phy))
  
  # 保存原始图形参数
  op <- par(
    mar = c(5, 4, 2, 2),  # ✅ 底部和右侧加大，防止裁剪
    xpd = NA             # ✅ 允许在图框外画刻度文字
  )
  
  # 画树，右侧人为多留一点空间
  plot(
    phy,
    show.tip.label = TRUE,
    cex = 0.7,
    x.lim = c(0, h * expand),
    use.edge.length = TRUE
  )
  
  # ✅ 用稳定的 axis 代替 axisPhylo（不容易被裁剪）
  ticks <- pretty(c(0, h))
  axis(1, at = ticks, labels = round(ticks, 2))
  
  # 恢复原始参数
  par(op)
}
