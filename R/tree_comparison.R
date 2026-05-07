#' Extract a tree object from a common result container
#'
#' Internal helper that lets comparison functions accept either a tree directly
#' or a common wrapper result such as `list(tree = ...)` or `list(phy = ...)`.
#'
#' @param x Tree-like object or result container.
#'
#' @return A tree-like object.
.extract_tree_object <- function(x) {
  if (inherits(x, c("phylo", "hclust", "dendrogram", "dist")) || is.matrix(x)) {
    return(x)
  }

  if (is.list(x)) {
    if (!is.null(x$tree)) {
      return(.extract_tree_object(x$tree))
    }
    if (!is.null(x$phy)) {
      return(.extract_tree_object(x$phy))
    }
    if (!is.null(x$phy1)) {
      return(.extract_tree_object(x$phy1))
    }
  }

  x
}

#' Coerce an object to a dendrogram
#'
#' Internal helper for topology comparison. Supported inputs are `dendrogram`,
#' `hclust`, `phylo`, `dist`, square distance matrices, and common result
#' containers such as `list(tree = ...)` or `list(phy = ...)`.
#'
#' @param x Tree-like object.
#'
#' @return A `dendrogram` object.
.coerce_to_dend <- function(x) {
  x <- .extract_tree_object(x)

  if (inherits(x, "dendrogram")) {
    return(x)
  }
  if (inherits(x, "hclust")) {
    return(stats::as.dendrogram(x))
  }
  if (inherits(x, "phylo")) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package `ape` is required to coerce `phylo` objects.")
    }
    if (!ape::is.binary(x)) {
      x <- ape::multi2di(x)
    }
    d <- ape::cophenetic.phylo(x)
    return(stats::as.dendrogram(stats::hclust(stats::as.dist(d), method = "average")))
  }
  if (inherits(x, "dist")) {
    return(stats::as.dendrogram(stats::hclust(x, method = "average")))
  }
  if (is.matrix(x)) {
    return(stats::as.dendrogram(stats::hclust(stats::as.dist(x), method = "average")))
  }

  stop("Unsupported tree object class: ", paste(class(x), collapse = "/"))
}

#' Coerce an object to phylo
#'
#' Internal helper for topology comparison. Supported inputs are `phylo`,
#' `dendrogram`, `hclust`, `dist`, square distance matrices, and common result
#' containers such as `list(tree = ...)` or `list(phy = ...)`.
#'
#' @param x Tree-like object.
#'
#' @return A `phylo` object.
.coerce_to_phylo <- function(x) {
  x <- .extract_tree_object(x)

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for phylo conversion.")
  }

  if (inherits(x, "phylo")) {
    if (!ape::is.binary(x)) {
      x <- ape::multi2di(x)
    }
    return(x)
  }
  if (inherits(x, "dendrogram")) {
    return(ape::as.phylo(x))
  }
  if (inherits(x, "hclust")) {
    return(ape::as.phylo(x))
  }
  if (inherits(x, "dist")) {
    return(ape::as.phylo(stats::hclust(x, method = "average")))
  }
  if (is.matrix(x)) {
    return(ape::as.phylo(stats::hclust(stats::as.dist(x), method = "average")))
  }

  stop("Unsupported tree object class: ", paste(class(x), collapse = "/"))
}

#' Align two dendrograms by shared labels
#'
#' Internal helper for topology comparison and tanglegram plotting.
#'
#' @param d1 First dendrogram.
#' @param d2 Second dendrogram.
#'
#' @return A `dendlist` containing aligned dendrograms.
.align_and_prune_dendrograms <- function(d1, d2) {
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Package `dendextend` is required for dendrogram alignment.")
  }

  labs1 <- labels(d1)
  labs2 <- labels(d2)
  if (is.null(labs1) || is.null(labs2)) {
    stop("Both trees must have leaf labels.")
  }

  labs <- intersect(labs1, labs2)
  if (length(labs) < 2L) {
    stop("The two trees share fewer than two leaf labels.")
  }

  d1p <- dendextend::prune(d1, setdiff(labs1, labs))
  d2p <- dendextend::prune(d2, setdiff(labs2, labs))

  target_order <- labels(d1p)
  d2p <- dendextend::rotate(d2p, order = target_order)

  dl <- dendextend::dendlist(d1p, d2p)
  tryCatch(
    dendextend::untangle(dl, method = "step1side"),
    error = function(e) dl
  )
}

#' Align two phylo trees by shared labels
#'
#' Internal helper for topology comparison.
#'
#' @param tree1 First `phylo` object.
#' @param tree2 Second `phylo` object.
#' @param labels_keep Optional labels to retain. If `NULL`, uses label
#'   intersection.
#'
#' @return A list with aligned `tree1`, `tree2`, and `labels`.
.align_and_prune_phylo <- function(tree1, tree2, labels_keep = NULL) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for phylo alignment.")
  }

  if (is.null(labels_keep)) {
    labels_keep <- intersect(tree1$tip.label, tree2$tip.label)
  }
  if (length(labels_keep) < 2L) {
    stop("The two trees share fewer than two tip labels.")
  }

  list(
    tree1 = ape::drop.tip(tree1, setdiff(tree1$tip.label, labels_keep)),
    tree2 = ape::drop.tip(tree2, setdiff(tree2$tip.label, labels_keep)),
    labels = labels_keep
  )
}

#' Compare two reconstructed trees by topology
#'
#' Computes common tree comparison metrics after pruning both trees to their
#' shared labels. This function accepts `phylo`, `hclust`, `dendrogram`, `dist`,
#' or square distance matrix inputs.
#'
#' @param tree1 First tree-like object.
#' @param tree2 Second tree-like object.
#' @param plot_tangle Logical; if `TRUE`, draw a tanglegram.
#' @param main_left Left-side title for the tanglegram.
#' @param main_right Right-side title for the tanglegram.
#'
#' @return A list containing similarity scores, distance scores, aligned
#'   dendrograms, and aligned `phylo` trees.
#'
#' @export
compare_tree_topology <- function(tree1,
                                  tree2,
                                  plot_tangle = FALSE,
                                  main_left = "Ground truth",
                                  main_right = "Predicted") {
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Package `dendextend` is required for topology comparison.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for topology comparison.")
  }

  d1 <- .coerce_to_dend(tree1)
  d2 <- .coerce_to_dend(tree2)
  p1 <- .coerce_to_phylo(tree1)
  p2 <- .coerce_to_phylo(tree2)

  dl <- .align_and_prune_dendrograms(d1, d2)
  d1a <- dl[[1]]
  d2a <- dl[[2]]
  labs <- labels(d1a)

  phy <- .align_and_prune_phylo(p1, p2, labels_keep = labs)
  p1a <- phy$tree1
  p2a <- phy$tree2

  cophenetic_cor <- dendextend::cor_cophenetic(d1a, d2a)
  bakers_gamma <- dendextend::cor_bakers_gamma(d1a, d2a)
  entangle <- dendextend::entanglement(dl)

  rf_similarity <- NA_real_
  rf_distance <- NA_real_
  nye_similarity <- NA_real_
  nye_distance <- NA_real_
  jrf_similarity <- NA_real_
  jrf_distance <- NA_real_
  tree_distance <- NA_real_

  if (requireNamespace("TreeDist", quietly = TRUE)) {
    rf_similarity <- TreeDist::RobinsonFoulds(p1a, p2a, similarity = TRUE, normalize = TRUE)
    rf_distance <- TreeDist::RobinsonFoulds(p1a, p2a, similarity = FALSE, normalize = TRUE)
    nye_similarity <- TreeDist::NyeSimilarity(p1a, p2a, similarity = TRUE, normalize = TRUE)
    nye_distance <- TreeDist::NyeSimilarity(p1a, p2a, similarity = FALSE, normalize = TRUE)
    jrf_similarity <- TreeDist::JaccardRobinsonFoulds(p1a, p2a, similarity = TRUE, normalize = TRUE)
    jrf_distance <- TreeDist::JaccardRobinsonFoulds(p1a, p2a, similarity = FALSE, normalize = TRUE)
    tree_distance <- TreeDist::TreeDistance(p1a, p2a)
  }

  triplet_distance_norm <- NA_real_
  if (requireNamespace("Quartet", quietly = TRUE) && length(p1a$tip.label) >= 3L) {
    f1 <- tempfile(fileext = ".nwk")
    f2 <- tempfile(fileext = ".nwk")
    on.exit(unlink(c(f1, f2)), add = TRUE)
    ape::write.tree(p1a, f1)
    ape::write.tree(p2a, f2)
    triplet_raw <- Quartet::TripletDistance(f1, f2)
    triplet_distance_norm <- triplet_raw / choose(length(p1a$tip.label), 3)
  }

  if (isTRUE(plot_tangle)) {
    dendextend::tanglegram(
      d1a, d2a,
      main_left = main_left,
      main_right = main_right,
      lab.cex = 0.8,
      edge.lwd = 1.2,
      highlight_distinct_edges = TRUE,
      common_subtrees_color_lines = TRUE
    )
  }

  similarity_scores <- data.frame(
    metric = c(
      "cophenetic_cor",
      "bakers_gamma",
      "1 - entanglement",
      "RF_similarity",
      "Nye_similarity",
      "JRF_similarity"
    ),
    value = c(
      cophenetic_cor,
      bakers_gamma,
      1 - entangle,
      rf_similarity,
      nye_similarity,
      jrf_similarity
    ),
    row.names = NULL
  )

  distance_scores <- data.frame(
    metric = c(
      "1 - cophenetic_cor",
      "1 - bakers_gamma",
      "entanglement",
      "RF_distance",
      "Nye_distance",
      "JRF_distance",
      "TreeDistance_distance",
      "TripletDistance_norm"
    ),
    value = c(
      1 - cophenetic_cor,
      1 - bakers_gamma,
      entangle,
      rf_distance,
      nye_distance,
      jrf_distance,
      tree_distance,
      triplet_distance_norm
    ),
    row.names = NULL
  )

  list(
    labels_used = labs,
    similarity_scores = similarity_scores,
    distance_scores = distance_scores,
    aligned_dendlist = dl,
    aligned_phylo = list(tree1 = p1a, tree2 = p2a)
  )
}

#' Compute MRCA height for one pair of tips
#'
#' @param tree Rooted `phylo` tree.
#' @param tip1 First tip label.
#' @param tip2 Second tip label.
#' @param from_present Logical; if `TRUE`, return distance from present back to
#'   MRCA. If `FALSE`, return root-to-MRCA depth.
#' @param T_max Optional total tree height. If `NULL`, uses maximum tip depth.
#'
#' @return Numeric MRCA height.
#'
#' @export
mrca_height_pair <- function(tree, tip1, tip2, from_present = TRUE, T_max = NULL) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for MRCA height calculation.")
  }
  stopifnot(inherits(tree, "phylo"))
  if (!ape::is.rooted(tree)) {
    stop("`tree` must be rooted.")
  }

  tip1 <- as.character(tip1)
  tip2 <- as.character(tip2)
  if (!tip1 %in% tree$tip.label) {
    stop("`tip1` not found in tree.")
  }
  if (!tip2 %in% tree$tip.label) {
    stop("`tip2` not found in tree.")
  }

  depths <- ape::node.depth.edgelength(tree)
  ntip <- length(tree$tip.label)
  total_height <- if (is.null(T_max)) max(depths[seq_len(ntip)]) else T_max
  mrca_node <- ape::getMRCA(tree, c(tip1, tip2))
  if (is.null(mrca_node) || is.na(mrca_node)) {
    stop("MRCA not found.")
  }

  depth_from_root <- depths[mrca_node]
  if (isTRUE(from_present)) {
    total_height - depth_from_root
  } else {
    depth_from_root
  }
}

#' Compute an MRCA height matrix
#'
#' @param tree Rooted `phylo` tree.
#' @param from_present Logical; if `TRUE`, return distance from present back to
#'   each MRCA. If `FALSE`, return root-to-MRCA depth.
#' @param T_max Optional total tree height. If `NULL`, uses maximum tip depth.
#'
#' @return Square matrix indexed by tip labels.
#'
#' @export
mrca_height_matrix <- function(tree, from_present = TRUE, T_max = NULL) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for MRCA height calculation.")
  }
  stopifnot(inherits(tree, "phylo"))
  if (!ape::is.rooted(tree)) {
    stop("`tree` must be rooted.")
  }

  tips <- tree$tip.label
  ntip <- length(tips)
  depths <- ape::node.depth.edgelength(tree)
  tip_depths <- depths[seq_len(ntip)]
  total_height <- if (is.null(T_max)) max(tip_depths) else T_max

  mrca_node_mat <- ape::mrca(tree)
  depth_mrca <- matrix(
    depths[mrca_node_mat],
    nrow = ntip,
    ncol = ntip,
    dimnames = list(tips, tips)
  )

  if (isTRUE(from_present)) {
    out <- total_height - depth_mrca
    diag(out) <- 0
  } else {
    out <- depth_mrca
    diag(out) <- tip_depths
  }

  out
}

#' Safely compute correlation
#'
#' Internal helper that ignores non-finite pairs and returns `NA` when fewer than
#' three complete observations are available.
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param method Correlation method.
#'
#' @return Numeric scalar.
.cor_safe <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

#' Extract first-merge height for hclust tips
#'
#' Internal helper used when a predicted tree is an `hclust` object.
#'
#' @param hc An `hclust` object.
#'
#' @return Named numeric vector indexed by leaf labels.
.leaf_first_merge_height <- function(hc) {
  n <- length(hc$labels)
  out <- numeric(n)
  names(out) <- hc$labels

  for (i in seq_len(n)) {
    row_idx <- which(hc$merge[, 1] == -i | hc$merge[, 2] == -i)[1]
    out[i] <- hc$height[row_idx]
  }

  out
}

#' Compare reconstructed tree against true time/division trees
#'
#' Computes correlations between predicted tree structure, probability distance,
#' and known true lineage quantities. This is useful for simulation benchmarks
#' where both a true time-scaled tree and a true division-scaled tree are known.
#'
#' @param true_tree_time True time-scaled rooted `phylo` tree.
#' @param true_tree_div True division-scaled rooted `phylo` tree.
#' @param pred_tree Reconstructed tree, either `phylo` or `hclust`.
#' @param dist_prob Optional probability distance matrix, usually the output of
#'   `compute_prob_distance()`.
#' @param true_divisions Optional named numeric vector of true division counts per
#'   cell. If `NULL`, it is computed from `true_history`.
#' @param true_history Optional list of terminal cell history matrices, such as
#'   `res$tree$history` from `simulate_lineage_tree()`. Used to compute
#'   `true_divisions = sapply(true_history, ncol) - 1`.
#' @param c0_name Name of artificial root/outgroup in `dist_prob`.
#' @param from_present Logical passed to `mrca_height_matrix()`.
#' @param cor_method Correlation method, either `"pearson"` or `"spearman"`.
#'
#' @return A list of correlations and diagnostic vectors.
#'
#' @export
compare_time_lineage <- function(true_tree_time,
                                 true_tree_div,
                                 pred_tree,
                                 dist_prob = NULL,
                                 true_divisions = NULL,
                                 true_history = NULL,
                                 c0_name = "C0",
                                 from_present = TRUE,
                                 cor_method = c("pearson", "spearman")) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package `ape` is required for time-lineage comparison.")
  }

  cor_method <- match.arg(cor_method)
  stopifnot(
    inherits(true_tree_time, "phylo"),
    inherits(true_tree_div, "phylo"),
    inherits(pred_tree, "phylo") || inherits(pred_tree, "hclust")
  )

  if (is.null(true_divisions)) {
    if (is.null(true_history)) {
      stop("Provide either `true_divisions` or `true_history`.")
    }
    if (!is.list(true_history) || length(true_history) == 0L) {
      stop("`true_history` must be a non-empty list of cell history matrices.")
    }
    true_divisions <- stats::setNames(
      vapply(true_history, ncol, integer(1L)) - 1L,
      names(true_history)
    )
  }

  if (!is.numeric(true_divisions) || is.null(names(true_divisions))) {
    stop("`true_divisions` must be a named numeric vector.")
  }

  has_dist <- !is.null(dist_prob)
  if (has_dist) {
    dist_prob <- as.matrix(dist_prob)
  }

  pred_type <- if (inherits(pred_tree, "phylo")) "phylo" else "hclust"
  pred_phy <- if (pred_type == "phylo") pred_tree else ape::as.phylo(pred_tree)

  parts <- list(
    true_tree_time$tip.label,
    true_tree_div$tip.label,
    pred_phy$tip.label,
    names(true_divisions)
  )

  if (has_dist) {
    dist_cells <- intersect(rownames(dist_prob), colnames(dist_prob))
    dist_cells <- setdiff(dist_cells, c0_name)
    parts <- c(parts, list(dist_cells))
  }

  cells <- Reduce(intersect, parts)
  if (length(cells) < 3L) {
    stop("Fewer than three shared cells are available for comparison.")
  }

  time_sub <- ape::drop.tip(true_tree_time, setdiff(true_tree_time$tip.label, cells))
  div_sub <- ape::drop.tip(true_tree_div, setdiff(true_tree_div$tip.label, cells))
  pred_sub <- ape::drop.tip(pred_phy, setdiff(pred_phy$tip.label, cells))
  true_div_sub <- true_divisions[cells]

  H_true_time <- mrca_height_matrix(time_sub, from_present = from_present)[cells, cells]
  H_true_div <- mrca_height_matrix(div_sub, from_present = from_present)[cells, cells]
  idx_lower <- lower.tri(H_true_time)

  v_true_time <- H_true_time[idx_lower]
  v_true_div <- H_true_div[idx_lower]

  cor_prob_vs_trueMRCA_time <- NA_real_
  cor_prob_vs_trueMRCA_div <- NA_real_
  if (has_dist) {
    v_prob <- dist_prob[cells, cells][idx_lower]
    cor_prob_vs_trueMRCA_time <- .cor_safe(v_prob, v_true_time, cor_method)
    cor_prob_vs_trueMRCA_div <- .cor_safe(v_prob, v_true_div, cor_method)
  }

  cor_predMRCA_vs_trueMRCA_time <- NA_real_
  cor_predMRCA_vs_trueMRCA_div <- NA_real_
  if (ape::is.rooted(pred_sub)) {
    H_pred <- mrca_height_matrix(pred_sub, from_present = from_present)[cells, cells]
    v_pred <- H_pred[idx_lower]
    cor_predMRCA_vs_trueMRCA_time <- .cor_safe(v_pred, v_true_time, cor_method)
    cor_predMRCA_vs_trueMRCA_div <- .cor_safe(v_pred, v_true_div, cor_method)
  }

  if (pred_type == "phylo") {
    pred_depth_all <- ape::node.depth.edgelength(pred_sub)
    pred_root2tip <- pred_depth_all[match(cells, pred_sub$tip.label)]
    names(pred_root2tip) <- cells
  } else {
    h_leaf <- .leaf_first_merge_height(pred_tree)
    pred_root2tip <- as.numeric(h_leaf[cells])
    names(pred_root2tip) <- cells
  }

  cor_treeHeight_vs_trueDivisions <- .cor_safe(pred_root2tip, true_div_sub, cor_method)

  dist_C0 <- NULL
  cor_C0ProbDist_vs_trueDivisions <- NA_real_
  if (has_dist && c0_name %in% rownames(dist_prob) && all(cells %in% colnames(dist_prob))) {
    dist_C0 <- as.numeric(dist_prob[c0_name, cells])
    names(dist_C0) <- cells
    cor_C0ProbDist_vs_trueDivisions <- .cor_safe(dist_C0, true_div_sub, cor_method)
  } else if (has_dist && c0_name %in% colnames(dist_prob) && all(cells %in% rownames(dist_prob))) {
    dist_C0 <- as.numeric(dist_prob[cells, c0_name])
    names(dist_C0) <- cells
    cor_C0ProbDist_vs_trueDivisions <- .cor_safe(dist_C0, true_div_sub, cor_method)
  }

  list(
    cells_used = cells,
    pred_tree_type = pred_type,
    correlation_method = cor_method,
    cor_probDist_vs_trueMRCA_time = cor_prob_vs_trueMRCA_time,
    cor_probDist_vs_trueMRCA_div = cor_prob_vs_trueMRCA_div,
    cor_predMRCA_vs_trueMRCA_time = cor_predMRCA_vs_trueMRCA_time,
    cor_predMRCA_vs_trueMRCA_div = cor_predMRCA_vs_trueMRCA_div,
    cor_treeHeight_vs_trueDivisions = cor_treeHeight_vs_trueDivisions,
    cor_C0ProbDist_vs_trueDivisions = cor_C0ProbDist_vs_trueDivisions,
    pred_root2tip_depth = pred_root2tip,
    C0_to_cell_prob_dist = dist_C0
  )
}

#' Compare a reconstructed tree with optional true references
#'
#' User-facing convenience wrapper that can run topology comparison, simulation
#' time/division comparison, or both depending on which reference objects are
#' provided.
#'
#' @param pred_tree Reconstructed tree.
#' @param true_tree Optional true tree for topology comparison.
#' @param true_tree_time Optional true time-scaled tree.
#' @param true_tree_div Optional true division-scaled tree.
#' @param dist_prob Optional probability distance matrix.
#' @param true_divisions Optional named numeric vector of true division counts.
#'   If `NULL`, it can be computed from `true_history`.
#' @param true_history Optional list of terminal cell history matrices, such as
#'   `res$tree$history` from `simulate_lineage_tree()`.
#' @param plot_tangle Logical; if `TRUE`, draw topology tanglegram.
#' @param main_left Left-side title for the tanglegram.
#' @param main_right Right-side title for the tanglegram.
#' @param c0_name Name of artificial root/outgroup in `dist_prob`.
#' @param from_present Logical passed to MRCA comparison.
#' @param cor_method Correlation method.
#'
#' @return A list with `topology` and/or `time_lineage` components.
#'
#' @export
compare_reconstructed_tree <- function(pred_tree,
                                       true_tree = NULL,
                                       true_tree_time = NULL,
                                       true_tree_div = NULL,
                                       dist_prob = NULL,
                                       true_divisions = NULL,
                                       true_history = NULL,
                                       plot_tangle = FALSE,
                                       main_left = "Ground truth",
                                       main_right = "Predicted",
                                       c0_name = "C0",
                                       from_present = TRUE,
                                       cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)
  out <- list()

  if (!is.null(true_tree)) {
    out$topology <- compare_tree_topology(
      tree1 = true_tree,
      tree2 = pred_tree,
      plot_tangle = plot_tangle,
      main_left = main_left,
      main_right = main_right
    )
  }

  has_time_inputs <- !is.null(true_tree_time) ||
    !is.null(true_tree_div) ||
    !is.null(true_divisions) ||
    !is.null(true_history)
  if (has_time_inputs) {
    if (is.null(true_tree_time) || is.null(true_tree_div)) {
      stop("`true_tree_time` and `true_tree_div` are required for time-lineage comparison.")
    }
    if (is.null(true_divisions) && is.null(true_history)) {
      stop("Provide either `true_divisions` or `true_history` for time-lineage comparison.")
    }

    out$time_lineage <- compare_time_lineage(
      true_tree_time = true_tree_time,
      true_tree_div = true_tree_div,
      pred_tree = pred_tree,
      dist_prob = dist_prob,
      true_divisions = true_divisions,
      true_history = true_history,
      c0_name = c0_name,
      from_present = from_present,
      cor_method = cor_method
    )
  }

  if (length(out) == 0L) {
    stop("Provide `true_tree` and/or the full time-lineage reference inputs.")
  }

  out
}
