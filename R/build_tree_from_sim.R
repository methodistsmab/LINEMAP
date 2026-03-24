#' Build `ape::phylo` objects from simulation output
#'
#' Converts edge tables and cell history from [simulate_lineage()] into one or
#' two trees (`phy` with branch lengths from recorded times, optional `phy1` with
#' unit edge lengths).
#'
#' @param sim Output list from [simulate_lineage()].
#' @param gRNA.num Number of gRNA positions (used for tip renaming).
#' @param which `"hg"` or `"sg"` lineage to use.
#' @param use_full_edges If `TRUE`, use full-state edge tables (`*_full`).
#' @param keep_alive_only If `TRUE`, prune edges with [alive_network_iterative()].
#' @param set_all_edge_time_to_one If `TRUE`, build a second tree with all edges length 1.
#' @param rename_cells If `TRUE`, map tip labels to `C1`, `C2`, ...
#' @param plot If `TRUE`, draw trees with [plot_time_scaled()].
#' @param plot_which Which tree(s) to plot: `"phy"`, `"phy1"`, `"both"`, or `"none"`.
#'
#' @return A list with `phy`, optional `phy1`, `alive.network`, `history`, and `diag`.
#' @export
build_tree_from_sim <- function(sim,
                                gRNA.num,
                                which = c("hg", "sg"),
                                use_full_edges = TRUE,
                                keep_alive_only = TRUE,
                                set_all_edge_time_to_one = FALSE,
                                rename_cells = TRUE,
                                plot = FALSE,
                                plot_which = c("phy", "phy1", "both", "none")) {
  which <- match.arg(which)

  hist_list <- NULL
  if (!is.null(sim$history[[which]])) {
    hist_list <- sim$history[[which]]
  } else if (which == "hg" && !is.null(sim$history$list_hg)) {
    hist_list <- sim$history$list_hg
  } else if (which == "sg" && !is.null(sim$history$list_sg)) {
    hist_list <- sim$history$list_sg
  } else {
    stop("Cannot find sim$history for: ", which)
  }

  names(hist_list) <- paste0("C", seq_along(hist_list))

  if (use_full_edges) {
    edge_mat <- if (which == "hg") sim$edges$hg_full else sim$edges$sg_full
  } else {
    edge_mat <- if (which == "hg") sim$edges$hg else sim$edges$sg
  }

  edge_df <- as.data.frame(edge_mat, stringsAsFactors = FALSE)
  colnames(edge_df) <- c("from", "to", "time")
  edge_df$time <- as.numeric(edge_df$time)

  alive_df <- edge_df
  if (keep_alive_only) {
    alive_mat <- alive_network_iterative(edge_mat, hist_list)
    colnames(alive_mat) <- c("from", "to", "time")
    alive_df <- as.data.frame(alive_mat, stringsAsFactors = FALSE)
    alive_df$time <- as.numeric(alive_df$time)
  }

  phy <- edges_to_phylo_from_table(alive_df)

  phy1 <- NULL
  if (set_all_edge_time_to_one) {
    alive_df_ones <- alive_df
    alive_df_ones$time <- 1
    phy1 <- edges_to_phylo_from_table(alive_df_ones)
  }

  if (rename_cells) {
    status_strings <- sapply(hist_list, function(status.tmp) {
      paste(status.tmp[1:(gRNA.num + 3), ncol(status.tmp)], collapse = "|")
    })
    cell_map <- stats::setNames(names(hist_list), status_strings)

    idx <- match(phy$tip.label, names(cell_map))
    phy$tip.label <- as.character(cell_map[idx])

    if (!is.null(phy1)) {
      idx <- match(phy1$tip.label, names(cell_map))
      phy1$tip.label <- as.character(cell_map[idx])
    }
  }

  plot_which <- match.arg(plot_which)

  if (plot) {
    if (plot_which %in% c("phy", "both")) {
      plot_time_scaled(phy)
      title("Time-scaled tree (phy)")
    }

    if (!is.null(phy1) && plot_which %in% c("phy1", "both")) {
      plot_time_scaled(phy1)
      title("Unit-length tree (phy1)")
    }
  }

  diag <- list(
    n_cell_division_hist = table(unlist(lapply(hist_list, ncol))),
    n_tips = length(phy$tip.label),
    n_edges_alive = nrow(alive_df)
  )

  list(
    phy = phy,
    phy1 = phy1,
    alive.network = alive_df,
    history = hist_list,
    diag = diag
  )
}

alive_network_iterative <- function(network, list) {
  if (is.null(network) || nrow(network) == 0) return(network)
  if (ncol(network) < 2) stop("network needs at least from, to columns")

  gRNA.num <- nrow(list[[1]]) - 3
  status_strings <- sapply(list, function(x) {
    paste(x[1:(gRNA.num + 3), ncol(x)], collapse = " ")
  })

  iter <- 0L
  removed_edges <- 0L

  repeat {
    iter <- iter + 1L
    froms <- network[, 1]
    tos <- network[, 2]

    leaves <- setdiff(tos, froms)
    if (length(leaves) == 0) break

    bad_leaves <- setdiff(leaves, status_strings)
    if (length(bad_leaves) == 0) {
      break
    }

    keep_mask <- !(tos %in% bad_leaves)
    removed_edges <- removed_edges + sum(!keep_mask)
    network <- network[keep_mask, , drop = FALSE]

    if (nrow(network) == 0) break
  }

  attr(network, "prune_iterations") <- iter - 1L
  attr(network, "removed_edges") <- removed_edges
  network
}


edges_to_phylo_from_table <- function(network,
                                      trailing_numeric = 3,
                                      time_pos = 1,
                                      use_edge_time = TRUE) {
  stopifnot(all(c("from", "to") %in% names(network)))
  parent_raw <- as.character(network$from)
  child_raw <- as.character(network$to)

  P <- lapply(parent_raw, parse_node_string, trailing_numeric = trailing_numeric, time_pos = time_pos)
  C <- lapply(child_raw, parse_node_string, trailing_numeric = trailing_numeric, time_pos = time_pos)

  node_time <- new.env(parent = emptyenv())
  add_node <- function(node) {
    if (!exists(node$key, envir = node_time, inherits = FALSE)) {
      assign(node$key, node$time, envir = node_time)
    } else {
      old <- get(node$key, envir = node_time, inherits = FALSE)
      if (is.finite(node$time) && abs(old - node$time) > 1e-9) {
        warning(
          "Inconsistent times for node ", node$key, ": first=", old,
          ", later=", node$time, ". Using first."
        )
      }
    }
  }
  for (i in seq_along(P)) {
    add_node(P[[i]])
    add_node(C[[i]])
  }

  df <- unique(data.frame(
    parent = vapply(P, `[[`, "", "key"),
    child = vapply(C, `[[`, "", "key"),
    stringsAsFactors = FALSE
  ))

  if (use_edge_time) {
    if (!"time" %in% names(network)) stop("network is missing required column `time`.")
    tmp <- data.frame(
      parent = vapply(P, `[[`, "", "key"),
      child = vapply(C, `[[`, "", "key"),
      elen = as.numeric(network$time),
      stringsAsFactors = FALSE
    )
    tmp <- stats::aggregate(elen ~ parent + child, data = tmp, FUN = function(x) x[1])
    df <- merge(df, tmp, by = c("parent", "child"), all.x = TRUE, sort = FALSE)
    edge_lengths <- df$elen
  } else {
    get_time <- function(k) get(k, envir = node_time, inherits = FALSE)
    edge_lengths <- mapply(function(pa, ch) {
      max(0, get_time(ch) - get_time(pa))
    }, df$parent, df$child)
    df$elen <- edge_lengths
  }

  if (any(!is.finite(edge_lengths) | edge_lengths < -1e-9)) {
    stop("Invalid branch lengths (NA or negative).")
  }
  edge_lengths[edge_lengths < 0] <- 0

  all_nodes <- unique(c(df$parent, df$child))

  indeg <- stats::setNames(integer(length(all_nodes)), all_nodes)
  for (ch in df$child) indeg[ch] <- indeg[ch] + 1
  roots <- names(indeg[indeg == 0])

  if (length(roots) > 1) {
    tmin <- min(sapply(all_nodes, function(k) get(k, envir = node_time, inherits = FALSE)))
    super_root <- sprintf("SUPER_ROOT@%.9f", tmin)
    assign(super_root, tmin, envir = node_time)
    df <- rbind(
      data.frame(parent = super_root, child = roots, elen = 0, stringsAsFactors = FALSE),
      df
    )
    all_nodes <- unique(c(super_root, all_nodes))
  } else if (length(roots) == 0) {
    stop("No root found (possible cycle in graph).")
  }

  outdeg <- stats::setNames(integer(length(all_nodes)), all_nodes)
  for (pa in df$parent) outdeg[pa] <- outdeg[pa] + 1
  leaves <- names(outdeg[outdeg == 0])

  tip.label <- leaves

  Ntip <- length(leaves)
  tip_index <- stats::setNames(seq_len(Ntip), leaves)
  internal_nodes <- setdiff(all_nodes, leaves)
  get_time <- function(k) get(k, envir = node_time, inherits = FALSE)
  internal_nodes <- internal_nodes[order(sapply(internal_nodes, get_time), decreasing = FALSE)]
  int_index <- stats::setNames(Ntip + seq_along(internal_nodes), internal_nodes)
  node_index <- c(tip_index, int_index)

  edge <- cbind(parent = node_index[df$parent], child = node_index[df$child])
  ord <- order(edge[, 1], edge[, 2])
  phy <- list(
    edge = edge[ord, , drop = FALSE],
    edge.length = df$elen[ord],
    tip.label = tip.label,
    Nnode = length(internal_nodes)
  )
  class(phy) <- "phylo"
  phy <- reorder.phylo(phy, "postorder")

  node_time_vec <- stats::setNames(
    sapply(c(leaves, internal_nodes), get_time),
    c(leaves, internal_nodes)
  )
  attr(phy, "node_time") <- node_time_vec

  phy
}
