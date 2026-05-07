# =========================
# Build phylo from sim history + edge table
# =========================
build_tree_from_sim <- function(sim,
                                gRNA.num,
                                which = c("hg", "sg"),
                                use_full_edges = TRUE,   # TRUE -> network.cell.single.*  (has timers in from/to strings)
                                keep_alive_only = TRUE,  # TRUE -> run alive_network_iterative(...)
                                set_all_edge_time_to_one = FALSE, # reproduce your phy1 behavior
                                rename_cells = TRUE,     # map tip labels to C1,C2,...
                                plot = FALSE,
                                plot_which = c("phy", "phy1", "both", "none")) {
  
  which <- match.arg(which)
  
  # ---- pick history list ----
  # supports sim$history$hg / sim$history$sg (recommended)
  # also supports older sim$history$list_hg naming if you used that
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
  
  # ---- name cells like your original list/list1 ----
  names(hist_list) <- paste0("C", seq_along(hist_list))
  
  # ---- choose edge table ----
  if (use_full_edges) {
    edge_mat <- if (which == "hg") sim$edges$hg_full else sim$edges$sg_full
  } else {
    edge_mat <- if (which == "hg") sim$edges$hg else sim$edges$sg
  }
  
  # safety
  edge_df <- as.data.frame(edge_mat, stringsAsFactors = FALSE)
  colnames(edge_df) <- c("from", "to", "time")
  edge_df$time <- as.numeric(edge_df$time)
  
  # ---- alive network (optional) ----
  alive_df <- edge_df
  if (keep_alive_only) {
    # alive_network_iterative expects (network, list_of_cells)
    alive_mat <- alive_network_iterative(edge_mat, hist_list)
    colnames(alive_mat) <- c("from", "to", "time")
    alive_df <- as.data.frame(alive_mat, stringsAsFactors = FALSE)
    alive_df$time <- as.numeric(alive_df$time)
  }
  
  # ---- build phylo ----
  phy <- edges_to_phylo_from_table(alive_df)
  
  # optional: build "all edge length = 1" tree, like your phy1
  phy1 <- NULL
  if (set_all_edge_time_to_one) {
    alive_df_ones <- alive_df
    alive_df_ones$time <- 1
    phy1 <- edges_to_phylo_from_table(alive_df_ones)
  }
  
  # ---- rename tip labels: map tip (status string) -> C#
  # In your code: status_strings derived from list1 last column and collapse with "|"
  if (rename_cells) {
    status_strings <- sapply(hist_list, function(status.tmp) {
      paste(status.tmp[1:(gRNA.num + 3), ncol(status.tmp)], collapse = "|")
    })
    # map: status_string -> cell_id (C#)
    cell_map <- setNames(names(hist_list), status_strings)
    
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
  
  
  # ---- quick diagnostics you used ----
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



# 迭代版：保留到“叶结点 ⊆ status_strings ∪ from列” 为止
alive_network_iterative <- function(network, list) {
  # 允许 matrix / data.frame
  if (is.null(network) || nrow(network) == 0) return(network)
  if (ncol(network) < 2) stop("network 至少需两列：from, to")
  
  gRNA.num = nrow(list[[1]]) -3
  status_strings <- sapply(list, function(x) {
    paste(x[1:(gRNA.num+3), ncol(x)], collapse = " ")
  })
  
  iter <- 0L
  removed_edges <- 0L
  
  repeat {
    iter <- iter + 1L
    froms <- network[, 1]
    tos   <- network[, 2]
    
    # 当前叶结点：只在 to 出现、不在 from 出现
    leaves <- setdiff(tos, froms)
    if (length(leaves) == 0) break
    
    # 不在 status_strings 中的 “坏叶子”
    bad_leaves <- setdiff(leaves, status_strings)
    if (length(bad_leaves) == 0) {
      # 所有叶子都在 status_strings 中，满足终止条件
      break
    }
    
    # 删掉指向坏叶子的边
    keep_mask <- !(tos %in% bad_leaves)
    removed_edges <- removed_edges + sum(!keep_mask)
    network <- network[keep_mask, , drop = FALSE]
    
    # 网络可能为空
    if (nrow(network) == 0) break
  }
  
  attr(network, "prune_iterations") <- iter - 1L
  attr(network, "removed_edges")    <- removed_edges
  network
}

# 从边表构建严格按时间的 phylo
# 参数：
# - trailing_numeric = 3 ：节点字符串末尾有3个数字
# - time_pos = 1 ：第1个数字是绝对时间
# - use_edge_time = TRUE ：优先使用 network$time 作为分支长度
edges_to_phylo_from_table <- function(network,
                                      trailing_numeric = 3,
                                      time_pos = 1,
                                      use_edge_time = TRUE) {
  stopifnot(all(c("from","to") %in% names(network)))
  parent_raw <- as.character(network$from)
  child_raw  <- as.character(network$to)
  
  P <- lapply(parent_raw, parse_node_string, trailing_numeric = trailing_numeric, time_pos = time_pos)
  C <- lapply(child_raw,  parse_node_string, trailing_numeric = trailing_numeric, time_pos = time_pos)
  
  # 汇总所有节点的绝对时间（key -> time）
  node_time <- new.env(parent = emptyenv())
  add_node <- function(node) {
    if (!exists(node$key, envir = node_time, inherits = FALSE)) {
      assign(node$key, node$time, envir = node_time)
    } else {
      # 若重复出现且时间不一致，可在此加一致性检查/警告
      old <- get(node$key, envir = node_time, inherits = FALSE)
      if (is.finite(node$time) && abs(old - node$time) > 1e-9) {
        warning("同一节点多次出现且时间不一致：", node$key,
                "；首次=", old, "，新见=", node$time, "。采用首次。")
      }
    }
  }
  for (i in seq_along(P)) { add_node(P[[i]]); add_node(C[[i]]) }
  
  # 边表（去重）
  df <- unique(data.frame(
    parent = vapply(P, `[[`, "", "key"),
    child  = vapply(C, `[[`, "", "key"),
    stringsAsFactors = FALSE
  ))
  
  # 分支长度
  if (use_edge_time) {
    if (!"time" %in% names(network)) stop("network 缺少 `time` 列。")
    edge_lengths <- as.numeric(network$time)
    # 与 df 对齐（上面 unique 可能改变顺序；做一次按 (parent,child) 合并）
    tmp <- data.frame(parent = vapply(P, `[[`, "", "key"),
                      child  = vapply(C, `[[`, "", "key"),
                      elen   = as.numeric(network$time),
                      stringsAsFactors = FALSE)
    tmp <- aggregate(elen ~ parent + child, data = tmp, FUN = function(x) x[1])  # 保留首个
    df <- merge(df, tmp, by = c("parent","child"), all.x = TRUE, sort = FALSE)
    edge_lengths <- df$elen
  } else {
    # 兜底：用子时刻 - 父时刻
    get_time <- function(k) get(k, envir = node_time, inherits = FALSE)
    edge_lengths <- mapply(function(pa, ch) {
      max(0, get_time(ch) - get_time(pa))
    }, df$parent, df$child)
  }
  
  # 检查非负
  if (any(!is.finite(edge_lengths) | edge_lengths < -1e-9)) {
    stop("检测到非法分支长度（NA/负数）.")
  }
  edge_lengths[edge_lengths < 0] <- 0
  
  # 所有节点
  all_nodes <- unique(c(df$parent, df$child))
  
  # 根（入度=0）；若多根，添加超根（时间=全图最小时间，边长=0）
  indeg <- setNames(integer(length(all_nodes)), all_nodes)
  for (ch in df$child) indeg[ch] <- indeg[ch] + 1
  roots <- names(indeg[indeg == 0])
  
  if (length(roots) > 1) {
    tmin <- min(sapply(all_nodes, function(k) get(k, envir = node_time, inherits = FALSE)))
    super_root <- sprintf("SUPER_ROOT@%.9f", tmin)
    assign(super_root, tmin, envir = node_time)
    df <- rbind(data.frame(parent = super_root, child = roots, elen = 0, stringsAsFactors = FALSE),
                df)
    all_nodes <- unique(c(super_root, all_nodes))
  } else if (length(roots) == 0) {
    stop("未找到根（可能存在环）。")
  }
  
  # 叶（出度=0）
  outdeg <- setNames(integer(length(all_nodes)), all_nodes)
  for (pa in df$parent) outdeg[pa] <- outdeg[pa] + 1
  leaves <- names(outdeg[outdeg == 0])
  
  # 叶子标签（直接用 key；若太长可自定义缩短策略）
  tip.label <- leaves
  
  # 映射到 phylo 索引：叶=1..Ntip；内部按时间升序 -> Ntip+1..Ntip+Nnode
  Ntip <- length(leaves)
  tip_index <- setNames(seq_len(Ntip), leaves)
  internal_nodes <- setdiff(all_nodes, leaves)
  get_time <- function(k) get(k, envir = node_time, inherits = FALSE)
  internal_nodes <- internal_nodes[order(sapply(internal_nodes, get_time), decreasing = FALSE)]
  int_index <- setNames(Ntip + seq_along(internal_nodes), internal_nodes)
  node_index <- c(tip_index, int_index)
  
  # 组装 edge 与长度
  edge <- cbind(parent = node_index[df$parent], child = node_index[df$child])
  ord <- order(edge[,1], edge[,2])
  phy <- list(
    edge = edge[ord, , drop = FALSE],
    edge.length = df$elen[ord],
    tip.label = tip.label,
    Nnode = length(internal_nodes)
  )
  class(phy) <- "phylo"
  phy <- reorder.phylo(phy, "postorder")
  
  # 附带节点绝对时间（按 key 命名，便于后续检索/作图）
  node_time_vec <- setNames(sapply(c(leaves, internal_nodes), get_time), c(leaves, internal_nodes))
  attr(phy, "node_time") <- node_time_vec
  
  phy
}