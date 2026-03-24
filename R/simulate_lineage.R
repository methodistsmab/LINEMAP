#' Simulate lineage under dual barcode models
#'
#' Runs the main branching process for `hg` and `sg` barcode dynamics, recording
#' edge tables and per-cell state history.
#'
#' @param Q.list List of transition matrices (one per gRNA), full model.
#' @param Q.list1 List of transition matrices (single-lineage style).
#' @param N Number of discrete barcode states per gRNA.
#' @param t_max Maximum simulation time.
#' @param cell.num Initial number of cells.
#' @param gRNA.num Number of gRNA positions.
#' @param death_rate_per_time Rate parameter for death time exponential.
#' @param p_quiesce_enter Probability of entering quiescence at division.
#' @param quiesce_mean,quiesce_sd Parameters for quiescence duration (positive normal).
#' @param p_senesce Probability of senescence at division.
#' @param div_mean,div_sd Parameters for inter-division time (positive normal).
#' @param verbose If `TRUE`, prints progress messages.
#'
#' @return A list with `edges`, `terminal`, and `history` components.
#' @export
simulate_lineage <- function(Q.list,
                             Q.list1,
                             N,
                             t_max = 500,
                             cell.num = 1,
                             gRNA.num = 200,
                             death_rate_per_time = 1 / 240,
                             p_quiesce_enter = 0.20,
                             quiesce_mean = 24,
                             quiesce_sd = 6,
                             p_senesce = 0.01,
                             div_mean = 36,
                             div_sd = 8,
                             verbose = TRUE) {
  list_hg <- list()
  list_sg <- list()

  for (i in seq_len(cell.num)) {
    status.matrix <- matrix(1, nrow = gRNA.num + 3, ncol = 1)
    status.matrix[gRNA.num + 1, 1] <- rposnorm(1, div_mean, div_sd) * runif(1, 0, 1)
    status.matrix[gRNA.num + 2, 1] <- 0
    status.matrix[gRNA.num + 3, 1] <- rexp(1, rate = death_rate_per_time)

    if (status.matrix[gRNA.num + 1, 1] >= status.matrix[gRNA.num + 3, 1]) {
      status.matrix[gRNA.num + 3, 1] <- status.matrix[gRNA.num + 1, 1] + 1
    }

    list_hg[[length(list_hg) + 1]] <- status.matrix
    list_sg[[length(list_sg) + 1]] <- status.matrix
  }

  network.cell.hg <- matrix(NA, nrow = 0, ncol = 3)
  network.cell.single.hg <- matrix(NA, nrow = 0, ncol = 3)
  network.cell.sg <- matrix(NA, nrow = 0, ncol = 3)
  network.cell.single.sg <- matrix(NA, nrow = 0, ncol = 3)

  for (time in 0:t_max) {
    if (verbose && (time %% 10 == 0)) message("time = ", time, " | n_cells = ", length(list_hg))

    remove.ind <- integer(0)

    for (i in seq_len(length(list_hg))) {
      status.vector_hg <- list_hg[[i]][, ncol(list_hg[[i]])]
      status.vector_sg <- list_sg[[i]][, ncol(list_sg[[i]])]

      next_div_time <- status.vector_hg[gRNA.num + 1]
      quiesce_end <- status.vector_hg[gRNA.num + 2]
      death_time <- status.vector_hg[gRNA.num + 3]

      if (time >= death_time) {
        remove.ind <- c(remove.ind, i)
        next
      }

      if (quiesce_end > time || !is.finite(quiesce_end)) next

      if (time >= next_div_time) {
        remove.ind <- c(remove.ind, i)

        status.update_hg <- matrix(0, nrow = gRNA.num + 3, ncol = 2)
        status.update_sg <- matrix(0, nrow = gRNA.num + 3, ncol = 2)

        for (j in seq_len(gRNA.num)) {
          status.update_hg[j, ] <- sample(
            x = 1:N, size = 2, replace = TRUE,
            prob = Q.list[[j]][status.vector_hg[j], ]
          )
          status.update_sg[j, ] <- sample(
            x = 1:N, size = 2, replace = TRUE,
            prob = Q.list1[[j]][status.vector_sg[j], ]
          )
        }

        next_divs <- rposnorm(2, div_mean, div_sd) + time
        status.update_hg[gRNA.num + 1, ] <- next_divs
        status.update_sg[gRNA.num + 1, ] <- next_divs

        for (k in 1:2) {
          if (runif(1) < p_senesce) {
            q_end <- Inf
          } else if (runif(1) < p_quiesce_enter) {
            q_end <- time + rposnorm(1, quiesce_mean, quiesce_sd)
          } else {
            q_end <- 0
          }
          status.update_hg[gRNA.num + 2, k] <- q_end
          status.update_sg[gRNA.num + 2, k] <- q_end
        }

        death_new <- time + rexp(2, rate = death_rate_per_time)
        status.update_hg[gRNA.num + 3, ] <- death_new
        status.update_sg[gRNA.num + 3, ] <- death_new

        list_hg[[length(list_hg) + 1]] <- cbind(list_hg[[i]], status.update_hg[, 1])
        list_hg[[length(list_hg) + 1]] <- cbind(list_hg[[i]], status.update_hg[, 2])
        list_sg[[length(list_sg) + 1]] <- cbind(list_sg[[i]], status.update_sg[, 1])
        list_sg[[length(list_sg) + 1]] <- cbind(list_sg[[i]], status.update_sg[, 2])

        parent_hg <- list_hg[[i]]
        parent_sg <- list_sg[[i]]

        from_hg <- paste(parent_hg[1:gRNA.num, ncol(parent_hg)], collapse = " ")
        to1_hg <- paste(status.update_hg[1:gRNA.num, 1], collapse = " ")
        to2_hg <- paste(status.update_hg[1:gRNA.num, 2], collapse = " ")
        dt1 <- min(t_max, max(status.update_hg[gRNA.num + 1, 1], status.update_hg[gRNA.num + 2, 1])) - time
        dt2 <- min(t_max, max(status.update_hg[gRNA.num + 1, 2], status.update_hg[gRNA.num + 2, 2])) - time
        network.cell.hg <- rbind(network.cell.hg, c(from_hg, to1_hg, dt1))
        network.cell.hg <- rbind(network.cell.hg, c(from_hg, to2_hg, dt2))

        from_full_hg <- paste(parent_hg[, ncol(parent_hg)], collapse = " ")
        to1_full_hg <- paste(status.update_hg[, 1], collapse = " ")
        to2_full_hg <- paste(status.update_hg[, 2], collapse = " ")
        network.cell.single.hg <- rbind(network.cell.single.hg, c(from_full_hg, to1_full_hg, dt1))
        network.cell.single.hg <- rbind(network.cell.single.hg, c(from_full_hg, to2_full_hg, dt2))

        from_sg <- paste(parent_sg[1:gRNA.num, ncol(parent_sg)], collapse = " ")
        to1_sg <- paste(status.update_sg[1:gRNA.num, 1], collapse = " ")
        to2_sg <- paste(status.update_sg[1:gRNA.num, 2], collapse = " ")
        dt1s <- min(t_max, max(status.update_sg[gRNA.num + 1, 1], status.update_sg[gRNA.num + 2, 1])) - time
        dt2s <- min(t_max, max(status.update_hg[gRNA.num + 1, 2], status.update_hg[gRNA.num + 2, 2])) - time
        network.cell.sg <- rbind(network.cell.sg, c(from_sg, to1_sg, dt1s))
        network.cell.sg <- rbind(network.cell.sg, c(from_sg, to2_sg, dt2s))

        from_full_sg <- paste(parent_sg[, ncol(parent_sg)], collapse = " ")
        to1_full_sg <- paste(status.update_sg[, 1], collapse = " ")
        to2_full_sg <- paste(status.update_sg[, 2], collapse = " ")
        network.cell.single.sg <- rbind(network.cell.single.sg, c(from_full_sg, to1_full_sg, dt1s))
        network.cell.single.sg <- rbind(network.cell.single.sg, c(from_full_sg, to2_full_sg, dt2s))
      }
    }

    if (length(remove.ind) > 0) {
      list_hg <- list_hg[-remove.ind]
      list_sg <- list_sg[-remove.ind]
    }
  }

  colnames(network.cell.hg) <- c("from", "to", "dt")
  colnames(network.cell.single.hg) <- c("from", "to", "dt")
  colnames(network.cell.sg) <- c("from", "to", "dt")
  colnames(network.cell.single.sg) <- c("from", "to", "dt")

  terminal_hg <- lapply(list_hg, function(m) m[, ncol(m)])
  terminal_sg <- lapply(list_sg, function(m) m[, ncol(m)])

  list(
    edges = list(
      hg = network.cell.hg,
      hg_full = network.cell.single.hg,
      sg = network.cell.sg,
      sg_full = network.cell.single.sg
    ),
    terminal = list(hg = terminal_hg, sg = terminal_sg),
    history = list(list_hg = list_hg, list_sg = list_sg)
  )
}
