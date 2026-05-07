#' Extract terminal barcode sequences from simulation histories
#'
#' Converts a list of terminal cell history matrices into a cell-by-gRNA sequence
#' matrix. Each row is one cell and each column is one gRNA/barcode position.
#' The returned sequence uses the last column of each cell history matrix.
#'
#' @param history A non-empty list of cell history matrices, such as
#'   `res$sim$history$list_hg`, `res$sim$history$list_sg`, or
#'   `res$tree$history` from `simulate_lineage_tree()`.
#' @param gRNA.num Optional number of gRNA/barcode rows to extract. If `NULL`,
#'   it is inferred as `nrow(history[[1]]) - n_extra_rows`.
#' @param n_extra_rows Number of non-barcode rows at the bottom of each history
#'   matrix. In the simulator this is `3`: next division time, quiescence end,
#'   and death time.
#'
#' @return A matrix with cells as rows and gRNA/barcode positions as columns.
#'
#' @export
sequence_from_history <- function(history,
                                  gRNA.num = NULL,
                                  n_extra_rows = 3) {
  if (!is.list(history) || length(history) == 0L) {
    stop("`history` must be a non-empty list of cell history matrices.")
  }
  if (!all(vapply(history, is.matrix, logical(1L)))) {
    stop("Every element of `history` must be a matrix.")
  }

  if (is.null(gRNA.num)) {
    gRNA.num <- nrow(history[[1]]) - n_extra_rows
  }
  if (!is.numeric(gRNA.num) || length(gRNA.num) != 1L || is.na(gRNA.num) || gRNA.num <= 0) {
    stop("`gRNA.num` must be a positive numeric scalar.")
  }
  gRNA.num <- as.integer(gRNA.num)

  n_rows <- vapply(history, nrow, integer(1L))
  if (any(n_rows < gRNA.num)) {
    stop("At least one history matrix has fewer rows than `gRNA.num`.")
  }

  sequence <- t(
    vapply(history, function(x) {
      x[seq_len(gRNA.num), ncol(x)]
    }, numeric(gRNA.num))
  )

  rownames(sequence) <- names(history)
  colnames(sequence) <- paste0("gRNA", seq_len(gRNA.num))
  sequence
}

#' Legacy alias for `sequence_from_history()`
#'
#' This alias is kept for compatibility with older analysis scripts.
#'
#' @param list A non-empty list of cell history matrices.
#'
#' @return A matrix with cells as rows and gRNA/barcode positions as columns.
#'
#' @export
sequence.creat <- function(list) {
  sequence_from_history(history = list)
}
