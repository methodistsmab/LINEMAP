#' Convert an observed state to a dictionary label
#'
#' Internal helper for indexing missing-aware distance dictionaries.
#'
#' @param x Observed state value.
#' @param missing_label Label used for missing observations.
#'
#' @return Character scalar used as row/column label.
#'
#' @noRd
.state_to_label <- function(x, missing_label = "MISSING") {
  if (is.na(x)) {
    missing_label
  } else {
    as.character(x)
  }
}

#' Get the dictionary for one barcode position
#'
#' Internal helper that supports either a shared matrix dictionary or a list of
#' per-position dictionaries.
#'
#' @param tabs_list Matrix dictionary or list of matrix dictionaries.
#' @param k Barcode position.
#' @param use_list Logical; whether `tabs_list` is a list.
#'
#' @return A matrix dictionary for position `k`.
#'
#' @noRd
.get_distance_tab <- function(tabs_list, k, use_list) {
  if (use_list) {
    tabs_list[[k]]
  } else {
    tabs_list
  }
}

#' Validate distance dictionaries
#'
#' Internal helper for checking matrix/list dictionary inputs before distance
#' calculation.
#'
#' @param tabs_list Matrix dictionary or list of matrix dictionaries.
#' @param n_positions Number of barcode positions to use.
#'
#' @return Logical scalar indicating whether `tabs_list` is a list.
#'
#' @noRd
.validate_distance_tabs <- function(tabs_list, n_positions) {
  if (is.list(tabs_list)) {
    if (length(tabs_list) < n_positions) {
      stop("When `tabs_list` is a list, its length must be at least the number of positions used.")
    }
    for (k in seq_len(n_positions)) {
      if (!is.matrix(tabs_list[[k]])) {
        stop(sprintf("`tabs_list[[%d]]` is not a matrix.", k))
      }
    }
    return(TRUE)
  }

  if (is.matrix(tabs_list)) {
    return(FALSE)
  }

  stop("`tabs_list` must be a matrix or a list of matrices.")
}

#' Check whether a state label exists in a dictionary
#'
#' Internal helper for clearer error messages during dictionary lookup.
#'
#' @param tabk Matrix dictionary.
#' @param x_lab Row label.
#' @param y_lab Column label.
#' @param k Barcode position.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @noRd
.check_dictionary_labels <- function(tabk, x_lab, y_lab, k) {
  if (!x_lab %in% rownames(tabk)) {
    stop(sprintf("State `%s` is not in rownames of dictionary %d.", x_lab, k))
  }
  if (!y_lab %in% colnames(tabk)) {
    stop(sprintf("State `%s` is not in colnames of dictionary %d.", y_lab, k))
  }
  invisible(TRUE)
}

#' Compute a missing-aware pairwise probability distance matrix
#'
#' Computes pairwise distances by summing per-position dictionary lookups. Missing
#' sequence entries are mapped to `missing_label`, so they contribute a modeled
#' distance instead of being skipped.
#'
#' @param sequence Cell-by-position matrix or data frame of observed barcode
#'   states. Missing values should be `NA`.
#' @param tabs_list A matrix dictionary shared by all positions, or a list of
#'   per-position dictionaries produced by `build_pair_loglik_tables()`.
#' @param gRNA.num Number of barcode positions to use from `sequence`.
#' @param barcode.label Root barcode state used when adding the artificial `C0`.
#' @param add_C0 Logical; if `TRUE`, prepend a root row named `"C0"`.
#' @param missing_label Label used for missing observations in `tabs_list`.
#' @param show_progress Logical; if `TRUE`, show a text progress bar.
#'
#' @return A symmetric pairwise distance matrix.
#'
#' @export
compute_prob_distance <- function(sequence,
                                  tabs_list,
                                  gRNA.num = 200,
                                  barcode.label = 1,
                                  add_C0 = TRUE,
                                  missing_label = "MISSING",
                                  show_progress = TRUE) {
  sequence <- as.matrix(sequence)
  if (ncol(sequence) < gRNA.num) {
    stop("`ncol(sequence)` must be at least `gRNA.num`.")
  }

  use_list <- .validate_distance_tabs(tabs_list, gRNA.num)
  seq_use <- sequence[, seq_len(gRNA.num), drop = FALSE]

  if (is.null(rownames(seq_use))) {
    rownames(seq_use) <- as.character(seq_len(nrow(seq_use)))
  }

  if (add_C0) {
    seq_use <- rbind(rep(barcode.label, gRNA.num), seq_use)
    rownames(seq_use)[1] <- "C0"
  }

  n <- nrow(seq_use)
  dist_mat <- matrix(0.0, nrow = n, ncol = n)
  rownames(dist_mat) <- rownames(seq_use)
  colnames(dist_mat) <- rownames(seq_use)

  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    for (j in i:n) {
      s <- 0.0
      for (k in seq_len(gRNA.num)) {
        tabk <- .get_distance_tab(tabs_list, k, use_list)
        x_lab <- .state_to_label(seq_use[i, k], missing_label)
        y_lab <- .state_to_label(seq_use[j, k], missing_label)

        .check_dictionary_labels(tabk, x_lab, y_lab, k)
        s <- s + tabk[x_lab, y_lab]
      }
      dist_mat[i, j] <- s
      dist_mat[j, i] <- s
    }

    if (!is.null(pb)) {
      setTxtProgressBar(pb, i)
    }
  }

  dist_mat
}

#' Compute pairwise distances from an existing distance dictionary
#'
#' Lower-level helper for cases where the sequence has already been subset and
#' no artificial root row should be added automatically.
#'
#' @param sequence Cell-by-position matrix of observed barcode states.
#' @param distance A shared matrix dictionary or a list of per-position
#'   dictionaries.
#' @param missing_label Label used for missing observations.
#' @param show_progress Logical; if `TRUE`, show a text progress bar.
#'
#' @return A symmetric pairwise distance matrix.
#'
#' @export
calculate_distance_matrix <- function(sequence,
                                      distance,
                                      missing_label = "MISSING",
                                      show_progress = TRUE) {
  sequence <- as.matrix(sequence)
  n <- nrow(sequence)
  m <- ncol(sequence)

  use_list <- .validate_distance_tabs(distance, m)

  if (is.null(rownames(sequence))) {
    rownames(sequence) <- as.character(seq_len(n))
  }

  dist_mat <- matrix(0.0, nrow = n, ncol = n)
  rownames(dist_mat) <- rownames(sequence)
  colnames(dist_mat) <- rownames(sequence)

  pb <- NULL
  if (isTRUE(show_progress)) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    for (j in i:n) {
      s <- 0.0
      for (k in seq_len(m)) {
        tabk <- .get_distance_tab(distance, k, use_list)
        x_lab <- .state_to_label(sequence[i, k], missing_label)
        y_lab <- .state_to_label(sequence[j, k], missing_label)

        .check_dictionary_labels(tabk, x_lab, y_lab, k)
        s <- s + tabk[x_lab, y_lab]
      }
      dist_mat[i, j] <- s
      dist_mat[j, i] <- s
    }

    if (!is.null(pb)) {
      setTxtProgressBar(pb, i)
    }
  }

  dist_mat
}

#' Compute distance from a reference barcode
#'
#' Computes each cell's distance from either the root barcode or a selected
#' reference row in the sequence matrix.
#'
#' @param sequence Cell-by-position matrix of observed barcode states.
#' @param distance A shared matrix dictionary or a list of per-position
#'   dictionaries.
#' @param ref Optional reference row name or row index. If `NULL`, a root barcode
#'   of `root_state` is used.
#' @param root_state Root barcode state used when `ref = NULL`.
#' @param missing_label Label used for missing observations.
#'
#' @return Named numeric vector of distances from the reference.
#'
#' @export
calculate_distance_from_ref <- function(sequence,
                                        distance,
                                        ref = NULL,
                                        root_state = 1,
                                        missing_label = "MISSING") {
  sequence <- as.matrix(sequence)
  n <- nrow(sequence)
  m <- ncol(sequence)

  use_list <- .validate_distance_tabs(distance, m)

  if (is.null(rownames(sequence))) {
    rownames(sequence) <- as.character(seq_len(n))
  }

  if (is.null(ref)) {
    ref_seq <- rep(root_state, m)
  } else if (is.character(ref)) {
    ref_idx <- match(ref, rownames(sequence))
    if (is.na(ref_idx)) {
      stop("Could not find `ref` in `rownames(sequence)`.")
    }
    ref_seq <- sequence[ref_idx, ]
  } else if (is.numeric(ref) && length(ref) == 1L) {
    if (ref < 1L || ref > n) {
      stop("`ref` is outside the row range of `sequence`.")
    }
    ref_seq <- sequence[ref, ]
  } else {
    stop("`ref` must be `NULL`, a row name, or a row index.")
  }

  d_vec <- numeric(n)
  names(d_vec) <- rownames(sequence)

  for (i in seq_len(n)) {
    s <- 0.0
    for (k in seq_len(m)) {
      tabk <- .get_distance_tab(distance, k, use_list)
      x_lab <- .state_to_label(ref_seq[k], missing_label)
      y_lab <- .state_to_label(sequence[i, k], missing_label)

      .check_dictionary_labels(tabk, x_lab, y_lab, k)
      s <- s + tabk[x_lab, y_lab]
    }
    d_vec[i] <- s
  }

  d_vec
}
