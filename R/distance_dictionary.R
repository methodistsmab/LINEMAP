#' Stable log-sum-exp
#'
#' Computes `log(sum(exp(v)))` without avoidable numerical overflow.
#'
#' @param v Numeric vector of log-values.
#'
#' @return A numeric scalar.
logsumexp <- function(v) {
  a <- max(v)
  if (!is.finite(a)) {
    return(a)
  }
  a + log(sum(exp(v - a)))
}

#' Validate a transition matrix
#'
#' Internal helper used before constructing pairwise likelihood dictionaries.
#'
#' @param P Square numeric transition matrix.
#'
#' @return The validated matrix with row and column names filled if missing.
.validate_transition_matrix <- function(P) {
  if (!is.matrix(P)) {
    stop("`P` must be a matrix.")
  }
  if (!is.numeric(P)) {
    stop("`P` must be numeric.")
  }
  if (nrow(P) != ncol(P)) {
    stop("`P` must be a square matrix.")
  }
  if (any(!is.finite(P))) {
    stop("`P` contains non-finite values.")
  }

  S <- nrow(P)
  if (is.null(rownames(P))) {
    rownames(P) <- as.character(seq_len(S))
  }
  if (is.null(colnames(P))) {
    colnames(P) <- rownames(P)
  }
  if (!identical(rownames(P), colnames(P))) {
    warning("`rownames(P)` and `colnames(P)` are not identical; downstream indexing uses column names.")
  }

  P
}

#' Matrix power helper
#'
#' Internal wrapper around `expm::%^%` so this file does not need to attach
#' the `expm` package globally.
#'
#' @param P Square numeric matrix.
#' @param power Non-negative numeric scalar.
#'
#' @return Matrix power `P ^ power`.
.matrix_power <- function(P, power) {
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("Package `expm` is required. Please install it before building distance dictionaries.")
  }
  expm::`%^%`(P, power)
}

#' Build one pairwise log-likelihood distance dictionary
#'
#' Builds a pairwise distance dictionary for one transition matrix. The output
#' includes one additional row and column for missing observations, so observed
#' states and missing states can be scored with the same lookup table.
#'
#' @param P Square transition matrix. Row and column names should be state labels.
#' @param division Number of divisions used to propagate the initial state.
#' @param idx0 Initial/root state index in `P`; usually `1`.
#' @param missing_label Label used for missing observations in the dictionary.
#'
#' @return A `(S + 1) x (S + 1)` matrix of negative log-likelihood distances,
#'   where `S = nrow(P)`. The final row and column are named `missing_label`.
#'
#' @keywords internal
.build_pair_loglik_single <- function(P,
                                      division,
                                      idx0 = 1,
                                      missing_label = "MISSING") {
  P <- .validate_transition_matrix(P)
  S <- nrow(P)

  if (!is.numeric(division) || length(division) != 1L || is.na(division) || division < 0) {
    stop("`division` must be a single non-negative number.")
  }
  if (!is.numeric(idx0) || length(idx0) != 1L || is.na(idx0) || idx0 < 1 || idx0 > S) {
    stop("`idx0` is out of range.")
  }
  idx0 <- as.integer(idx0)

  Pd <- .matrix_power(P, division)
  r <- as.vector(Pd[idx0, ])
  margin <- as.vector(P %*% r)

  states <- c(colnames(P), missing_label)
  LL <- matrix(NA_real_, nrow = S + 1L, ncol = S + 1L)
  rownames(LL) <- states
  colnames(LL) <- states

  logP <- matrix(-Inf, nrow = S, ncol = S)
  logP[P > 0] <- log(P[P > 0])
  rownames(logP) <- rownames(P)
  colnames(logP) <- colnames(P)

  logPd <- matrix(-Inf, nrow = S, ncol = S)
  logPd[Pd > 0] <- log(Pd[Pd > 0])
  rownames(logPd) <- rownames(P)
  colnames(logPd) <- colnames(P)

  log_margin <- rep(-Inf, S)
  log_margin[margin > 0] <- log(margin[margin > 0])

  log_from_root <- logPd[idx0, ]

  for (ix in seq_len(S)) {
    for (iy in seq_len(S)) {
      log_terms <- log_from_root + logP[, ix] + logP[, iy]
      LL[ix, iy] <- -logsumexp(as.numeric(log_terms))
    }
  }

  for (ix in seq_len(S)) {
    log_terms <- log_from_root + logP[, ix] + log_margin
    LL[ix, S + 1L] <- -logsumexp(as.numeric(log_terms))
  }
  LL[S + 1L, seq_len(S)] <- LL[seq_len(S), S + 1L]

  log_terms <- log_from_root + 2 * log_margin
  LL[S + 1L, S + 1L] <- -logsumexp(as.numeric(log_terms))

  diag(LL)[seq_len(S)] <- 0
  LL
}

#' Build pairwise log-likelihood distance dictionaries
#'
#' Builds missing-aware pairwise distance dictionaries from either one
#' transition matrix or a list of transition matrices. This is the recommended
#' dictionary builder for barcode data with possible missing observations.
#'
#' @param P_or_list A square transition matrix, or a list of square transition
#'   matrices. A list is typically used when each barcode position has its own
#'   transition matrix.
#' @param division Number of divisions used to propagate the initial state.
#' @param idx0 Initial/root state index in each transition matrix; usually `1`.
#' @param missing_label Label used for missing observations in each dictionary.
#'
#' @return If `P_or_list` is a matrix, returns one `(S + 1) x (S + 1)` distance
#'   dictionary. If `P_or_list` is a list, returns a list of dictionaries with
#'   matching names.
#'
#' @export
build_pair_loglik_tables <- function(P_or_list,
                                     division,
                                     idx0 = 1,
                                     missing_label = "MISSING") {
  if (is.matrix(P_or_list)) {
    return(.build_pair_loglik_single(
      P = P_or_list,
      division = division,
      idx0 = idx0,
      missing_label = missing_label
    ))
  }

  if (is.list(P_or_list)) {
    if (length(P_or_list) == 0L) {
      stop("`P_or_list` is an empty list.")
    }

    out <- lapply(seq_along(P_or_list), function(i) {
      Pi <- P_or_list[[i]]
      if (!is.matrix(Pi)) {
        stop(sprintf("Element %d in `P_or_list` is not a matrix.", i))
      }
      .build_pair_loglik_single(
        P = Pi,
        division = division,
        idx0 = idx0,
        missing_label = missing_label
      )
    })

    if (!is.null(names(P_or_list))) {
      names(out) <- names(P_or_list)
    }
    return(out)
  }

  stop("`P_or_list` must be a matrix or a list of matrices.")
}
