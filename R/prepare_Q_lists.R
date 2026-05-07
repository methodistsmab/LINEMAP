#' Prepare transition matrix lists for hgRNA and sgRNA simulations
#'
#' Normalizes an input transition matrix and expands it into two per-gRNA
#' transition matrix lists. `Q.list` is used for the hgRNA-like model, where all
#' rows follow the transition matrix. `Q.list1` is used for the sgRNA-like model,
#' where only the initial-state row can mutate and all other rows are identity.
#'
#' @param Q Square transition matrix.
#' @param gRNA.num Integer. Number of gRNA/barcode positions.
#' @param mutation.rate Numeric vector of mutation rates assigned to the first
#'   row of the base transition matrix.
#' @param recycle_blocks Integer. Number of times to repeat `mutation.rate`
#'   blocks. Must satisfy `gRNA.num == length(mutation.rate) * recycle_blocks`.
#'
#' @return A list with normalized `Q`, hgRNA matrices `Q.list`, sgRNA matrices
#'   `Q.list1`, and the number of states `N`.
#'
#' @export
prepare_Q_lists <- function(Q,
                            gRNA.num = 200,
                            mutation.rate = seq(0.05, 0.25, 0.05),
                            recycle_blocks = 40) {
  Q <- row.normalize(Q)
  N <- nrow(Q)
  
  # 5 base Q's differing only in row 1 mutation rate
  Q.base <- vector("list", length(mutation.rate))
  for (k in seq_along(mutation.rate)) {
    Qk <- Q
    mr <- mutation.rate[k]
    Qk[1, -1] <- Qk[1, -1] / (sum(Qk[1, -1]) / mr)
    Qk[1, 1]  <- 1 - mr
    Q.base[[k]] <- Qk
  }
  
  # Your original code makes Q.list length = gRNA.num = 200 by repeating 5 patterns
  # recycle_blocks=40 -> 5*(40)=200
  stopifnot(gRNA.num == length(mutation.rate) * recycle_blocks)
  
  Q.list <- vector("list", gRNA.num)
  idx <- 1
  for (b in 1:recycle_blocks) {
    for (k in seq_along(mutation.rate)) {
      Q.list[[idx]] <- Q.base[[k]]
      idx <- idx + 1
    }
  }
  
  # Q.list1: identity except first row from Q.list[[i]] first row
  Q.list1 <- vector("list", gRNA.num)
  for (i in 1:gRNA.num) {
    Qi <- diag(N)
    Qi[1, ] <- Q.list[[i]][1, ]
    Q.list1[[i]] <- Qi
  }
  
  list(Q = Q, Q.list = Q.list, Q.list1 = Q.list1, N = N)
}

#' Row-normalize a matrix
#'
#' @param M Numeric matrix.
#'
#' @return Row-normalized matrix.
#'
#' @noRd
row.normalize <- function(M) {
  rs <- rowSums(M)
  rs[rs == 0] <- 1
  M / rs
}

#' Draw from a non-negative normal distribution
#'
#' @param n Number of observations.
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation of the normal distribution.
#'
#' @return Numeric vector with negative values truncated to zero.
#'
#' @noRd
rposnorm <- function(n, mean, sd) pmax(0, rnorm(n, mean, sd))


#' Parse a simulation node string
#'
#' @param s Node string containing barcode states followed by numeric metadata.
#' @param trailing_numeric Number of numeric metadata fields at the end of `s`.
#' @param time_pos Which trailing numeric field stores absolute time.
#'
#' @return A list with node `key` and absolute `time`.
#'
#' @noRd
parse_node_string <- function(s, trailing_numeric = 3, time_pos = 1) {
  toks <- strsplit(trimws(s), "\\s+")[[1]]
  n <- length(toks)
  if (n < trailing_numeric + 1) stop("节点字符串太短：", s)
  num_tail <- suppressWarnings(as.numeric(toks[(n - trailing_numeric + 1):n]))
  if (all(is.na(num_tail))) stop("尾部数字全为 NA，无法解析：", s)
  key_tokens <- toks[1:(n )]
  key <- paste(key_tokens, collapse = "|")
  t_abs <- as.numeric(num_tail[time_pos])
  list(key = key, time = t_abs)
}
