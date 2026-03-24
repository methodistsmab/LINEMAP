#' Prepare per-gRNA transition matrix lists
#'
#' Builds `Q.list` (full row-wise dynamics) and `Q.list1` (single-lineage style:
#' only the first row varies per gRNA) from a base transition matrix `Q`.
#'
#' @param Q A square transition matrix (rows will be normalized).
#' @param gRNA.num Number of gRNA positions (must equal `length(mutation.rate) * recycle_blocks`).
#' @param mutation.rate Numeric vector of mutation rates for the first row pattern.
#' @param recycle_blocks Number of blocks; total gRNAs = `length(mutation.rate) * recycle_blocks`.
#'
#' @return A list with components `Q`, `Q.list`, `Q.list1`, and `N` (state space size).
#' @export
prepare_Q_lists <- function(Q,
                            gRNA.num = 200,
                            mutation.rate = seq(0.05, 0.25, 0.05),
                            recycle_blocks = 40) {
  Q <- row.normalize(Q)
  N <- nrow(Q)

  Q.base <- vector("list", length(mutation.rate))
  for (k in seq_along(mutation.rate)) {
    Qk <- Q
    mr <- mutation.rate[k]
    Qk[1, -1] <- Qk[1, -1] / (sum(Qk[1, -1]) / mr)
    Qk[1, 1] <- 1 - mr
    Q.base[[k]] <- Qk
  }

  stopifnot(gRNA.num == length(mutation.rate) * recycle_blocks)

  Q.list <- vector("list", gRNA.num)
  idx <- 1
  for (b in seq_len(recycle_blocks)) {
    for (k in seq_along(mutation.rate)) {
      Q.list[[idx]] <- Q.base[[k]]
      idx <- idx + 1
    }
  }

  Q.list1 <- vector("list", gRNA.num)
  for (i in seq_len(gRNA.num)) {
    Qi <- diag(N)
    Qi[1, ] <- Q.list[[i]][1, ]
    Q.list1[[i]] <- Qi
  }

  list(Q = Q, Q.list = Q.list, Q.list1 = Q.list1, N = N)
}

row.normalize <- function(M) {
  rs <- rowSums(M)
  rs[rs == 0] <- 1
  M / rs
}

rposnorm <- function(n, mean, sd) pmax(0, rnorm(n, mean, sd))


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

