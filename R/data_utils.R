#' @noRd
row.normalize <- function(M) {
  M <- as.matrix(M)
  rs <- rowSums(M)
  rs[!is.finite(rs) | rs == 0] <- 1
  sweep(M, 1, rs, "/")
}

#' @noRd
rposnorm <- function(n, mean, sd) {
  x <- numeric(n)
  for (i in seq_len(n)) {
    repeat {
      z <- stats::rnorm(1, mean = mean, sd = sd)
      if (is.finite(z) && z > 0) {
        x[i] <- z
        break
      }
    }
  }
  x
}

#' Load example transition matrix
#'
#' Loads the bundled example transition matrix stored in `inst/extdata/X.Rdata`.
#'
#' @return The object `X` from the example data file (a square transition matrix).
#' @export
load_example_Q <- function() {
  path <- system.file("extdata", "X.Rdata", package = "LINEMAP")
  if (path == "") {
    stop("Example data file `X.Rdata` not found.")
  }

  e <- new.env(parent = emptyenv())
  load(path, envir = e)

  if (!exists("X", envir = e)) {
    stop("Object `X` not found in example data.")
  }

  get("X", envir = e)
}
