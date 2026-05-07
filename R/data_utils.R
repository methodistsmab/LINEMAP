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
