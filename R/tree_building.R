#' Validate a named distance matrix
#'
#' Internal helper for tree reconstruction. It checks that the distance matrix is
#' square, named, symmetric by labels, and contains the requested outgroup when
#' needed.
#'
#' @param distance_matrix Square numeric distance matrix.
#' @param root_outgroup Optional outgroup label used for rooting NJ trees.
#' @param require_outgroup Logical; whether `root_outgroup` must be present.
#'
#' @return The validated distance matrix.
.validate_distance_matrix_for_tree <- function(distance_matrix,
                                               root_outgroup = "C0",
                                               require_outgroup = TRUE) {
  d <- as.matrix(distance_matrix)

  if (!is.numeric(d)) {
    stop("`distance_matrix` must be numeric.")
  }
  if (nrow(d) != ncol(d)) {
    stop("`distance_matrix` must be a square matrix.")
  }
  if (is.null(rownames(d)) || is.null(colnames(d))) {
    stop("`distance_matrix` must have row names and column names.")
  }
  if (!identical(rownames(d), colnames(d))) {
    stop("Row names and column names of `distance_matrix` must be exactly identical.")
  }
  if (any(!is.finite(d))) {
    stop("`distance_matrix` contains non-finite values.")
  }
  if (any(d < 0)) {
    stop("`distance_matrix` contains negative distances.")
  }
  if (isTRUE(require_outgroup) && !root_outgroup %in% rownames(d)) {
    stop("`root_outgroup` is not present in `distance_matrix`.")
  }

  d
}

#' Build a tree from a distance matrix
#'
#' Reconstructs a tree from a named distance matrix using either neighbor joining
#' or UPGMA. For NJ, the tree can be rooted by an outgroup such as `"C0"` and the
#' outgroup can optionally be dropped after rooting. For UPGMA, the outgroup is
#' removed before clustering by default.
#'
#' @param distance_matrix Square named distance matrix.
#' @param method Tree reconstruction method, either `"NJ"` or `"UPGMA"`.
#' @param root_outgroup Outgroup label, usually `"C0"`.
#' @param drop_outgroup Logical; if `TRUE`, remove `root_outgroup` from the final
#'   tree/clustering.
#' @param resolve.root Logical passed to `ape::root()` for NJ rooting.
#'
#' @return A rooted `phylo` object for `"NJ"` or an `hclust` object for
#'   `"UPGMA"`.
#'
#' @export
build_tree_from_distance <- function(distance_matrix,
                                     method = c("NJ", "UPGMA"),
                                     root_outgroup = "C0",
                                     drop_outgroup = TRUE,
                                     resolve.root = TRUE) {
  method <- match.arg(method)
  d <- .validate_distance_matrix_for_tree(
    distance_matrix = distance_matrix,
    root_outgroup = root_outgroup,
    require_outgroup = !is.null(root_outgroup)
  )

  if (method == "NJ") {
    if (!requireNamespace("phangorn", quietly = TRUE)) {
      stop("Package `phangorn` is required for NJ tree reconstruction.")
    }
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package `ape` is required for NJ tree reconstruction.")
    }

    tree <- phangorn::NJ(stats::as.dist(d))

    if (!is.null(root_outgroup)) {
      tree <- ape::root(tree, outgroup = root_outgroup, resolve.root = resolve.root)
      if (isTRUE(drop_outgroup)) {
        tree <- ape::drop.tip(tree, root_outgroup)
      }
    }

    tree$tip.label <- as.character(tree$tip.label)
    return(tree)
  }

  if (method == "UPGMA") {
    keep <- rownames(d)
    if (!is.null(root_outgroup) && isTRUE(drop_outgroup)) {
      keep <- setdiff(keep, root_outgroup)
    }
    if (length(keep) < 2L) {
      stop("At least two samples are required to build a UPGMA tree.")
    }

    d_use <- d[keep, keep, drop = FALSE]
    tree <- stats::hclust(stats::as.dist(d_use), method = "average")
    return(tree)
  }

  stop("Unsupported tree reconstruction method.")
}

#' Build a tree from barcode sequences
#'
#' High-level user-facing wrapper for barcode tree reconstruction. Depending on
#' the inputs provided, it can start from transition matrices, precomputed
#' distance dictionaries, or a precomputed distance matrix.
#'
#' The full workflow is:
#' transition matrix/list -> distance dictionary -> distance matrix -> tree.
#'
#' @param sequence Cell-by-position matrix or data frame of observed barcode
#'   states. Missing values should be `NA`.
#' @param P_or_list Optional square transition matrix, or list of transition
#'   matrices. Used to build `tabs_list` when `tabs_list` and `distance_matrix`
#'   are not provided.
#' @param division Number of divisions used to build the distance dictionary.
#'   Required when `P_or_list` is provided.
#' @param tabs_list Optional matrix dictionary or list of per-position
#'   dictionaries from `build_pair_loglik_tables()`.
#' @param distance_matrix Optional precomputed distance matrix. If provided,
#'   dictionary and distance calculation are skipped.
#' @param gRNA.num Number of barcode positions to use from `sequence`.
#' @param barcode.label Root barcode state used when adding `"C0"`.
#' @param method Tree reconstruction method, either `"NJ"` or `"UPGMA"`.
#' @param add_C0 Logical; if `TRUE`, add root row `"C0"` before computing the
#'   distance matrix.
#' @param root_outgroup Outgroup label, usually `"C0"`.
#' @param drop_outgroup Logical; if `TRUE`, remove `root_outgroup` from the final
#'   tree/clustering.
#' @param resolve.root Logical passed to `ape::root()` for NJ rooting.
#' @param idx0 Initial/root state index in each transition matrix; usually `1`.
#' @param missing_label Label used for missing observations in `tabs_list`.
#' @param show_progress Logical; if `TRUE`, show progress while computing
#'   distances.
#'
#' @return A list with three components:
#' \describe{
#'   \item{tree}{The reconstructed tree.}
#'   \item{distance}{The distance matrix used to reconstruct the tree.}
#'   \item{dictionary}{The distance dictionary used to compute the matrix, or
#'     `NULL` if `distance_matrix` was provided directly.}
#' }
#'
#' @export
build_tree_from_sequence <- function(sequence = NULL,
                                     P_or_list = NULL,
                                     division = NULL,
                                     tabs_list = NULL,
                                     distance_matrix = NULL,
                                     gRNA.num = 200,
                                     barcode.label = 1,
                                     method = c("NJ", "UPGMA"),
                                     add_C0 = TRUE,
                                     root_outgroup = "C0",
                                     drop_outgroup = TRUE,
                                     resolve.root = TRUE,
                                     idx0 = 1,
                                     missing_label = "MISSING",
                                     show_progress = TRUE) {
  method <- match.arg(method)

  if (is.null(distance_matrix)) {
    if (is.null(sequence)) {
      stop("`sequence` is required when `distance_matrix` is not provided.")
    }

    if (is.null(tabs_list)) {
      if (is.null(P_or_list) || is.null(division)) {
        stop("Provide `distance_matrix`, `tabs_list`, or both `P_or_list` and `division`.")
      }
      if (!exists("build_pair_loglik_tables", mode = "function")) {
        stop("`build_pair_loglik_tables()` is required. Please source `distance_dictionary.R` first.")
      }

      tabs_list <- build_pair_loglik_tables(
        P_or_list = P_or_list,
        division = division,
        idx0 = idx0,
        missing_label = missing_label
      )
    }

    if (!exists("compute_prob_distance", mode = "function")) {
      stop("`compute_prob_distance()` is required. Please source `distance_matrix.R` first.")
    }

    distance_matrix <- compute_prob_distance(
      sequence = sequence,
      tabs_list = tabs_list,
      gRNA.num = gRNA.num,
      barcode.label = barcode.label,
      add_C0 = add_C0,
      missing_label = missing_label,
      show_progress = show_progress
    )
  }

  tree <- build_tree_from_distance(
    distance_matrix = distance_matrix,
    method = method,
    root_outgroup = root_outgroup,
    drop_outgroup = drop_outgroup,
    resolve.root = resolve.root
  )

  list(
    tree = tree,
    distance = distance_matrix,
    dictionary = tabs_list
  )
}
