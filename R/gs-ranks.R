#' Build a ranked gene vector for gene-set testing
#'
#' Turns a differential-expression table into the named, decreasing-sorted
#' numeric vector that [gs_test()] feeds to fgsea. This is the whole of the old
#' `run_gsea()` preamble, split out so ranking is inspectable and testable
#' independently of the test itself.
#'
#' Missing and infinite values are dropped with a warning, because a rank
#' vector containing them is silently truncated by fgsea anyway.
#'
#' @param x A data frame of DE results, or an already-named numeric vector.
#' @param metric Character. Column of `x` to rank by. Default `"t"`.
#' @param genes Optional character vector of gene identifiers, or the name of a
#'   column of `x` holding them. Defaults to `rownames(x)`.
#' @param collapse How to resolve duplicated gene identifiers: `"max_abs"`
#'   keeps the row with the largest absolute metric, `"mean"` averages, and
#'   `"none"` (default) leaves duplicates in place, matching the legacy
#'   `run_gsea()` behaviour.
#' @return A named numeric vector, sorted decreasing.
#' @examples
#' de <- data.frame(t = c(3, -1, 2), row.names = c("A", "B", "C"))
#' gs_ranks(de)
#' @export
gs_ranks <- function(x,
                     metric = "t",
                     genes = NULL,
                     collapse = c("none", "max_abs", "mean")) {
  collapse <- match.arg(collapse)

  if (is.numeric(x) && !is.data.frame(x)) {
    if (is.null(names(x))) {
      stop("`x` must be a *named* numeric vector of gene-level statistics.",
           call. = FALSE)
    }
    v <- x
  } else {
    if (!is.data.frame(x)) {
      stop("`x` must be a data frame of DE results or a named numeric vector.",
           call. = FALSE)
    }
    if (!metric %in% names(x)) {
      stop("`metric` column \"", metric, "\" not found in `x`.", call. = FALSE)
    }
    ids <- if (is.null(genes)) {
      rownames(x)
    } else if (length(genes) == 1L && genes %in% names(x)) {
      as.character(x[[genes]])
    } else {
      as.character(genes)
    }
    if (is.null(ids)) {
      stop("`x` has no rownames and `genes` was not supplied, so genes cannot ",
           "be identified.", call. = FALSE)
    }
    if (length(ids) != nrow(x)) {
      stop("`genes` has ", length(ids), " entries but `x` has ", nrow(x),
           " rows.", call. = FALSE)
    }
    v <- stats::setNames(as.numeric(x[[metric]]), ids)
  }

  if (anyNA(v) || anyNA(names(v))) {
    warning("Dropping ", sum(is.na(v) | is.na(names(v))),
            " gene(s) with a missing rank metric or identifier.",
            call. = FALSE)
    v <- v[!is.na(v) & !is.na(names(v))]
  }
  v <- sort(v, decreasing = TRUE)
  if (any(is.infinite(v))) {
    warning("Dropping ", sum(is.infinite(v)),
            " gene(s) with an infinite rank metric.", call. = FALSE)
    v <- v[!is.infinite(v)]
  }

  if (collapse != "none" && anyDuplicated(names(v))) {
    v <- .gs_collapse_ranks(v, collapse)
    v <- sort(v, decreasing = TRUE)
  }
  v
}

#' Collapse duplicated gene identifiers in a rank vector
#'
#' @param v A named numeric vector.
#' @param how `"max_abs"` or `"mean"`.
#' @return A named numeric vector with unique names.
#' @keywords internal
.gs_collapse_ranks <- function(v, how) {
  sp <- split(unname(v), names(v))
  out <- vapply(sp, function(z) {
    if (how == "mean") mean(z) else z[which.max(abs(z))]
  }, numeric(1L))
  out
}
