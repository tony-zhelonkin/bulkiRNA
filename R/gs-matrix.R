#' The `gs_matrix` contract
#'
#' `gs_matrix` is the shape a per-sample scoring method returns: a numeric
#' matrix of **pathways (rows) x samples (columns)**, carrying enough metadata
#' to be tested or plotted without the caller re-supplying it.
#'
#' Attributes:
#' \describe{
#'   \item{`database`}{character. The `gsdb_*` provider the sets came from.}
#'   \item{`method`}{character, e.g. `"gsva"`.}
#'   \item{`score_type`}{character. What a cell is, e.g. `"gsva"`, `"ssgsea"`,
#'     `"zscore"`, `"plage"`. Renderers label the colour bar from this.}
#'   \item{`pathway_names`}{named character vector mapping `rownames()` to
#'     human-readable labels. Renderers use it for axis text.}
#'   \item{`sample_data`}{data frame of sample metadata, one row per column of
#'     the matrix, in column order -- this is what `gs_test()` builds a design
#'     from.}
#' }
#'
#' `gs_test()` dispatching on a `gs_matrix` is the GSVA -> limma pipeline; it
#' returns a [gs_result] with `stat_type = "t"`.
#'
#' @name gs_matrix-class
NULL

#' Construct a `gs_matrix`
#'
#' @param x A numeric matrix, pathways x samples, with row and column names.
#' @param database Character. Source database label.
#' @param method Character. Scoring method, e.g. `"gsva"`.
#' @param score_type Character. What a cell holds; defaults to `method`.
#' @param pathway_names Optional named character vector mapping row names to
#'   display labels. Defaults to the row names themselves.
#' @param sample_data Optional data frame of sample metadata with one row per
#'   column of `x`. Row names, if present, must match `colnames(x)`; otherwise
#'   the rows are taken to be in column order.
#' @return A `gs_matrix`.
#' @keywords internal
gs_matrix <- function(x,
                      database,
                      method,
                      score_type = method,
                      pathway_names = NULL,
                      sample_data = NULL) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix (pathways x samples).", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("`x` must have both row names (pathways) and column names (samples).",
         call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) {
    stop("`x` has duplicated pathway row names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) {
    stop("`x` has duplicated sample column names.", call. = FALSE)
  }
  for (nm in c("database", "method", "score_type")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1L || is.na(val)) {
      stop("`", nm, "` must be a length-1 character value.", call. = FALSE)
    }
  }

  pathway_names <- pathway_names %||% stats::setNames(rownames(x), rownames(x))
  if (is.null(names(pathway_names))) {
    stop("`pathway_names` must be a *named* character vector.", call. = FALSE)
  }
  unmapped <- setdiff(rownames(x), names(pathway_names))
  if (length(unmapped)) {
    pathway_names <- c(pathway_names, stats::setNames(unmapped, unmapped))
  }
  pathway_names <- as.character(pathway_names[rownames(x)])
  names(pathway_names) <- rownames(x)

  if (!is.null(sample_data)) {
    sample_data <- .gs_align_sample_data(sample_data, colnames(x))
  }

  new_gs_matrix(
    x,
    database = database, method = method, score_type = score_type,
    pathway_names = pathway_names, sample_data = sample_data
  )
}

#' Low-level `gs_matrix` constructor
#'
#' @inheritParams gs_matrix
#' @return A `gs_matrix`.
#' @keywords internal
new_gs_matrix <- function(x, database, method, score_type,
                          pathway_names, sample_data) {
  structure(
    x,
    database = database,
    method = method,
    score_type = score_type,
    pathway_names = pathway_names,
    sample_data = sample_data,
    class = c("gs_matrix", "matrix", "array")
  )
}

#' Validate a `gs_matrix`
#'
#' @param x Object to check.
#' @return `x`, invisibly on success; errors otherwise.
#' @keywords internal
validate_gs_matrix <- function(x) {
  if (!inherits(x, "gs_matrix")) {
    stop("`x` must be a gs_matrix; construct it with `gs_score()`.",
         call. = FALSE)
  }
  if (!is.numeric(x) || length(dim(x)) != 2L) {
    stop("`x` must be a 2-dimensional numeric gs_matrix.", call. = FALSE)
  }
  pn <- attr(x, "pathway_names")
  if (!identical(names(pn), rownames(x))) {
    stop("`pathway_names` must be named by, and ordered as, `rownames(x)`.",
         call. = FALSE)
  }
  sd <- attr(x, "sample_data")
  if (!is.null(sd) && nrow(sd) != ncol(x)) {
    stop("`sample_data` must have one row per column of the matrix.",
         call. = FALSE)
  }
  invisible(x)
}

#' Align sample metadata to matrix columns
#'
#' @param sample_data Data frame of sample metadata.
#' @param samples Character vector of column names to align to.
#' @return `sample_data` reordered to `samples`.
#' @keywords internal
.gs_align_sample_data <- function(sample_data, samples) {
  if (!is.data.frame(sample_data)) {
    stop("`sample_data` must be a data frame.", call. = FALSE)
  }
  rn <- rownames(sample_data)
  if (!is.null(rn) && all(samples %in% rn)) {
    return(sample_data[samples, , drop = FALSE])
  }
  if (nrow(sample_data) != length(samples)) {
    stop(
      "`sample_data` has ", nrow(sample_data), " rows but the matrix has ",
      length(samples), " samples, and its row names do not identify them.",
      call. = FALSE
    )
  }
  rownames(sample_data) <- samples
  sample_data
}

# ---- accessors --------------------------------------------------------------

#' Metadata accessors for a `gs_matrix`
#'
#' @param x A `gs_matrix`.
#' @return The corresponding attribute: a length-1 character for
#'   `gs_database()`, `gs_method()` and `gs_score_type()`; a named character
#'   vector for `gs_pathway_names()`; a data frame or `NULL` for
#'   `gs_sample_data()`.
#' @name gs_matrix-accessors
#' @keywords internal
NULL

#' @rdname gs_matrix-accessors
#' @keywords internal
gs_database <- function(x) attr(x, "database")

#' @rdname gs_matrix-accessors
#' @keywords internal
gs_method <- function(x) attr(x, "method")

#' @rdname gs_matrix-accessors
#' @keywords internal
gs_score_type <- function(x) attr(x, "score_type")

#' @rdname gs_matrix-accessors
#' @keywords internal
gs_pathway_names <- function(x) attr(x, "pathway_names")

#' @rdname gs_matrix-accessors
#' @keywords internal
gs_sample_data <- function(x) attr(x, "sample_data")

# ---- S3 methods -------------------------------------------------------------

#' Print a `gs_matrix`
#'
#' @param x A `gs_matrix`.
#' @param n Number of rows to show.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.gs_matrix <- function(x, n = 6L, ...) {
  cat(sprintf(
    "# gs_matrix: %d pathways x %d samples [%s / %s]\n",
    nrow(x), ncol(x), attr(x, "database"), attr(x, "score_type")
  ))
  sd <- attr(x, "sample_data")
  if (!is.null(sd)) {
    cat("# sample_data: ", paste(names(sd), collapse = ", "), "\n", sep = "")
  }
  show <- min(n, nrow(x))
  print(unclass(x)[seq_len(show), , drop = FALSE])
  if (nrow(x) > show) {
    cat(sprintf("# ... %d more pathways\n", nrow(x) - show))
  }
  invisible(x)
}

#' Summarise a `gs_matrix`
#'
#' @param object A `gs_matrix`.
#' @param ... Ignored.
#' @return A one-row tibble of dimensions and score range.
#' @export
summary.gs_matrix <- function(object, ...) {
  v <- as.numeric(unclass(object))
  tibble::tibble(
    database = attr(object, "database"),
    method = attr(object, "method"),
    score_type = attr(object, "score_type"),
    n_pathways = nrow(object),
    n_samples = ncol(object),
    min = min(v, na.rm = TRUE),
    median = stats::median(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE),
    n_na = sum(is.na(v))
  )
}

#' Subset a `gs_matrix`
#'
#' Keeps the class and carries `pathway_names` / `sample_data` along with the
#' subset. Drops to a plain vector or matrix when the result is no longer a
#' 2-dimensional matrix.
#'
#' @param x A `gs_matrix`.
#' @param i,j Row (pathway) and column (sample) subscripts.
#' @param ... Ignored.
#' @param drop Passed to the matrix method.
#' @return A `gs_matrix`, or a plain vector when `drop` collapses a dimension.
#'   A subscript matching no pathway or no sample yields an empty `gs_matrix`,
#'   not an error.
#' @export
`[.gs_matrix` <- function(x, i, j, ..., drop = TRUE) {
  out <- unclass(x)
  attributes(out) <- attributes(out)[c("dim", "dimnames")]
  # A character subscript that misses gives a raw "subscript out of bounds";
  # name the offenders instead, as the renderers' `samples` argument does.
  if (!missing(i) && is.character(i)) {
    .gs_check_subscript(i, rownames(out), "pathway")
  }
  if (!missing(j) && is.character(j)) {
    .gs_check_subscript(j, colnames(out), "sample")
  }
  out <- if (missing(i) && missing(j)) {
    out[, , drop = drop]
  } else if (missing(i)) {
    out[, j, drop = drop]
  } else if (missing(j)) {
    out[i, , drop = drop]
  } else {
    out[i, j, drop = drop]
  }
  if (length(dim(out)) != 2L) {
    return(out)
  }
  # This goes through the low-level constructor, not `gs_matrix()`: a subset is
  # well-formed by construction, whereas the full constructor rejects two things
  # that are legitimate here. Zero-extent results -- `gm[rownames(gm) %in% keep, ]`
  # where `keep` matches nothing -- aborted the script instead of yielding an
  # empty matrix, even though "an empty result is a valid answer" is the
  # contract on the `gs_result` side; and a repeated subscript (`gm[c(1, 1), ]`)
  # tripped the duplicate-row-name check.
  rn <- rownames(out) %||% character(0)
  cn <- colnames(out) %||% character(0)
  dimnames(out) <- list(rn, cn)
  pn <- attr(x, "pathway_names")[rn]
  names(pn) <- rn

  sd <- attr(x, "sample_data")
  new_gs_matrix(
    out,
    database = attr(x, "database"),
    method = attr(x, "method"),
    score_type = attr(x, "score_type"),
    pathway_names = pn,
    sample_data = if (is.null(sd)) NULL else sd[cn, , drop = FALSE]
  )
}

#' Report character subscripts that match nothing
#'
#' @param idx Character subscript.
#' @param have Available names.
#' @param what `"pathway"` or `"sample"`, used in the message.
#' @return `NULL`, invisibly; errors when `idx` has unmatched entries.
#' @keywords internal
.gs_check_subscript <- function(idx, have, what) {
  unknown <- setdiff(idx, have)
  if (length(unknown)) {
    stop("`i` selects ", length(unknown), " ", what,
         if (length(unknown) == 1L) "" else "s",
         " not found in the gs_matrix: ",
         paste(utils::head(unknown, 5L), collapse = ", "), ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Convert a `gs_matrix` to long form
#'
#' @param x A `gs_matrix`.
#' @param ... Ignored.
#' @return A tibble with `pathway_id`, `pathway_name`, `sample`, `score`,
#'   `database`, `method`, `score_type`, plus any `sample_data` columns.
#' @importFrom tibble as_tibble
#' @exportS3Method tibble::as_tibble
as_tibble.gs_matrix <- function(x, ...) {
  m <- unclass(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]
  out <- tibble::tibble(
    pathway_id = rep(rownames(m), times = ncol(m)),
    sample = rep(colnames(m), each = nrow(m)),
    score = as.numeric(m)
  )
  out$pathway_name <- unname(attr(x, "pathway_names")[out$pathway_id])
  out$database <- attr(x, "database")
  out$method <- attr(x, "method")
  out$score_type <- attr(x, "score_type")
  out <- out[, c(
    "pathway_id", "pathway_name", "sample", "score",
    "database", "method", "score_type"
  )]
  sd <- attr(x, "sample_data")
  if (!is.null(sd)) {
    extra <- sd[out$sample, setdiff(names(sd), names(out)), drop = FALSE]
    rownames(extra) <- NULL
    out <- dplyr::bind_cols(out, tibble::as_tibble(extra))
  }
  out
}
