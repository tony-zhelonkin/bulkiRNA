#' Test gene-set coregulation with GESECA
#'
#' `gs_coregulation()` tests whether genes in each set vary together across
#' samples. Unlike [gs_test()], it takes a genes x samples expression matrix
#' and has no contrast. The result uses `contrast = "coregulation"` as the
#' required non-contrast grouping label in the [gs_result] contract.
#'
#' Rows are centered by default because GESECA's variance-along-a-direction
#' score assumes zero-mean gene profiles. Scaling remains off: unit-variance
#' scaling would change the question from coregulation in the supplied
#' expression units to equal contribution from every gene. Set `center = FALSE`
#' only when the rows have already been centered, for example after a PCA
#' reduction or for a CoReSh expression chunk.
#'
#' The statistic is unsigned, so `direction` is `NA` on every row. Filtering
#' this result to `direction = "up"`, `"down"`, or `"ns"` therefore returns no
#' rows; rank it by `stat`, `p_value`, or `padj` instead. Its
#' [summary.gs_result()] reports `n_up` and `n_down` as unavailable (`NA`).
#' Direction-restricted or direction-faceted plots are therefore not
#' meaningful; use the renderers' default `direction = "both"` and an
#' unsigned facet such as database or contrast.
#'
#' A `pct_var` result cannot enter [gs_to_master()] by default because the
#' master table's `nes` column is reserved for NES. Pass
#' `stat_as_nes = TRUE` only when deliberately asserting that percentage of
#' variance explained may populate that column, or keep the result out of the
#' master table.
#'
#' @param expr Numeric matrix of expression values, genes x samples. Rows must
#'   have unique, non-missing gene identifiers matching `db`; columns are
#'   samples. A [gs_matrix] is rejected because its rows are pathways, not
#'   genes.
#' @param db A `gs_db` (see `gsdb_*`), or a named list of `gs_db` objects; the
#'   list name becomes the `database` column.
#' @param center Logical. Center every gene across samples before testing.
#'   Defaults to `TRUE`, as required for a general expression matrix.
#' @param scale Logical. Scale every gene to unit variance before testing.
#'   Defaults to `FALSE` so genes retain their relative variation.
#' @param min_size,max_size Positive whole-number set-size bounds after
#'   intersecting each set with the rows of `expr`.
#' @param sample_size Positive whole number controlling the precision of the
#'   adaptive multilevel estimator. The default `101L` is fgsea's GESECA
#'   default; larger values improve precision at additional cost.
#' @param eps Non-negative finite p-value boundary passed to
#'   [fgsea::geseca()]. Use `0` to remove the boundary.
#' @param seed Finite whole-number RNG seed, or `NULL` to use the caller's RNG
#'   stream. The default `123L` matches the package's general fgsea adapter;
#'   CoReSh keeps its separate historical default of `1L`.
#' @param verbose Logical. Emit one progress message per database.
#' @return A [gs_result] with `method = "geseca"`,
#'   `stat_type = "pct_var"`, `stat` equal to percentage of total variance
#'   explained, `direction = NA`, and GESECA's `log2err` retained.
#' @examples
#' expr <- rbind(
#'   A = c(-2, -1, 1, 2),
#'   B = c(-1.8, -0.9, 0.9, 1.8),
#'   C = c(1, -1, 1, -1),
#'   D = c(-1, 1, -1, 1)
#' )
#' db <- gsdb_register(
#'   list(COREGULATED = c("A", "B"), SCRAMBLED = c("A", "C")),
#'   database = "demo", species = "Homo sapiens"
#' )
#' gs_coregulation(expr, db, min_size = 2L, max_size = 3L)
#' @export
gs_coregulation <- function(expr, db, center = TRUE, scale = FALSE,
                            min_size = 10L, max_size = 500L,
                            sample_size = 101L, eps = 1e-50,
                            seed = 123L, verbose = FALSE) {
  .gs_coregulation_validate_expr(expr)
  .gs_coregulation_logical(center, "center")
  .gs_coregulation_logical(scale, "scale")
  min_size <- .gs_coregulation_positive_integer(min_size, "min_size")
  max_size <- .gs_coregulation_positive_integer(max_size, "max_size")
  if (min_size > max_size) {
    stop("`min_size` must not exceed `max_size`.", call. = FALSE)
  }
  sample_size <- .gs_coregulation_positive_integer(sample_size, "sample_size")
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) ||
      !is.finite(eps) || eps < 0) {
    stop("`eps` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
       !is.finite(seed) || seed != floor(seed))) {
    stop("`seed` must be one finite whole number or `NULL`.", call. = FALSE)
  }
  .gs_coregulation_logical(verbose, "verbose")

  dbs <- .gs_db_list(db)
  parts <- lapply(dbs, function(d) {
    if (verbose) {
      message("Testing coregulation of ", length(d$sets),
              " sets from ", d$database)
    }
    out <- .gs_geseca(
      expr, d$sets,
      center = center, scale = scale,
      min_size = min_size, max_size = max_size,
      sample_size = sample_size, eps = eps, seed = seed
    )
    res <- gs_result(
      out,
      database = d$database,
      contrast = "coregulation",
      method = "geseca",
      stat_type = "pct_var"
    )
    pathway_names <- d$pathway_names
    if (!is.null(pathway_names)) {
      hit <- res$pathway_id %in% names(pathway_names)
      res$pathway_name[hit] <- unname(pathway_names[res$pathway_id[hit]])
    }
    res
  })
  do.call(rbind, parts)
}

#' Validate a GESECA expression matrix
#'
#' @param expr Object supplied as `expr` to [gs_coregulation()].
#' @return `expr`, invisibly.
#' @keywords internal
.gs_coregulation_validate_expr <- function(expr) {
  if (inherits(expr, "gs_matrix")) {
    stop("`expr` must be a gene x sample expression matrix, not a `gs_matrix`; ",
         "a `gs_matrix` has pathways in its rows.", call. = FALSE)
  }
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr` must be a numeric matrix of genes x samples.", call. = FALSE)
  }
  if (nrow(expr) < 2L || ncol(expr) < 2L) {
    stop("`expr` must contain at least two genes and two samples.",
         call. = FALSE)
  }
  genes <- rownames(expr)
  if (is.null(genes)) {
    stop("`expr` must have gene row names matching `db`.", call. = FALSE)
  }
  if (anyNA(genes) || any(!nzchar(genes))) {
    stop("`expr` gene row names must be non-missing and non-empty.",
         call. = FALSE)
  }
  if (anyDuplicated(genes)) {
    stop("`expr` has duplicated gene row names; aggregate or select one row ",
         "per gene before calling `gs_coregulation()`.", call. = FALSE)
  }
  if (any(!is.finite(expr))) {
    stop("`expr` must contain only finite values.", call. = FALSE)
  }
  invisible(expr)
}

#' Validate a logical GESECA argument
#'
#' @param x Value to validate.
#' @param name Argument name without backticks.
#' @return `x`, invisibly.
#' @keywords internal
.gs_coregulation_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  invisible(x)
}

#' Validate a positive whole-number GESECA argument
#'
#' @param x Value to validate.
#' @param name Argument name without backticks.
#' @return The value as an integer.
#' @keywords internal
.gs_coregulation_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop("`", name, "` must be a positive whole number.", call. = FALSE)
  }
  as.integer(x)
}

#' GESECA adapter
#'
#' @param expr Numeric genes x samples matrix.
#' @param sets Named list of gene sets.
#' @inheritParams gs_coregulation
#' @return A data frame of core [gs_result] columns plus `log2err`.
#' @keywords internal
.gs_geseca <- function(expr, sets, center, scale, min_size, max_size,
                       sample_size, eps, seed) {
  res <- .with_pinned_seed(
    seed,
    fgsea::geseca(
      pathways = sets,
      E = expr,
      minSize = min_size,
      maxSize = max_size,
      center = center,
      scale = scale,
      sampleSize = sample_size,
      eps = eps,
      nproc = 1L
    )
  )
  if (!is.data.frame(res)) {
    stop("GESECA returned a non-tabular result.", call. = FALSE)
  }
  expected <- c("pathway", "pctVar", "pval", "padj", "log2err", "size")
  missing <- setdiff(expected, names(res))
  if (length(missing)) {
    stop("GESECA result is missing column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!nrow(res)) {
    out <- .gs_empty_core(numeric_cols = "log2err")
    out$direction <- character(0)
    return(out)
  }
  pct_var <- as.numeric(res$pctVar)
  if (anyNA(pct_var) || any(!is.finite(pct_var)) || any(pct_var < 0)) {
    stop("GESECA returned a missing, non-finite, or negative `pctVar`; ",
         "check that `expr` retains non-zero variation after centering and ",
         "scaling.", call. = FALSE)
  }
  data.frame(
    pathway_id = as.character(res$pathway),
    pathway_name = as.character(res$pathway),
    n_genes = as.integer(lengths(sets[as.character(res$pathway)])),
    n_genes_tested = as.integer(res$size),
    stat = pct_var,
    direction = rep(NA_character_, nrow(res)),
    p_value = as.numeric(res$pval),
    padj = as.numeric(res$padj),
    log2err = as.numeric(res$log2err),
    stringsAsFactors = FALSE
  )
}
