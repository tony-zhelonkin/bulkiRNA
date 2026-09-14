#' Score gene sets per sample
#'
#' Turns an expression matrix into a pathway x sample [gs_matrix] with GSVA.
#' Feed the result to [gs_test()] to get a [gs_result] with
#' `stat_type = "t"` -- that pair is the GSVA -> limma pipeline, written once.
#'
#' GSVA is a `Suggests`; the call is guarded.
#'
#' @param expr A numeric matrix of expression values (genes x samples), rows
#'   named by the same identifier type as `db`. Log-scale values (voom logCPM,
#'   vst, logTPM) are what `kcdf = "Gaussian"` expects.
#' @param db A `gs_db`.
#' @param method One of `"gsva"`, `"ssgsea"`, `"zscore"`, `"plage"`. Selects
#'   the GSVA parameter object; it also becomes the matrix's `score_type`.
#' @param min_size,max_size Integer set-size bounds after intersecting with the
#'   rows of `expr`.
#' @param kcdf Kernel for the `"gsva"` method: `"Gaussian"` for log-scale
#'   continuous data, `"Poisson"` for integer counts.
#' @param sample_data Optional data frame of sample metadata, one row per
#'   column of `expr`; [gs_test()] builds its design from this.
#' @param verbose Logical. Passed to GSVA and used for progress messages.
#' @param ... Additional arguments passed to the GSVA parameter constructor.
#' @return A [gs_matrix].
#' @examples
#' \dontrun{
#' sc  <- gs_score(logcpm, db, method = "gsva", sample_data = meta)
#' res <- gs_test(sc, design = ~ 0 + group, contrast = "groupKO-groupWT")
#' }
#' @export
gs_score <- function(expr, db,
                     method = c("gsva", "ssgsea", "zscore", "plage"),
                     min_size = 10L, max_size = 500L,
                     kcdf = c("Gaussian", "Poisson", "none"),
                     sample_data = NULL, verbose = FALSE, ...) {
  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr` must be a numeric matrix of genes x samples.", call. = FALSE)
  }
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("`expr` must have gene row names and sample column names.",
         call. = FALSE)
  }
  dbs <- .gs_db_list(db)
  if (length(dbs) != 1L) {
    stop("`gs_score()` scores one database at a time; got ", length(dbs),
         ". Call it once per database and combine the matrices yourself.",
         call. = FALSE)
  }
  d <- dbs[[1L]]

  scores <- .gs_gsva(
    expr, d$sets, method = method,
    min_size = min_size, max_size = max_size, kcdf = kcdf,
    verbose = verbose, ...
  )

  gs_matrix(
    scores,
    database = d$database,
    method = method,
    score_type = method,
    pathway_names = d$pathway_names,
    sample_data = sample_data
  )
}

#' GSVA adapter
#'
#' GSVA 2.x uses a parameter object: build it, then call [GSVA::gsva()].
#'
#' @inheritParams gs_score
#' @param sets Named list of character vectors.
#' @return A numeric matrix, pathways x samples.
#' @keywords internal
.gs_gsva <- function(expr, sets, method = "gsva",
                     min_size = 10L, max_size = 500L, kcdf = "Gaussian",
                     verbose = FALSE, ...) {
  .require_pkg("GSVA", "`gs_score()`", 'BiocManager::install("GSVA")')
  param <- switch(
    method,
    gsva = GSVA::gsvaParam(
      expr, sets, minSize = min_size, maxSize = max_size, kcdf = kcdf, ...
    ),
    ssgsea = GSVA::ssgseaParam(
      expr, sets, minSize = min_size, maxSize = max_size, ...
    ),
    zscore = GSVA::zscoreParam(
      expr, sets, minSize = min_size, maxSize = max_size, ...
    ),
    plage = GSVA::plageParam(
      expr, sets, minSize = min_size, maxSize = max_size, ...
    ),
    stop("Unknown scoring `method`: ", encodeString(method, quote = "\""),
         ".", call. = FALSE)
  )
  out <- GSVA::gsva(param, verbose = verbose)
  if (!is.matrix(out)) out <- as.matrix(out)
  out
}
