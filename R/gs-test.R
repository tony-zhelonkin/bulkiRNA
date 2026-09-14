#' Test gene sets
#'
#' One verb for every gene-set test in the package. The method is chosen by the
#' class of `x`, so the common cases need no `method` argument at all:
#'
#' \describe{
#'   \item{named numeric vector}{a ranked statistic (see [gs_ranks()]) ->
#'     preranked GSEA via [fgsea::fgseaMultilevel()], `stat_type = "NES"`.}
#'   \item{character vector}{a gene list -> over-representation analysis via
#'     [fgsea::fora()], `stat_type = "log2_fold_enrichment"`.}
#'   \item{`gs_matrix`}{per-sample scores (see [gs_score()]) -> limma over the
#'     score matrix, `stat_type = "t"`.}
#' }
#'
#' Every branch returns the same [gs_result]. Pooling across contrasts is
#' `rbind()` on those results, because `contrast` is a column -- there is no
#' separate pooled code path.
#'
#' @param x A named numeric vector, a character vector of gene identifiers, or
#'   a `gs_matrix`.
#' @param db A `gs_db` (see `gsdb_*`), or a named list of `gs_db` objects; the
#'   list name becomes the `database` column. Not used when `x` is a
#'   `gs_matrix`, whose sets were already applied by [gs_score()].
#' @param method Optional character, `"fgsea"`, `"ora"` or `"limma"`. Only
#'   needed to override the class-based choice.
#' @param contrast Character label recorded in the `contrast` column.
#' @param ... Passed to the adapter: see [gs_test_fgsea_params] for the fgsea
#'   and ORA knobs, and the limma arguments `design`, `coef`, `block`.
#' @param verbose Logical. Emit progress messages.
#' @return A [gs_result].
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' ranks <- stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6])
#' gs_test(ranks, db, min_size = 1, max_size = 10)
#' @export
gs_test <- function(x, db, method = NULL, contrast = "contrast", ...,
                    verbose = FALSE) {
  UseMethod("gs_test")
}

#' @rdname gs_test
#' @export
gs_test.default <- function(x, db, method = NULL, contrast = "contrast", ...,
                            verbose = FALSE) {
  stop(
    "`gs_test()` has no method for an object of class ",
    paste(encodeString(class(x), quote = "\""), collapse = "/"),
    ". Supply a named numeric vector (fgsea), a character vector (ORA), ",
    "or a gs_matrix (limma).",
    call. = FALSE
  )
}

#' @rdname gs_test
#' @export
gs_test.numeric <- function(x, db, method = NULL, contrast = "contrast", ...,
                            verbose = FALSE) {
  method <- method %||% "fgsea"
  if (is.null(names(x))) {
    stop("A numeric `x` must be named by gene identifier; see `gs_ranks()`.",
         call. = FALSE)
  }
  .gs_over_dbs(db, function(sets, database) {
    switch(
      method,
      fgsea = .gs_fgsea(x, sets, ...),
      stop("`method` must be \"fgsea\" for a ranked numeric `x`; got ",
           encodeString(method, quote = "\""), ".", call. = FALSE)
    ) |>
      gs_result(database = database, contrast = contrast,
                method = "fgsea", stat_type = "NES")
  }, verbose = verbose)
}

#' @rdname gs_test
#' @export
gs_test.integer <- function(x, db, method = NULL, contrast = "contrast", ...,
                            verbose = FALSE) {
  gs_test.numeric(as.numeric(x) |> stats::setNames(names(x)), db,
                  method = method, contrast = contrast, ...,
                  verbose = verbose)
}

#' @rdname gs_test
#' @export
gs_test.character <- function(x, db, method = NULL, contrast = "contrast", ...,
                              verbose = FALSE) {
  method <- method %||% "ora"
  if (!identical(method, "ora")) {
    stop("`method` must be \"ora\" for a character `x`; got ",
         encodeString(method, quote = "\""),
         ".", call. = FALSE)
  }
  .gs_over_dbs(db, function(sets, database) {
    .gs_ora(x, sets, ...) |>
      gs_result(database = database, contrast = contrast,
                method = "ora", stat_type = "log2_fold_enrichment")
  }, verbose = verbose)
}

#' @rdname gs_test
#' @export
gs_test.gs_matrix <- function(x, db = NULL, method = NULL,
                              contrast = "contrast", ..., verbose = FALSE) {
  method <- method %||% "limma"
  if (!identical(method, "limma")) {
    stop("`method` must be \"limma\" for a gs_matrix `x`; got ",
         encodeString(method, quote = "\""),
         ".", call. = FALSE)
  }
  out <- .gs_limma(x, contrast = contrast, ...)
  gs_result(
    out,
    database = attr(x, "database"), contrast = out$contrast[1L] %||% contrast,
    # `gs_score()` records which of gsva/ssgsea/zscore/plage produced the
    # scores; hard-coding "gsva" here mislabelled three of the four in the
    # result's own provenance -- and that column is what gets exported.
    method = gs_method(x) %||% "gsva", stat_type = "t"
  )
}

# ---- db plumbing ------------------------------------------------------------

#' Normalise the `db` argument to a named list of gene-set lists
#'
#' Accepts a single `gs_db` or a named list of them, and returns a list of
#' `list(sets =, database =, pathway_names =)`.
#'
#' @param db A `gs_db` or named list of `gs_db`.
#' @return A list of one entry per database.
#' @keywords internal
.gs_db_list <- function(db) {
  if (inherits(db, "gs_db")) {
    return(list(list(
      sets = .gs_db_sets(db),
      database = attr(db, "database") %||% "database",
      pathway_names = attr(db, "pathway_names")
    )))
  }
  if (is.list(db) && length(db) && all(vapply(db, inherits, logical(1L), "gs_db"))) {
    # `%||%` only fires when `names(db)` is wholly NULL. A *partially* named
    # list returns "" for the unnamed slots, so the fallback never ran and half
    # the rows carried `database = ""` -- which then groups plots and
    # `gs_top(per = "database")` under a blank label.
    own <- vapply(db, function(d) attr(d, "database") %||% NA_character_,
                  character(1L), USE.NAMES = FALSE)
    nms <- names(db)
    if (is.null(nms)) nms <- rep(NA_character_, length(db))
    fill <- is.na(nms) | !nzchar(nms)
    nms[fill] <- own[fill]
    nms[is.na(nms)] <- "database"
    return(Map(function(d, nm) {
      list(
        sets = .gs_db_sets(d),
        database = nm,
        pathway_names = attr(d, "pathway_names")
      )
    }, db, nms))
  }
  stop("`db` must be a `gs_db` or a named list of `gs_db` objects. ",
       "Build one with `gsdb_load()`, `gsdb_msigdb()` or `gsdb_from_file()`.",
       call. = FALSE)
}

#' Strip the `gs_db` class down to the bare named list of character vectors
#'
#' @param db A `gs_db`.
#' @return A named list of character vectors.
#' @keywords internal
.gs_db_sets <- function(db) {
  sets <- unclass(db)
  attributes(sets) <- list(names = names(sets))
  sets
}

#' Run an adapter over one or several databases and row-bind the results
#'
#' @param db A `gs_db` or named list of `gs_db`.
#' @param fn Function of `(sets, database)` returning a [gs_result].
#' @param verbose Logical, emit per-database progress.
#' @return A [gs_result].
#' @keywords internal
.gs_over_dbs <- function(db, fn, verbose = FALSE) {
  dbs <- .gs_db_list(db)
  parts <- lapply(dbs, function(d) {
    if (verbose) message("Testing ", length(d$sets), " sets from ", d$database)
    res <- fn(d$sets, d$database)
    pn <- d$pathway_names
    if (!is.null(pn)) {
      hit <- res$pathway_id %in% names(pn)
      res$pathway_name[hit] <- unname(pn[res$pathway_id[hit]])
    }
    res
  })
  do.call(rbind, parts)
}

# ---- adapters ---------------------------------------------------------------

#' Shared fgsea / ORA parameters
#'
#' @param min_size,max_size Integer set-size bounds, applied after intersecting
#'   with the data. Defaults `10` and `500` -- the values
#'   `clusterProfiler::GSEA()` forwarded to fgsea, kept so results are
#'   unchanged.
#' @param eps Numeric. fgsea's p-value accuracy floor; `0` means "no floor",
#'   the setting the legacy pipeline used.
#' @param n_perm_simple Integer, fgsea's `nPermSimple`.
#' @param score_type One of `"std"`, `"pos"`, `"neg"`.
#' @param seed Integer seed set immediately before the fgsea call, or `NULL`
#'   to leave the RNG untouched. When set, the generator is pinned as well.
#' @name gs_test_fgsea_params
#' @keywords internal
NULL

#' fgsea adapter
#'
#' @param ranks Named numeric vector, decreasing.
#' @param sets Named list of character vectors.
#' @inheritParams gs_test_fgsea_params
#' @return A data frame of core [gs_result] columns plus the fgsea-specific
#'   `log2err`, `es`, and `leading_edge` columns.
#' @keywords internal
.gs_fgsea <- function(ranks, sets, min_size = 10L, max_size = 500L,
                      eps = 0, n_perm_simple = 100000L,
                      score_type = c("std", "pos", "neg"), seed = 123L) {
  score_type <- match.arg(score_type)
  # `set.seed()` mutates the caller's global stream, so a script that seeds
  # itself and then bootstraps after a `gs_test()` call was silently drawing
  # from the seed-123 stream instead of its own -- once per database, at that.
  # Pin the generator too because `bplapply()` changes it inside tasks. Seed
  # for fgsea's benefit, then restore the caller's stream exactly, including
  # "there was no stream yet".
  res <- .with_pinned_seed(
    seed,
    fgsea::fgseaMultilevel(
      pathways = sets, stats = ranks,
      minSize = min_size, maxSize = max_size,
      eps = eps, nPermSimple = n_perm_simple, scoreType = score_type
    )
  )
  res <- as.data.frame(res)
  if (!nrow(res)) {
    return(.gs_empty_core(
      numeric_cols = c("es", "log2err"), list_cols = "leading_edge"
    ))
  }
  data.frame(
    pathway_id = res$pathway,
    pathway_name = res$pathway,
    n_genes = lengths(sets[res$pathway]),
    n_genes_tested = res$size,
    stat = res$NES,
    p_value = res$pval,
    padj = stats::p.adjust(res$pval, method = "BH"),
    es = res$ES,
    log2err = res$log2err,
    stringsAsFactors = FALSE
  ) |>
    .gs_attach_list_col("leading_edge", res$leadingEdge)
}

#' ORA adapter
#'
#' @param genes Character vector of query genes.
#' @param sets Named list of character vectors.
#' @param universe Character vector of background genes. Defaults to the union
#'   of `genes` and every gene in `sets` -- state it explicitly whenever the
#'   real background is the set of detected genes, which it usually is.
#' @inheritParams gs_test_fgsea_params
#' @return A data frame of core [gs_result] columns.
#' @keywords internal
.gs_ora <- function(genes, sets, universe = NULL,
                    min_size = 10L, max_size = 500L) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  universe <- universe %||% union(genes, unlist(sets, use.names = FALSE))
  universe <- unique(universe[!is.na(universe)])
  outside <- setdiff(genes, universe)
  if (length(outside)) {
    stop(length(outside), " gene(s) in `x` are absent from `universe`, ",
         "e.g. ", paste(utils::head(outside, 3L), collapse = ", "),
         ". Fix `universe` rather than the gene list.", call. = FALSE)
  }
  res <- as.data.frame(fgsea::fora(
    pathways = sets, genes = genes, universe = universe,
    minSize = min_size, maxSize = max_size
  ))
  res <- res[res$overlap > 0L, , drop = FALSE]
  if (!nrow(res)) {
    return(.gs_empty_core(numeric_cols = "fold_enrichment",
                          integer_cols = "overlap",
                          list_cols = "leading_edge"))
  }
  fe <- if (!is.null(res$foldEnrichment)) {
    res$foldEnrichment
  } else {
    (res$overlap / length(genes)) / (res$size / length(universe))
  }
  data.frame(
    pathway_id = res$pathway,
    pathway_name = res$pathway,
    n_genes = lengths(sets[res$pathway]),
    n_genes_tested = res$size,
    stat = log2(fe),
    p_value = res$pval,
    padj = res$padj,
    fold_enrichment = fe,
    overlap = res$overlap,
    stringsAsFactors = FALSE
  ) |>
    .gs_attach_list_col("leading_edge", res$overlapGenes)
}

#' limma-over-scores adapter
#'
#' The GSVA half of the GSVA -> limma pipeline, run on a [gs_matrix] as if its
#' pathway scores were genes.
#'
#' @param x A `gs_matrix`.
#' @param design A model matrix with one row per sample, or a one-sided formula
#'   evaluated in the matrix's `sample_data`.
#' @param contrast A contrast specification: a string such as `"KO-WT"` passed
#'   to [limma::makeContrasts()], a numeric contrast vector, or a plain label
#'   when `coef` is used instead.
#' @param coef Coefficient of `design` to test when `contrast` is only a label.
#'   Defaults to the last column.
#' @param block Optional blocking factor for [limma::duplicateCorrelation()].
#' @param trend,robust Passed to [limma::eBayes()].
#' @return A data frame of core [gs_result] columns.
#' @keywords internal
.gs_limma <- function(x, design = NULL, contrast = "contrast", coef = NULL,
                      block = NULL, trend = FALSE, robust = FALSE) {
  .require_pkg("limma", "Testing a gs_matrix",
               'BiocManager::install("limma")')
  sd <- attr(x, "sample_data")
  if (inherits(design, "formula")) {
    if (is.null(sd)) {
      stop("`design` is a formula but the gs_matrix carries no `sample_data` ",
           "to evaluate it in.", call. = FALSE)
    }
    design <- stats::model.matrix(design, data = sd)
  }
  if (is.null(design)) {
    stop("`design` is required: supply a model matrix, or a one-sided formula ",
         "over the gs_matrix's `sample_data` columns.", call. = FALSE)
  }
  if (nrow(design) != ncol(x)) {
    stop("`design` has ", nrow(design), " rows but the gs_matrix has ",
         ncol(x), " samples.", call. = FALSE)
  }

  m <- unclass(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]

  cor <- NULL
  if (!is.null(block)) {
    cor <- limma::duplicateCorrelation(m, design, block = block)$consensus
  }
  fit <- limma::lmFit(m, design, block = block, correlation = cor)

  contrast_label <- if (is.character(contrast) && length(contrast) == 1L) {
    contrast
  } else {
    "contrast"
  }
  cvec <- NULL
  if (is.numeric(contrast)) {
    cvec <- contrast
  } else if (is.character(contrast) && length(contrast) == 1L &&
             grepl("[-+]", contrast) &&
             all(.gs_contrast_terms(contrast) %in% colnames(design))) {
    cvec <- limma::makeContrasts(contrasts = contrast, levels = design)
  }
  if (!is.null(cvec)) {
    fit <- limma::contrasts.fit(fit, cvec)
    coef <- 1L
  } else if (is.null(coef)) {
    coef <- ncol(design)
  }
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")

  data.frame(
    pathway_id = rownames(tt),
    pathway_name = unname(attr(x, "pathway_names")[rownames(tt)]),
    contrast = contrast_label,
    n_genes = NA_integer_,
    n_genes_tested = NA_integer_,
    stat = tt$t,
    p_value = tt$P.Value,
    padj = tt$adj.P.Val,
    es = tt$logFC,
    stringsAsFactors = FALSE
  )
}

#' Split a contrast string into its design-column terms
#'
#' @param s A contrast string such as `"KO-WT"`.
#' @return Character vector of trimmed terms.
#' @keywords internal
.gs_contrast_terms <- function(s) {
  parts <- strsplit(s, "[-+]")[[1L]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

#' A zero-row core-column data frame
#'
#' An empty result is a valid answer; `NULL` is not.
#'
#' Adapters name their optional columns here so an empty result has the same
#' *shape* as a non-empty one from the same method. Otherwise a contrast that
#' happened to yield no pathways lost `leading_edge`, and a per-contrast branch
#' calling [gs_leading_edge()] got "carries no `leading_edge` column" -- a usage
#' error reported for what is simply an empty answer.
#'
#' @param numeric_cols,integer_cols,list_cols Optional column names to append,
#'   zero-length and of the named type.
#' @return A zero-row data frame with the core [gs_result] columns.
#' @keywords internal
.gs_empty_core <- function(numeric_cols = character(0),
                           integer_cols = character(0),
                           list_cols = character(0)) {
  out <- data.frame(
    pathway_id = character(0), pathway_name = character(0),
    n_genes = integer(0), n_genes_tested = integer(0),
    stat = numeric(0), p_value = numeric(0), padj = numeric(0),
    stringsAsFactors = FALSE
  )
  for (nm in numeric_cols) out[[nm]] <- numeric(0)
  for (nm in integer_cols) out[[nm]] <- integer(0)
  for (nm in list_cols) out[[nm]] <- list()
  out
}

#' Attach a list column to a data frame without `I()` surprises
#'
#' @param df A data frame.
#' @param name Column name.
#' @param value A list, one element per row.
#' @return `df` with the list column attached.
#' @keywords internal
.gs_attach_list_col <- function(df, name, value) {
  if (is.null(value)) {
    return(df)
  }
  df[[name]] <- lapply(value, function(z) as.character(unlist(z)))
  df
}
