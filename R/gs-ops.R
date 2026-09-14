#' Filter a `gs_result`
#'
#' The named-argument filter every downstream script re-wrote by hand. Every
#' argument is optional and `NULL` means "no constraint", so
#' `gs_filter(res, padj = 0.05, direction = "up")` reads as the sentence it is.
#'
#' @param x A [gs_result].
#' @param padj Numeric. Keep rows with `padj < padj`.
#' @param p_value Numeric. Keep rows with `p_value < p_value`.
#' @param stat Numeric. Keep rows with `abs(stat) >= stat`.
#' @param direction Character subset of `"up"`, `"down"`, `"ns"`.
#' @param database,contrast Character vectors of values to keep.
#' @param pathway_id Character vector of pathway identifiers to keep.
#' @param pattern Regular expression matched against `pathway_name` and
#'   `pathway_id`.
#' @param min_genes,max_genes Integer bounds on `n_genes_tested`.
#' @return A [gs_result] with the surviving rows.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_filter(res, padj = 0.5, direction = "up")
#' @export
gs_filter <- function(x,
                      padj = NULL,
                      p_value = NULL,
                      stat = NULL,
                      direction = NULL,
                      database = NULL,
                      contrast = NULL,
                      pathway_id = NULL,
                      pattern = NULL,
                      min_genes = NULL,
                      max_genes = NULL) {
  .gs_check_result(x)
  keep <- rep(TRUE, nrow(x))
  if (!is.null(padj)) keep <- keep & !is.na(x$padj) & x$padj < padj
  if (!is.null(p_value)) {
    keep <- keep & !is.na(x$p_value) & x$p_value < p_value
  }
  if (!is.null(stat)) keep <- keep & !is.na(x$stat) & abs(x$stat) >= stat
  if (!is.null(direction)) {
    bad <- setdiff(direction, c("up", "down", "ns"))
    if (length(bad)) {
      stop("`direction` must be \"up\", \"down\" or \"ns\"; got ",
           paste(encodeString(bad, quote = "\""), collapse = ", "), ".",
           call. = FALSE)
    }
    keep <- keep & x$direction %in% direction
  }
  if (!is.null(database)) keep <- keep & x$database %in% database
  if (!is.null(contrast)) keep <- keep & x$contrast %in% contrast
  if (!is.null(pathway_id)) keep <- keep & x$pathway_id %in% pathway_id
  if (!is.null(pattern)) {
    keep <- keep & (grepl(pattern, x$pathway_name, ignore.case = TRUE) |
                      grepl(pattern, x$pathway_id, ignore.case = TRUE))
  }
  if (!is.null(min_genes)) {
    keep <- keep & !is.na(x$n_genes_tested) & x$n_genes_tested >= min_genes
  }
  if (!is.null(max_genes)) {
    keep <- keep & !is.na(x$n_genes_tested) & x$n_genes_tested <= max_genes
  }
  x[keep, , drop = FALSE]
}

#' Take the top pathways of a `gs_result`
#'
#' Selection, not filtering: `gs_top()` decides *which rows are shown*, leaving
#' significance to the renderer's highlighting. That separation is the one the
#' old `gsea_dotplot()` documented and everything else forgot.
#'
#' @param x A [gs_result].
#' @param n Integer, how many rows to keep per group.
#' @param by What to rank by: `"padj"`, `"p_value"` (smallest first) or
#'   `"stat"` (largest `abs(stat)` first).
#' @param by_direction Logical. Take `n` up and `n` down separately, which is
#'   what a balanced dotplot needs.
#' @param per Character vector of grouping columns; `n` applies within each
#'   group. Defaults to `c("database", "contrast")`.
#' @return A [gs_result] with at most `n` rows per group (per direction, if
#'   `by_direction`), ordered by the ranking column.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_top(res, n = 1)
#' @export
gs_top <- function(x, n = 10L,
                   by = c("padj", "p_value", "stat"),
                   by_direction = FALSE,
                   per = c("database", "contrast")) {
  .gs_check_result(x)
  by <- match.arg(by)
  if (!length(per)) per <- character(0)
  missing_per <- setdiff(per, names(x))
  if (length(missing_per)) {
    stop("`per` names column(s) not in the result: ",
         paste(missing_per, collapse = ", "), ".", call. = FALSE)
  }

  score <- if (by == "stat") -abs(x$stat) else x[[by]]
  score[is.na(score)] <- Inf

  grp_cols <- per
  if (by_direction) grp_cols <- c(grp_cols, "direction")
  key <- if (length(grp_cols)) {
    group_factors <- lapply(
      unname(as.list(as.data.frame(x)[grp_cols])),
      factor,
      exclude = NULL
    )
    do.call(interaction, c(group_factors, drop = TRUE, sep = "\r"))
  } else {
    rep("", nrow(x))
  }

  idx <- unlist(lapply(unique(key), function(k) {
    rows <- which(key == k)
    rows[order(score[rows])][seq_len(min(n, length(rows)))]
  }), use.names = FALSE)
  idx <- idx[order(score[idx])]
  x[idx, , drop = FALSE]
}

#' Split a `gs_result` into a named list
#'
#' The inverse of the `rbind()` that pools contrasts.
#'
#' @param x A [gs_result].
#' @param by Character vector of columns to split on. Default
#'   `c("database", "contrast")`.
#' @param drop_empty Logical. Drop combinations with no rows.
#' @return A named list of [gs_result] objects, named by the `by` values joined
#'   with `"."`.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_split(res, by = "contrast")
#' @export
gs_split <- function(x, by = c("database", "contrast"), drop_empty = TRUE) {
  .gs_check_result(x)
  missing_by <- setdiff(by, names(x))
  if (length(missing_by)) {
    stop("`by` names column(s) not in the result: ",
         paste(missing_by, collapse = ", "), ".", call. = FALSE)
  }
  f <- interaction(as.data.frame(x)[by], drop = drop_empty, sep = ".")
  lapply(split(seq_len(nrow(x)), f), function(i) x[i, , drop = FALSE])
}

#' Extract leading-edge (or overlap) genes
#'
#' Replaces `get_pathway_genes()`, `get_pathway_genes_all()` and
#' `get_significant_pathways()`: one function over the normalized result,
#' instead of three over an S4 `gseaResult`. ORA results carry their overlap
#' genes in the same `leading_edge` column.
#'
#' @param x A [gs_result].
#' @param padj Numeric. Keep only pathways below this FDR; `NULL` keeps all.
#' @param top_n Integer. Keep only the top `top_n` pathways by `padj`.
#' @param unique_genes Logical. Return one pooled, de-duplicated character
#'   vector instead of a per-pathway list.
#' @return A named list of character vectors (names are `pathway_id`), or a
#'   single character vector when `unique_genes = TRUE`.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_leading_edge(res)
#' @export
gs_leading_edge <- function(x, padj = NULL, top_n = NULL,
                            unique_genes = FALSE) {
  .gs_check_result(x)
  if (is.null(x[["leading_edge"]])) {
    stop("This result carries no `leading_edge` column, so there are no genes ",
         "to extract. fgsea and ORA results have one; GSVA -> limma results ",
         "do not.", call. = FALSE)
  }
  if (!is.null(padj)) x <- gs_filter(x, padj = padj)
  if (!is.null(top_n)) {
    x <- gs_top(x, n = top_n, by = "padj", per = character(0))
  }
  out <- stats::setNames(
    lapply(x[["leading_edge"]], function(z) as.character(unlist(z))),
    x$pathway_id
  )
  if (unique_genes) {
    return(unique(unlist(out, use.names = FALSE)))
  }
  out
}

#' Error unless `x` is a `gs_result`
#'
#' @param x Object to check.
#' @return `x`, invisibly.
#' @keywords internal
.gs_check_result <- function(x) {
  if (!inherits(x, "gs_result")) {
    stop("`x` must be a `gs_result`, as returned by `gs_test()` or ",
         "`gs_coregulation()`; got ",
         paste(encodeString(class(x), quote = "\""), collapse = "/"), ".",
         call. = FALSE)
  }
  invisible(x)
}
