# Internal helpers shared by the gs_plot_* renderers. Nothing here is
# exported and nothing here computes statistics: every value a renderer draws
# already exists on the gs_result / gs_matrix it was handed.

#' Colour ramp for signed statistics
#'
#' Blue-white-orange, kept as the default because it survives colour-vision
#' deficiency better than the red-ended alternatives. Measured as the CIE Lab
#' distance between the two ends under simulated deficiency, which is what
#' decides whether a reader can still tell up from down:
#'
#' | ramp | normal | deutan | protan | tritan | worst |
#' |---|---|---|---|---|---|
#' | blue-orange (this) | 103.1 | 103.6 | 90.4 | 85.4 | **85.4** |
#' | `heatmap3` RdYlBu | 113.4 | 95.9 | 72.1 | 123.5 | 72.1 |
#' | `heatmap2` RdBu | 110.2 | 93.2 | 67.2 | 119.8 | 67.2 |
#'
#' Red loses the most lightness under protanopia, so a red-ended ramp scores
#' worst exactly where it matters. Both alternatives remain available through
#' [bulki_palettes()]; pass `bulki_palettes("heatmap2")[c(5, 3, 1)]` or
#' `bulki_palettes("heatmap3")[c(5, 3, 1)]` as `palette` to use one.
#'
#' Keeping this default also means no already-published figure changes colour.
#'
#' @return A length-3 character vector: low, mid, high.
#' @keywords internal
.gs_diverging_colours <- function() {
  c(low = "#2166AC", mid = "#F7F7F7", high = "#B35806")
}

#' Is a display name a placeholder rather than a name?
#'
#' A `pathway_name` is meant to be the reader-facing label, so a renderer
#' treats it as authoritative. That breaks when a provider wrote a placeholder:
#' `coresh_derived_sets.gmt` ships `-` in every description field, and the
#' resulting legend read `-` five times. Anything with no letter or digit in it
#' cannot be a name, so it falls back to the machine id.
#'
#' Applied to already-cached results too, which is the point: a cache built by
#' an older provider must still render, without being rebuilt.
#'
#' @param x Character vector of candidate display names.
#' @return A logical vector, `TRUE` where `x` carries no usable name.
#' @keywords internal
.gs_placeholder_name <- function(x) {
  x <- as.character(x)
  is.na(x) | !grepl("[[:alnum:]]", x)
}

#' Check that an object is a usable `gs_result`
#'
#' @param x Object to check.
#' @param arg Name of the argument being checked, for the error message.
#' @return `x`, invisibly.
#' @keywords internal
.gs_plot_check_result <- function(x, arg = "x") {
  if (!inherits(x, "gs_result")) {
    stop("`", arg, "` must be a gs_result; got ",
         paste(class(x), collapse = "/"), ".", call. = FALSE)
  }
  invisible(x)
}

#' Resolve display labels for the `database` key
#'
#' `gs_result$database` is the stable snake_case registry key and must never be
#' shown to a reader. The display string is looked up, in order, from an
#' explicit `labels` map, a `database_label` column, or a `database_label`
#' attribute; unresolved keys fall back to the key itself.
#'
#' @param x A `gs_result`.
#' @param labels Named character vector mapping keys to display labels, or
#'   `NULL`.
#' @return A character vector, one label per row of `x`.
#' @keywords internal
.gs_database_label <- function(x, labels = NULL) {
  keys <- as.character(x[["database"]])
  if (is.null(labels)) {
    col <- x[["database_label"]]
    if (!is.null(col)) return(as.character(col))
    labels <- attr(x, "database_label")
  }
  if (is.null(labels)) return(keys)
  if (is.null(names(labels))) {
    stop("`database_labels` must be a named character vector mapping ",
         "`database` keys to display labels.", call. = FALSE)
  }
  out <- unname(as.character(labels[keys]))
  out[is.na(out)] <- keys[is.na(out)]
  out
}

#' Display form of the `direction` column
#'
#' @param direction Character vector of `"up"` / `"down"` / `"ns"`.
#' @return A factor with levels `Up`, `Down`, `NS`, dropping unused levels.
#' @keywords internal
.gs_direction_factor <- function(direction) {
  pretty <- c(up = "Up", down = "Down", ns = "NS")
  out <- unname(pretty[as.character(direction)])
  out[is.na(out)] <- "NS"
  factor(out, levels = intersect(c("Up", "Down", "NS"), unique(out)))
}

#' Leading-edge gene ratio
#'
#' `length(leading_edge) / n_genes` for fgsea results, `overlap / n_genes` for
#' ORA. Read off the object; nothing is recomputed.
#'
#' @param x A `gs_result`.
#' @return A numeric vector, or `NULL` when the object carries neither column.
#' @keywords internal
.gs_gene_ratio <- function(x) {
  n <- as.numeric(x[["n_genes"]])
  if (!is.null(x[["leading_edge"]])) {
    return(lengths(x[["leading_edge"]]) / n)
  }
  if (!is.null(x[["overlap"]])) {
    return(as.numeric(x[["overlap"]]) / n)
  }
  NULL
}

#' Symmetric limits for a diverging scale
#'
#' @param values Numeric vector.
#' @return A length-2 numeric vector centred on zero.
#' @keywords internal
.gs_symmetric_limits <- function(values) {
  m <- suppressWarnings(max(abs(values), na.rm = TRUE))
  if (!is.finite(m) || m == 0) m <- 1
  c(-m, m)
}

#' Build the plotting frame for a `gs_result`
#'
#' Selects rows, resolves display labels, and returns a plain data frame with
#' the columns every dot/bar/tile renderer needs. Selection is by `sort_by`
#' within `group_by`; ties keep the object's own row order.
#'
#' @param x A `gs_result`.
#' @param top_n Integer. Rows kept per selection group; `NULL` or `Inf` keeps
#'   all.
#' @param sort_by One of `"padj"`, `"p_value"`, `"stat"` (absolute magnitude),
#'   `"stat_signed"`.
#' @param direction One of `"both"`, `"up"`, `"down"`.
#' @param group_by Character vector of columns defining the selection groups,
#'   or `NULL`.
#' @param keep_ids Character vector of `pathway_id`s. When supplied, rows are
#'   restricted to these ids and `top_n` selection is skipped -- this is how
#'   `compare =` keeps every panel of a comparison grid on the same pathways.
#' @param highlight Numeric FDR threshold for the `significant` flag, or
#'   `NULL`.
#' @param wrap_width Soft character width for label wrapping.
#' @param strip_prefix Logical, passed to [format_pathway_name()].
#' @param database_labels Named character vector of display labels, or `NULL`.
#' @return A data frame with columns `pathway_id`, `label`, `stat`, `p_value`,
#'   `padj`, `neg_log_padj`, `direction`, `database`, `database_label`,
#'   `contrast`, `significant`, and `gene_ratio` when available.
#' @keywords internal
.gs_plot_frame <- function(x,
                           top_n = NULL,
                           sort_by = "padj",
                           direction = "both",
                           group_by = NULL,
                           keep_ids = NULL,
                           highlight = 0.05,
                           wrap_width = 50,
                           strip_prefix = TRUE,
                           database_labels = NULL) {
  .gs_plot_check_result(x)

  df <- data.frame(
    pathway_id = as.character(x[["pathway_id"]]),
    pathway_name = as.character(x[["pathway_name"]]),
    stat = as.numeric(x[["stat"]]),
    p_value = as.numeric(x[["p_value"]]),
    padj = as.numeric(x[["padj"]]),
    direction = as.character(x[["direction"]]),
    database = as.character(x[["database"]]),
    database_label = .gs_database_label(x, database_labels),
    contrast = as.character(x[["contrast"]]),
    stringsAsFactors = FALSE
  )
  ratio <- .gs_gene_ratio(x)
  if (!is.null(ratio)) df$gene_ratio <- as.numeric(ratio)

  if (!identical(direction, "both")) {
    if (!direction %in% c("up", "down")) {
      stop("`direction` must be one of \"both\", \"up\", \"down\"; got ",
           dQuote(direction), ".", call. = FALSE)
    }
    df <- df[!is.na(df$direction) & df$direction == direction, , drop = FALSE]
  }

  if (is.null(keep_ids)) {
    df <- .gs_select_top(df, top_n = top_n, sort_by = sort_by,
                         group_by = group_by)
  } else {
    df <- df[df$pathway_id %in% keep_ids, , drop = FALSE]
    df <- .gs_select_top(
      df, top_n = NULL, sort_by = sort_by, group_by = NULL
    )
  }

  # Two different pathways can format to the same display label -- the same set
  # name in two collections, or two ids that differ only in a part
  # `format_pathway_name()` drops. `label` becomes the y-axis factor, so equal
  # labels collapse onto ONE row and the bars/points stack on top of each other:
  # the figure silently shows n-1 of n results. Disambiguate the collisions only
  # (with the machine id, the one thing guaranteed unique), before wrapping, so
  # the appended text is wrapped too.
  # Formatted selectively, not unconditionally. format_pathway_name() is built
  # for ALL_CAPS_SNAKE ids and is not idempotent on prose: it turns "." into a
  # space and title-cases every word, so a real display name arrives mangled --
  # "GSE174808 · 40.1% var" came out as "Gse174808 · 40 1% Var". A provider that
  # supplied pathway_names, or a caller who passed real labels, has already
  # decided how they read. Only names still equal to their machine id are
  # formatted. gs_plot_running() has done this since 1.0.0; the bar, dot and
  # heatmap renderers had not.
  lbl <- as.character(df$pathway_name)
  raw <- .gs_placeholder_name(lbl) | lbl == df$pathway_id
  lbl[raw] <- format_pathway_name(df$pathway_id[raw],
                                  strip_prefix = strip_prefix)
  dup <- lbl %in% lbl[duplicated(lbl)]
  if (any(dup)) {
    lbl[dup] <- paste0(lbl[dup], " (", df$pathway_id[dup], ")")
  }
  df$label <- .gs_wrap_label(lbl, width = wrap_width)

  # Carry through any column the result has that the plotting contract does not
  # use. Nothing here plots them: they ride along so that `gs_save()`'s
  # same-stem `.tsv` keeps them. That is what makes an enriched result -- a
  # CoReSh table joined to its GEO titles, accessions and percent of variation
  # -- arrive beside the figure without the caller passing `data =` by hand.
  # The axis stays compact and the provenance is one lookup away.
  extra <- setdiff(names(x), names(df))
  for (nm in extra) {
    col <- x[[nm]][match(df$pathway_id, x[["pathway_id"]])]
    # A list column would break the y-axis factor if it ever reached an
    # aesthetic; .gs_save_tsv() collapses it on the way out.
    df[[nm]] <- col
  }

  # `padj == 0` is reachable whenever a p-value underflows -- routine for
  # permutation-free methods and large n -- and `-log10(0)` is `Inf`, which no
  # size or colour scale can map. ggplot2 drops such a point with no warning
  # and no error, so the pathway silently missing from the figure is the single
  # STRONGEST hit in the analysis. Clamp to the smallest non-zero `padj` in the
  # frame, which keeps the scale's spread meaningful, rather than to the double
  # floor, which would peg one point at -log10 ~= 308 and flatten the rest.
  pos <- df$padj[is.finite(df$padj) & df$padj > 0]
  padj_floor <- if (length(pos)) min(pos) else .Machine$double.xmin
  df$neg_log_padj <- -log10(pmax(df$padj, padj_floor))
  df$significant <- if (is.null(highlight)) {
    rep(FALSE, nrow(df))
  } else {
    !is.na(df$padj) & df$padj < highlight
  }
  df
}

#' Select the top rows of a plotting frame
#'
#' @param df A data frame from [.gs_plot_frame()].
#' @param top_n Integer or `NULL`.
#' @param sort_by One of `"padj"`, `"p_value"`, `"stat"`, `"stat_signed"`.
#' @param group_by Character vector of grouping columns, or `NULL`.
#' @return A data frame, sorted and truncated.
#' @keywords internal
.gs_select_top <- function(df, top_n = NULL, sort_by = "padj",
                           group_by = NULL) {
  allowed <- c("padj", "p_value", "stat", "stat_signed")
  if (!sort_by %in% allowed) {
    stop("`sort_by` must be one of ",
         paste(dQuote(allowed), collapse = ", "), "; got ",
         dQuote(sort_by), ".", call. = FALSE)
  }
  if (nrow(df) == 0L) return(df)

  key <- switch(
    sort_by,
    padj = df$padj,
    p_value = df$p_value,
    stat = -abs(df$stat),
    stat_signed = -df$stat
  )
  ord <- order(key, seq_len(nrow(df)), na.last = TRUE)
  df <- df[ord, , drop = FALSE]

  if (is.null(top_n) || !is.finite(top_n)) return(df)
  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1) {
    stop("`top_n` must be a single positive number, or NULL for all rows.",
         call. = FALSE)
  }
  top_n <- as.integer(top_n)

  if (is.null(group_by)) {
    return(utils::head(df, top_n))
  }
  grp <- interaction(df[, group_by, drop = FALSE], drop = TRUE)
  keep <- unlist(lapply(split(seq_len(nrow(df)), grp), utils::head, top_n),
                 use.names = FALSE)
  df[sort(keep), , drop = FALSE]
}

#' Order pathway labels for a categorical axis
#'
#' @param df A plotting frame.
#' @param by Column name whose value orders the labels.
#' @param decreasing Logical; `TRUE` puts the largest value at the top.
#' @return `df` with `label` converted to an ordered factor.
#' @keywords internal
.gs_order_labels <- function(df, by = "stat", decreasing = TRUE) {
  if (nrow(df) == 0L) return(df)
  score <- stats::ave(df[[by]], df$label, FUN = function(v) {
    mean(v, na.rm = TRUE)
  })
  lev <- unique(df$label[order(score, decreasing = !decreasing)])
  df$label <- factor(df$label, levels = lev)
  df
}

#' Diverging fill scale for a signed statistic
#'
#' @param name Legend title, normally `gs_stat_label(x)`.
#' @param limits Length-2 numeric limits.
#' @param palette Length-3 character vector: low, mid, high.
#' @return A ggplot2 scale.
#' @keywords internal
#' @importFrom scales squish
.gs_fill_scale <- function(name, limits, palette = .gs_diverging_colours()) {
  scale_fill_gradient2(
    low = palette[[1L]], mid = palette[[2L]], high = palette[[3L]],
    midpoint = 0, limits = limits, oob = scales::squish, name = name
  )
}

#' An empty placeholder plot
#'
#' Returned instead of erroring when a selection leaves nothing to draw -- a
#' zero-row `gs_result` is a valid answer, so a zero-row plot must be too.
#'
#' @param title Plot title.
#' @param subtitle Explanatory subtitle.
#' @return A ggplot object.
#' @keywords internal
.gs_empty_plot <- function(title = NULL, subtitle = "No pathways to plot") {
  # `ggplot()` with no data leaves `$data` a `waiver`, so `nrow(p$data)` is
  # NULL and any caller writing `if (nrow(p$data) > 0L)` throws "argument is of
  # length zero". That is not hypothetical: it is why `gsea_barplot()` errored
  # outright whenever nothing passed `padj_cutoff` -- the normal outcome at the
  # default 0.05. An explicit zero-row frame makes the empty plot answer
  # `nrow()`, `ncol()` and `$data` like every non-empty one.
  p <- ggplot(data.frame()) +
    labs(title = title, subtitle = subtitle) +
    theme_bulki()
  attr(p, "gs_source") <- data.frame()
  p
}

#' Attach the source table a plot was drawn from
#'
#' [gs_save()] writes this alongside the figure, so "a figure's source table is
#' its same-stem neighbour" is enforced by code rather than remembered.
#'
#' @param p A ggplot object.
#' @param data The data frame the plot was built from.
#' @return `p`, with a `gs_source` attribute.
#' @keywords internal
.gs_attach_source <- function(p, data) {
  attr(p, "gs_source") <- data
  p
}
