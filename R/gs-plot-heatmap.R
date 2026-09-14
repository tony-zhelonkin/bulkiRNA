#' Gene-set heatmap
#'
#' Pathways down the y axis, and along x either the contrasts/databases of a
#' [gs_result] or the samples of a [gs_matrix]. One generic covers what used to
#' be six heatmap functions built on \pkg{pheatmap}; the return value is a
#' ggplot like every other renderer, so [gs_save()] handles it.
#'
#' @param x A [gs_result] (tiles are the signed statistic) or a [gs_matrix]
#'   (tiles are the per-sample scores).
#' @param ... Passed to the method.
#' @return A ggplot object, carrying its source table in the `gs_source`
#'   attribute.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "SET_A", SET_B = "SET_B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_plot_heatmap(res)
#' @export
gs_plot_heatmap <- function(x, ...) {
  UseMethod("gs_plot_heatmap")
}

#' @export
gs_plot_heatmap.default <- function(x, ...) {
  stop("`x` must be a gs_result or a gs_matrix; got ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}

#' Heatmap of a `gs_result`
#'
#' @param x A [gs_result].
#' @param top_n Number of pathways to display, selected over the whole object.
#' @param sort_by Selection metric: `"padj"` (default), `"p_value"`, `"stat"`
#'   or `"stat_signed"`.
#' @param by What the x axis carries: `"contrast"` (default) or `"database"`.
#' @param highlight FDR threshold below which a tile is marked with `*`, or
#'   `NULL` for no marks.
#' @param limits Length-2 numeric fill limits. `NULL` (default) derives
#'   symmetric limits **from this figure's own data**, so the same statistic can
#'   take a different colour in two figures; pass an explicit `limits` when
#'   panels are meant to be compared. The old renderer used a fixed
#'   `c(-3.5, 3.5)`. Values outside are squished, not dropped.
#' @param palette Length-3 character vector -- low, mid, high fill colours.
#' @param wrap_width Soft character width for wrapping pathway labels.
#' @param strip_prefix Logical, passed to [format_pathway_name()].
#' @param title Plot title, or `NULL`.
#' @param database_labels Named character vector mapping `database` keys to
#'   display labels.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
gs_plot_heatmap.gs_result <- function(x,
                                      top_n = 20,
                                      sort_by = c("padj", "p_value", "stat",
                                                  "stat_signed"),
                                      by = c("contrast", "database"),
                                      highlight = 0.05,
                                      limits = NULL,
                                      palette = .gs_diverging_colours(),
                                      wrap_width = 40,
                                      strip_prefix = TRUE,
                                      title = NULL,
                                      database_labels = NULL,
                                      ...) {
  sort_by <- match.arg(sort_by)
  by <- match.arg(by)

  # `top_n` counts pathways, not rows: a pathway tested in three contrasts is
  # one row of the heatmap, not three.
  sel <- .gs_plot_frame(
    x, top_n = NULL, sort_by = sort_by, highlight = highlight,
    wrap_width = wrap_width, strip_prefix = strip_prefix,
    database_labels = database_labels
  )
  if (nrow(sel) == 0L) return(.gs_empty_plot(title))
  ids <- unique(sel$pathway_id)
  if (!is.null(top_n) && is.finite(top_n)) {
    ids <- utils::head(ids, as.integer(top_n))
  }

  df <- .gs_plot_frame(
    x, keep_ids = ids, sort_by = sort_by,
    highlight = highlight, wrap_width = wrap_width,
    strip_prefix = strip_prefix, database_labels = database_labels
  )
  df$column <- if (by == "database") df$database_label else df[[by]]
  if (by == "contrast") {
    # `gs_test()`'s `contrast` formal defaults to the literal placeholder string
    # "contrast", so the common single-contrast case -- including this
    # function's own roxygen example -- printed an axis tick reading the word
    # "contrast". Blank only the placeholder and missing values; a caller who
    # supplied a real contrast name keeps it.
    df$column[is.na(df$column) | df$column == "contrast"] <- ""
  }
  if (is.null(limits)) limits <- .gs_symmetric_limits(df$stat)
  df <- .gs_order_labels(df, by = "stat", decreasing = TRUE)

  p <- ggplot(df, aes(x = .data$column, y = .data$label)) +
    geom_tile(aes(fill = .data$stat), colour = "white", linewidth = 0.3)

  marks <- df[df$significant, , drop = FALSE]
  if (nrow(marks) > 0L) {
    p <- p + geom_text(data = marks, label = "*", vjust = 0.75, size = 5)
  }

  p <- p +
    .gs_fill_scale(gs_stat_label(x), limits, palette) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bulki() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())

  .gs_attach_source(p, df)
}

#' Heatmap of a `gs_matrix`
#'
#' @param x A [gs_matrix] of pathway scores.
#' @param top_n Number of pathways to display, chosen by score variance across
#'   samples. `NULL` shows all.
#' @param samples Character vector selecting and ordering the samples, or
#'   `NULL` for the matrix order.
#' @param group Name of a column in the matrix's `sample_data` to facet the
#'   columns by, or `NULL`.
#' @param limits Length-2 numeric fill limits. `NULL` (default) derives
#'   symmetric limits **from this figure's own data**, so the same statistic can
#'   take a different colour in two figures; pass an explicit `limits` when
#'   panels are meant to be compared. The old renderer used a fixed
#'   `c(-3.5, 3.5)`. Values outside are squished, not dropped.
#' @param palette Length-3 character vector -- low, mid, high fill colours.
#' @param wrap_width Soft character width for wrapping pathway labels.
#' @param strip_prefix Logical, passed to [format_pathway_name()].
#' @param title Plot title, or `NULL`.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
gs_plot_heatmap.gs_matrix <- function(x,
                                      top_n = 30,
                                      samples = NULL,
                                      group = NULL,
                                      limits = NULL,
                                      palette = .gs_diverging_colours(),
                                      wrap_width = 40,
                                      strip_prefix = TRUE,
                                      title = NULL,
                                      ...) {
  m <- x
  if (!is.null(samples)) {
    unknown <- setdiff(samples, colnames(m))
    if (length(unknown)) {
      stop("`samples` not found in the gs_matrix: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    m <- m[, samples, drop = FALSE]
  }
  if (!is.null(top_n) && is.finite(top_n) && nrow(m) > top_n) {
    spread <- apply(as.matrix(m), 1L, stats::var, na.rm = TRUE)
    keep <- utils::head(order(spread, decreasing = TRUE), as.integer(top_n))
    m <- m[sort(keep), , drop = FALSE]
  }
  if (nrow(m) == 0L || ncol(m) == 0L) return(.gs_empty_plot(title))

  names_map <- gs_pathway_names(x)
  ids <- rownames(m)
  pretty <- if (is.null(names_map)) ids else unname(names_map[ids])
  pretty[is.na(pretty)] <- ids[is.na(pretty)]

  df <- data.frame(
    pathway_id = rep(ids, times = ncol(m)),
    label = rep(
      .gs_wrap_label(format_pathway_name(pretty, strip_prefix = strip_prefix),
                     width = wrap_width),
      times = ncol(m)
    ),
    sample = rep(colnames(m), each = nrow(m)),
    score = as.numeric(as.matrix(m)),
    stringsAsFactors = FALSE
  )
  df$sample <- factor(df$sample, levels = colnames(m))

  if (!is.null(group)) {
    sd <- gs_sample_data(x)
    if (is.null(sd) || is.null(sd[[group]])) {
      stop("`group` column ", dQuote(group),
           " is not in the gs_matrix sample_data.", call. = FALSE)
    }
    df$.facet_col <- sd[[group]][match(as.character(df$sample), rownames(sd))]
  }

  if (is.null(limits)) limits <- .gs_symmetric_limits(df$score)
  df <- .gs_order_labels(df, by = "score", decreasing = TRUE)

  p <- ggplot(df, aes(x = .data$sample, y = .data$label)) +
    geom_tile(aes(fill = .data$score), colour = "white", linewidth = 0.3) +
    .gs_fill_scale(paste0(gs_score_type(x), " score"), limits, palette) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bulki() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())

  if (!is.null(group)) {
    p <- p + facet_grid(cols = vars(.data$.facet_col),
                        scales = "free_x", space = "free_x")
  }
  .gs_attach_source(p, df)
}
