#' Gene-set barplot
#'
#' Horizontal bars of the signed statistic, one per pathway, filled by the same
#' diverging scale the dotplot uses. The axis label comes from the object's
#' `stat_type`, so the bars are labelled `NES`, `t statistic` or
#' `log2 fold enrichment` as appropriate.
#'
#' Unlike the old `gsea_barplot()`, non-significant pathways are **not silently
#' dropped**: `top_n` and `sort_by` choose what is shown and `highlight` outlines
#' the significant ones. Pass `padj_max` if you do want a hard significance
#' filter.
#'
#' @param x A [gs_result].
#' @param top_n Number of pathways to display, per selection group.
#' @param sort_by Selection metric: `"stat"` (default, largest absolute value),
#'   `"padj"`, `"p_value"` or `"stat_signed"`.
#' @param direction Restrict to `"both"` (default), `"up"` or `"down"`.
#' @param facet Row facets: `"none"` (default), `"direction"`, `"database"` or
#'   `"contrast"`.
#' @param padj_max Optional hard FDR filter applied before selection.
#' @param highlight FDR threshold for the black bar outline, or `NULL`.
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
#' gs_plot_bar(res, top_n = 5)
#' @export
gs_plot_bar <- function(x,
                        top_n = 20,
                        sort_by = c("stat", "padj", "p_value", "stat_signed"),
                        direction = c("both", "up", "down"),
                        facet = c("none", "direction", "database", "contrast"),
                        padj_max = NULL,
                        highlight = 0.05,
                        limits = NULL,
                        palette = .gs_diverging_colours(),
                        wrap_width = 50,
                        strip_prefix = TRUE,
                        title = NULL,
                        database_labels = NULL) {
  .gs_plot_check_result(x)
  sort_by <- match.arg(sort_by)
  direction <- match.arg(direction)
  facet <- match.arg(facet)

  if (!is.null(padj_max)) {
    if (!is.numeric(padj_max) || length(padj_max) != 1L) {
      stop("`padj_max` must be a single number, or NULL.", call. = FALSE)
    }
    x <- x[!is.na(x[["padj"]]) & x[["padj"]] < padj_max, , drop = FALSE]
  }

  df <- .gs_plot_frame(
    x, top_n = top_n, sort_by = sort_by, direction = direction,
    group_by = if (facet == "none") NULL else facet,
    highlight = highlight, wrap_width = wrap_width,
    strip_prefix = strip_prefix, database_labels = database_labels
  )
  if (nrow(df) == 0L) return(.gs_empty_plot(title))

  if (is.null(limits)) limits <- .gs_symmetric_limits(df$stat)
  df <- .gs_order_labels(df, by = "stat", decreasing = TRUE)
  df <- .gs_facet_columns(df, facet, NULL)
  df$outline <- ifelse(df$significant, "black", "transparent")

  p <- ggplot(df, aes(x = .data$label, y = .data$stat)) +
    geom_col(aes(fill = .data$stat, colour = .data$outline),
             linewidth = 0.4) +
    scale_colour_identity() +
    .gs_fill_scale(gs_stat_label(x), limits, palette) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    coord_flip() +
    labs(title = title, x = NULL, y = gs_stat_label(x)) +
    theme_bulki() +
    theme(legend.position = "none")

  p <- .gs_add_facets(p, facet, NULL)
  .gs_attach_source(p, df)
}
