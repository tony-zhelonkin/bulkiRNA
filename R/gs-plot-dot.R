#' Gene-set dotplot
#'
#' One dotplot with arguments where the toolkit used to have eight functions.
#' Pathways run down the y axis; the statistic (or the leading-edge gene ratio)
#' runs along x; dot area is `-log10(padj)`; dot fill is the signed statistic on
#' a colorblind-safe diverging scale.
#'
#' **Selection and highlighting are separate**, which is the behaviour the old
#' `gsea_dotplot()` established and worth keeping: `top_n` and `sort_by` decide
#' *what is shown*, `highlight` decides *what gets a black outline*. A pathway
#' can be displayed without being significant.
#'
#' The axis label for the statistic comes from the object's `stat_type`, so a
#' GSVA t-statistic is never mislabelled `NES`. Facet strips and legends show
#' the database *label*, never the snake_case `database` key.
#'
#' @param x A [gs_result].
#' @param top_n Number of pathways to display, per selection group. `NULL` shows
#'   all.
#' @param aes_x What the x axis carries: `"stat"` (default) or `"gene_ratio"`
#'   (leading-edge size over set size; requires a `leading_edge` or `overlap`
#'   column).
#' @param sort_by Selection metric: `"padj"` (default), `"p_value"`, `"stat"`
#'   (largest absolute value) or `"stat_signed"` (most positive).
#' @param direction Restrict to `"both"` (default), `"up"` or `"down"`.
#' @param facet Row facets: `"none"` (default), `"direction"` -- the old
#'   `gsea_dotplot_facet()` -- `"database"` or `"contrast"`. `top_n` applies
#'   within each facet.
#' @param compare Column facets for side-by-side comparison: `NULL` (default),
#'   `"contrast"` or `"database"`. Pathways are selected once over the whole
#'   object and every panel then shows the same set, so the grid is not ragged.
#' @param highlight FDR threshold for the black significance outline, or `NULL`
#'   for no outline.
#' @param label Logical. Annotate each dot with its FDR, placed with
#'   \pkg{ggrepel}.
#' @param size_range Length-2 numeric, the dot size range.
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
#'   display labels. Only needed when the object carries neither a
#'   `database_label` column nor attribute.
#' @return A ggplot object, carrying its source table in the `gs_source`
#'   attribute for [gs_save()].
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "SET_A", SET_B = "SET_B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' gs_plot_dot(res, top_n = 5)
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales label_scientific
#' @export
gs_plot_dot <- function(x,
                        top_n = 10,
                        aes_x = c("stat", "gene_ratio"),
                        sort_by = c("padj", "p_value", "stat", "stat_signed"),
                        direction = c("both", "up", "down"),
                        facet = c("none", "direction", "database", "contrast"),
                        compare = NULL,
                        highlight = 0.05,
                        label = FALSE,
                        size_range = c(2, 10),
                        limits = NULL,
                        palette = .gs_diverging_colours(),
                        wrap_width = 50,
                        strip_prefix = TRUE,
                        title = NULL,
                        database_labels = NULL) {
  .gs_plot_check_result(x)
  aes_x <- match.arg(aes_x)
  sort_by <- match.arg(sort_by)
  direction <- match.arg(direction)
  facet <- match.arg(facet)
  compare <- .gs_check_compare(compare)

  group_by <- if (facet == "none") NULL else facet

  df <- .gs_plot_frame(
    x, top_n = top_n, sort_by = sort_by, direction = direction,
    group_by = group_by, highlight = highlight, wrap_width = wrap_width,
    strip_prefix = strip_prefix, database_labels = database_labels
  )
  if (!is.null(compare) && nrow(df) > 0L) {
    df <- .gs_plot_frame(
      x, keep_ids = unique(df$pathway_id), sort_by = sort_by,
      direction = direction, highlight = highlight, wrap_width = wrap_width,
      strip_prefix = strip_prefix, database_labels = database_labels
    )
  }
  if (nrow(df) == 0L) return(.gs_empty_plot(title))

  if (aes_x == "gene_ratio" && is.null(df$gene_ratio)) {
    stop("`aes_x = \"gene_ratio\"` needs a `leading_edge` or `overlap` ",
         "column; this gs_result has neither. Use `aes_x = \"stat\"`.",
         call. = FALSE)
  }
  x_lab <- if (aes_x == "gene_ratio") "Gene ratio" else gs_stat_label(x)
  fill_lab <- gs_stat_label(x)
  if (is.null(limits)) limits <- .gs_symmetric_limits(df$stat)

  df <- .gs_order_labels(df, by = aes_x, decreasing = TRUE)
  df <- .gs_facet_columns(df, facet, compare)

  p <- ggplot(df, aes(x = .data[[aes_x]], y = .data$label)) +
    geom_point(
      aes(size = .data$neg_log_padj, fill = .data$stat),
      shape = 21, stroke = 0, colour = "transparent"
    )

  sig <- df[df$significant, , drop = FALSE]
  if (nrow(sig) > 0L) {
    p <- p + geom_point(
      data = sig,
      aes(size = .data$neg_log_padj),
      shape = 21, colour = "black", fill = "transparent", stroke = 1.2
    )
  }

  if (isTRUE(label)) {
    df$point_label <- scales::label_scientific(digits = 2)(df$padj)
    p <- p + ggrepel::geom_text_repel(
      data = df, aes(label = .data$point_label),
      size = 3, min.segment.length = 0, max.overlaps = Inf,
      segment.colour = "grey60"
    )
  }

  p <- p +
    .gs_fill_scale(fill_lab, limits, palette) +
    scale_size_continuous(
      name = expression(-log[10] ~ FDR), range = size_range
    ) +
    guides(size = guide_legend(
      override.aes = list(shape = 16, fill = "black", colour = "black")
    )) +
    labs(title = title, x = x_lab, y = NULL) +
    theme_bulki() +
    theme(legend.position = "right")

  p <- .gs_add_facets(p, facet, compare)
  .gs_attach_source(p, df)
}

#' Validate the `compare` argument
#'
#' @param compare `NULL`, `"contrast"` or `"database"`.
#' @return `compare`, unchanged.
#' @keywords internal
.gs_check_compare <- function(compare) {
  if (is.null(compare)) return(NULL)
  allowed <- c("contrast", "database")
  if (!is.character(compare) || length(compare) != 1L ||
        !compare %in% allowed) {
    stop("`compare` must be NULL, \"contrast\" or \"database\"; got ",
         dQuote(paste(compare, collapse = ", ")), ".", call. = FALSE)
  }
  compare
}

#' Add the `.facet_row` / `.facet_col` columns a facet spec needs
#'
#' `database` facets read `database_label`, never the key.
#'
#' @param df A plotting frame.
#' @param facet Row facet name, or `"none"`.
#' @param compare Column facet name, or `NULL`.
#' @return `df`, with the facet columns added where required.
#' @keywords internal
.gs_facet_columns <- function(df, facet = "none", compare = NULL) {
  pick <- function(nm) {
    if (nm == "direction") .gs_direction_factor(df$direction)
    else if (nm == "database") df$database_label
    else df[[nm]]
  }
  if (!identical(facet, "none")) df$.facet_row <- pick(facet)
  if (!is.null(compare)) df$.facet_col <- pick(compare)
  df
}

#' Attach the facet layer matching a facet spec
#'
#' @param p A ggplot object.
#' @param facet Row facet name, or `"none"`.
#' @param compare Column facet name, or `NULL`.
#' @return `p`, faceted where required.
#' @keywords internal
.gs_add_facets <- function(p, facet = "none", compare = NULL) {
  has_row <- !identical(facet, "none")
  has_col <- !is.null(compare)
  if (!has_row && !has_col) return(p)
  p + facet_grid(
    rows = if (has_row) vars(.data$.facet_row) else NULL,
    cols = if (has_col) vars(.data$.facet_col) else NULL,
    scales = if (has_row) "free_y" else "fixed",
    space = if (has_row) "free_y" else "fixed"
  ) + theme(panel.spacing.y = unit(1, "lines"))
}
