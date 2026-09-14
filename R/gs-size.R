#' Scale a figure's text to its canvas
#'
#' Point sizes are absolute; canvas sizes are not. A theme built at 14 pt reads
#' well on a 7 x 5 in figure and looks oversized on a 9 x 10 in one, which is
#' how a renderer that never sees the canvas ends up producing labels that
#' crowd out the panel. This applies one scale factor,
#' `sqrt((width * height) / (7 * 5))`, to every text element, so the *relative*
#' typography a renderer chose survives and only the absolute size moves.
#'
#' It works on a plain `ggplot` and on a `patchwork`, and it does not care who
#' built the object. Use it on figures this package did not draw.
#'
#' The scaled sizes are absolute, so they supersede any base-size floor,
#' including [theme_bulki()]'s 14 pt.
#'
#' @param plot A `ggplot` or `patchwork`.
#' @param width,height Canvas size in inches -- the values you will pass to
#'   [ggplot2::ggsave()].
#' @param base_font_size Numeric. Reference point size at the 7 x 5 in
#'   reference canvas.
#' @return The plot with a text-size theme applied. For a `patchwork` the
#'   theme is applied to every panel with patchwork's `&`.
#' @seealso [gs_save()], [gs_plot_size()]
#' @examples
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
#' gs_scale_fonts(p, width = 12, height = 9)
#' @export
gs_scale_fonts <- function(plot, width, height, base_font_size = 10) {
  if (!inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    stop("`plot` must be a ggplot or patchwork object; got ",
         paste(class(plot), collapse = "/"), ".", call. = FALSE)
  }
  for (arg in c("width", "height", "base_font_size")) {
    value <- get(arg)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
          !is.finite(value) || value <= 0) {
      stop("`", arg, "` must be one finite positive number.", call. = FALSE)
    }
  }

  th <- .gs_font_theme(width, height, base_font_size)
  # `&` distributes a theme across a patchwork's panels; `+` would attach it to
  # the composition only and leave the panels at their original sizes.
  if (inherits(plot, "patchwork")) {
    return(plot & th)
  }
  plot + th
}

#' The text-size theme for a canvas
#'
#' Ratios ported from the pre-package `save_gsea_plot()`, which is what made
#' the original figures readable across databases of very different density.
#'
#' @param width,height Canvas size in inches.
#' @param base_font_size Reference point size.
#' @return A ggplot2 theme.
#' @keywords internal
.gs_font_theme <- function(width, height, base_font_size) {
  scale_fac <- sqrt((width * height) / (7 * 5))
  size <- base_font_size * scale_fac
  theme(
    text = element_text(size = size),
    axis.title = element_text(size = size),
    axis.text = element_text(size = size * 0.9),
    plot.title = element_text(size = size * 1.2),
    plot.subtitle = element_text(size = size * 0.9),
    legend.title = element_text(size = size * 0.9),
    legend.text = element_text(size = size * 0.8),
    strip.text = element_text(size = size)
  )
}

#' Suggested canvas size for a figure
#'
#' A *suggestion*, not a policy. Panel density depends on how many sets a
#' database returns and how long their names are, which the package cannot
#' predict -- so this returns numbers you are expected to override, and every
#' renderer works at whatever canvas you actually pass.
#'
#' The table is ported from the pre-package `get_db_plot_params()`, including
#' its one structural rule: a running-sum figure gets 20% more height than the
#' other plot types for the same database, because it stacks three panels.
#'
#' @param type One of `"bar"`, `"dot"`, `"facet"`, `"running"`, `"heatmap"`.
#' @param database Optional database name, matched case-insensitively against
#'   the known collections (`hallmark`, `gobp`, `gomf`, `gocc`, `kegg`,
#'   `reactome`, `biocarta`, `wikipathways`). Anything unrecognised, including
#'   `NULL`, gets the default.
#' @return A list with `width`, `height` (inches) and `base_font_size`.
#' @seealso [gs_save()], [gs_scale_fonts()]
#' @examples
#' gs_plot_size("running", "GO_BP")
#' gs_plot_size("bar", "Hallmark")
#'
#' # Override freely -- these are starting points.
#' size <- gs_plot_size("dot", "Reactome")
#' size$height <- 12
#' @export
gs_plot_size <- function(type = c("bar", "dot", "facet", "running", "heatmap"),
                         database = NULL) {
  type <- match.arg(type)
  base <- .gs_size_table(database)
  # Height multipliers, one per plot type: a running-sum stacks three panels
  # and a faceted dotplot doubles its rows.
  mult <- switch(type, running = 1.2, facet = 1.4, 1.0)
  list(
    width = base$width,
    height = base$height * mult,
    base_font_size = base$base_font_size
  )
}

#' Per-database canvas defaults
#'
#' @param database Optional database name.
#' @return A list with `width`, `height` and `base_font_size`.
#' @keywords internal
.gs_size_table <- function(database) {
  default <- list(width = 7, height = 5, base_font_size = 10)
  if (is.null(database) || !is.character(database) ||
        length(database) != 1L || is.na(database)) {
    return(default)
  }
  # Matched on a normalised key, so "GO_BP", "gobp" and "GO BP" agree.
  key <- gsub("[^a-z0-9]", "", tolower(database))
  table <- list(
    hallmark = list(width = 7, height = 5, base_font_size = 10),
    gobp = list(width = 8, height = 8, base_font_size = 9),
    gomf = list(width = 8, height = 7, base_font_size = 9),
    gocc = list(width = 7, height = 6, base_font_size = 9),
    kegg = list(width = 8, height = 7, base_font_size = 9),
    reactome = list(width = 9, height = 8, base_font_size = 8),
    biocarta = list(width = 7, height = 6, base_font_size = 9),
    wikipathways = list(width = 8, height = 7, base_font_size = 9)
  )
  table[[key]] %||% default
}
