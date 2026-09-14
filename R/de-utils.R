#' Theme used by the `de_*` renderers
#'
#' Resolves to the shared `theme_bulki()` look when that renderer-layer theme is
#' present in the package namespace, and otherwise falls back to a
#' self-contained minimal theme with the same visual contract (white
#' background, visible axis lines, no grid). The indirection exists because the
#' DE module and the theme live in different modules of the refactor.
#'
#' Two details are load-bearing, both of which this function got wrong when the
#' modules were first merged:
#'
#' * `base_size` is **forwarded**. The original called `theme_bulki()` with no
#'   arguments, so a caller asking for `base_size = 20` silently got 14.
#' * The floor is **skipped** (`floor = NULL`). The DE renderers' documented
#'   default is 12 pt and their geometry derives from it, so `theme_bulki()`'s
#'   deliberate 14 pt floor scaled every volcano and MD plot by 14/12 -- point
#'   stroke 1.5 -> 1.75, line widths to match. The floor still applies to the
#'   `gs_plot_*` layer it was designed for.
#'
#' @param base_size Base font size.
#' @return A `ggplot2` theme object.
#' @keywords internal
.de_theme <- function(base_size = 12) {
  ns <- asNamespace("bulkiRNA")
  if (exists(".theme_bulki", envir = ns, inherits = FALSE)) {
    fn <- get(".theme_bulki", envir = ns)
    return(fn(base_size = base_size, floor = NULL))
  }
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", colour = "transparent"),
      plot.background  = element_rect(fill = "white", colour = "transparent"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line   = element_line(colour = "black", linewidth = 0.5),
      axis.ticks  = element_line(colour = "black", linewidth = 0.5),
      panel.border = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text.x  = element_text(margin = margin(t = 5)),
      axis.text.y  = element_text(margin = margin(r = 5)),
      plot.title   = element_text(hjust = 0.5, margin = margin(b = 10)),
      plot.margin  = margin(10, 10, 10, 10)
    )
}

#' Darken a colour
#'
#' Multiplies each RGB channel by `factor`. Used to derive label text colours
#' that stay legible against the point colour they belong to.
#'
#' @param hex Character vector of colours.
#' @param factor Numeric multiplier in `[0, 1]`.
#' @return A character vector of hex colours, same length as `hex`.
#' @keywords internal
.de_shade <- function(hex, factor = 0.6) {
  vapply(hex, function(h) {
    rgb <- grDevices::col2rgb(h) / 255 * factor
    rgb <- pmax(pmin(rgb, 1), 0)
    grDevices::rgb(rgb[1], rgb[2], rgb[3])
  }, character(1))
}

