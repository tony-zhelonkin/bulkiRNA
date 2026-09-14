#' The bulkiRNA plot theme
#'
#' A clean, publication-oriented theme: white background, visible axis lines
#' and ticks, no panel border, and -- by default -- no grid. Every `gs_plot_*`
#' and `de_*` renderer ends with it, so figure styling is one function rather
#' than a `theme()` call an agent might forget.
#'
#' The **14 pt base-size floor** lives here on purpose. A figure destined for a
#' figure panel is unreadable below it, so `base_size` values under 14 are
#' raised to 14 rather than honoured. Pass `base_size = 14` explicitly if you
#' want the floor to be visible in your code.
#'
#' @param base_size Base font size in points. Values below `14` are raised to
#'   `14`.
#' @param base_family Base font family. Default `""` (the device default).
#' @param grid Logical. `FALSE` (default) removes all grid lines; `TRUE` keeps
#'   a light major grid on both axes.
#' @return A ggplot2 theme object.
#' @examples
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'   ggplot2::geom_point() +
#'   theme_bulki()
#' @export
theme_bulki <- function(base_size = 14, base_family = "", grid = FALSE) {
  .theme_bulki(base_size = base_size, base_family = base_family,
               grid = grid, floor = 14)
}

#' `theme_bulki()` with the base-size floor as a parameter
#'
#' The floor is deliberate for the public theme, but the deprecated
#' `custom_minimal_theme_with_grid()` shim has to reproduce the old theme
#' exactly, and the old theme had no floor -- its documented default was 12 pt
#' and it rendered at 12 pt. Routing both through this helper lets the public
#' function keep its floor while the shim passes `floor = NULL`, instead of
#' either silently enlarging every legacy figure or clamping in two places.
#'
#' @inheritParams theme_bulki
#' @param floor Numeric(1) smallest honoured `base_size`, or `NULL` for none.
#' @return A ggplot2 theme object.
#' @keywords internal
#' @name dot-theme_bulki
.theme_bulki <- function(base_size = 14, base_family = "", grid = FALSE,
                         floor = 14) {
  if (!is.numeric(base_size) || length(base_size) != 1L || is.na(base_size)) {
    stop("`base_size` must be a single number.", call. = FALSE)
  }
  if (!is.logical(grid) || length(grid) != 1L || is.na(grid)) {
    stop("`grid` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(floor)) base_size <- max(base_size, floor)

  grid_major <- if (grid) {
    element_line(colour = "grey92", linewidth = 0.3)
  } else {
    element_blank()
  }

  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "white", colour = "transparent"),
      plot.background = element_rect(fill = "white", colour = "transparent"),
      panel.grid.major = grid_major,
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )
}
