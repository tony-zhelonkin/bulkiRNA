#' GSEA running-sum (enrichment curve) plot
#'
#' The canonical three-panel GSEA figure, built on
#' [fgsea::plotEnrichmentData()]:
#'
#' 1. the running enrichment score (ES) curve,
#' 2. one lane of gene-hit ticks per pathway,
#' 3. the ranked metric itself.
#'
#' # Why this returns a patchwork
#'
#' Because the three panels have unrelated y units, and a figure whose panels
#' share one y scale cannot be post-styled. An earlier version drew all three
#' as rows of a single [ggplot2::facet_grid()], rescaling each panel's data
#' into a synthetic window so that `space = "free_y"` would size the panels in
#' the requested ratio. That worked in isolation and broke on contact: any
#' downstream `coord_cartesian(ylim = )`, `scale_y_continuous()` or patchwork
#' `&` forced every panel to a common range, which flattened the height ratio
#' to 1:1:1, crushed each curve into a fraction of its panel and pushed the
#' inverse-mapped axis labels outside the visible window. Real panels have real
#' y scales, so a caller can restyle one without destroying the others.
#'
#' The return value therefore has class `patchwork`, and carries a
#' `grs_restyle` attribute -- a function that rebuilds the composition with new
#' layout arguments, reusing the already-computed curve data. It exists so a
#' project theme can be applied *through* the composer rather than on top of a
#' finished figure.
#'
#' # Why colours cannot be permuted
#'
#' The colour aesthetic is mapped to **`pathway_id`**, never to a display
#' label, and the scale is an explicit
#' `scale_colour_manual(values = <named vector>, breaks = names(values))` keyed
#' by that same id. Legend *text* is applied through the scale's `labels`, so
#' renaming a pathway can never move its colour.
#'
#' Colours are stable *within* a figure, not across figures: with no `palette`
#' the default is assigned in `sort(ids)` order, which is deterministic for a
#' given set of ids but still changes when the plotted set changes. For one
#' colour per pathway across a whole project, freeze a lookup once with
#' [gs_palette()] and pass it to every call.
#'
#' @param x A [gs_result][gs_result-class] (the usual case), a `gs_db`, or a
#'   named list of character vectors. A `gs_result` supplies pathway ids,
#'   display names and the `stat` used to pick the default pathways; it does
#'   **not** carry gene-set membership, so `db` is required with it.
#' @param ranks Named numeric vector of gene-level statistics, decreasing, as
#'   returned by [gs_ranks()]. If `NULL`, `attr(x, "ranks")` is used -- the hook
#'   a deprecation shim can fill in.
#' @param db A `gs_db` or named list of character vectors giving set
#'   membership. If `NULL`, `attr(x, "gene_sets")` is used, or `x` itself when
#'   `x` is already a `gs_db` / named list.
#' @param pathways Character vector of pathway ids, or (for a `gs_result`)
#'   integer row indices. `NULL` picks the top `top_n` pathways by `abs(stat)`
#'   for a `gs_result`, or the first `top_n` sets otherwise. The order given is
#'   the order legend keys are assigned in.
#' @param top_n Integer. How many pathways to pick when `pathways` is `NULL`.
#' @param labels Optional character vector of legend labels, named by pathway
#'   id (a partial map degrades gracefully -- unnamed ids keep their
#'   `pathway_name`), or unnamed and zipped to `pathways` in the given order.
#' @param palette Optional colours, named by pathway id or unnamed and zipped
#'   to `pathways` in the given order (recycled). `NULL` uses [gs_palette()].
#' @param panel_heights Length-3 numeric, the ES : ticks : metric height ratio,
#'   passed to [patchwork::wrap_plots()] as `heights`.
#' @param es_ylim Length-2 numeric y limit for the **ES panel only**, applied
#'   with [ggplot2::coord_cartesian()] so nothing is dropped. Defaults to
#'   `c(-1, 1)`, the mathematical range of a running enrichment score, so that
#'   curves from different contrasts and different databases are directly
#'   comparable without the reader checking each axis. `NULL` lets each figure
#'   pick its own range. The tick and metric panels are never clamped.
#' @param gsea_param Numeric exponent passed to
#'   [fgsea::plotEnrichmentData()]'s `gseaParam`.
#' @param linewidth Numeric. Weight of the ES curves; `1.6` by default, chosen
#'   so the figure still reads on a projected slide. Gene ticks are drawn at
#'   `0.45 *` this, thinner on purpose: at curve weight a dense set merges into
#'   a solid block. Lower it to about `1.1` for a print-only figure.
#' @param metric_label Character. Axis label for the ranked-metric panel; name
#'   the statistic you ranked by (e.g. `"t statistic"`, `"log2 FC"`).
#' @param title Optional plot title, drawn over the ES panel.
#' @param base_size Base font size. Note [theme_bulki()] applies its own 14 pt
#'   floor to the default foundation theme, and
#'   `gs_save(rescale_fonts = TRUE)` sets absolute sizes that supersede both.
#' @param max_name_length Integer. Legend labels longer than this are *wrapped*
#'   onto several lines, never truncated.
#' @param legend_position One of `"right"` (default), `"inside"`, `"bottom"`
#'   or `"none"`. The default keeps the legend out of the panel, so a long
#'   pathway name cannot sit over the curves; it costs figure width, which
#'   `"inside"` does not.
#' @param legend_pos Length-2 numeric, the inside-legend position in npc units,
#'   used only when `legend_position = "inside"`.
#' @param xticks `"bottom"` (default) draws x-axis text and ticks on the
#'   bottom panel only; `"all"` draws them on every panel.
#' @param rug_ylabels Logical. Show the tick panel's lane indices on its y
#'   axis. `FALSE` (default) hides them, since the index carries no meaning.
#' @param base_theme Optional complete [ggplot2::theme()] used as the
#'   foundation for every panel; the panel chrome is re-applied on top so a
#'   project theme cannot clobber the layout. `NULL` uses [theme_bulki()].
#' @return A `patchwork` of three panels, carrying a `grs_restyle` attribute:
#'   `function(es_ylim, legend_position, xticks, rug_ylabels, panel_heights,
#'   base_theme, base_size, ...)` returning a freshly composed figure.
#' @seealso [gs_ranks()], [gs_test()], [gs_leading_edge()], [gs_palette()]
#' @importFrom rlang .data
#' @importFrom fgsea plotEnrichmentData
#' @examples
#' sets <- list(SET_A = c("A", "B", "C"), SET_B = c("D", "E", "F"))
#' ranks <- stats::setNames(seq(3, -3, length.out = 8), LETTERS[1:8])
#' gs_plot_running(sets, ranks = ranks)
#' @export
gs_plot_running <- function(x,
                            ranks = NULL,
                            db = NULL,
                            pathways = NULL,
                            top_n = 5L,
                            labels = NULL,
                            palette = NULL,
                            panel_heights = c(2.4, 0.7, 0.9),
                            es_ylim = c(-1, 1),
                            gsea_param = 1,
                            linewidth = 1.6,
                            metric_label = "Ranked metric",
                            title = NULL,
                            base_size = 14,
                            max_name_length = 40,
                            legend_position = c("right", "inside", "bottom",
                                                "none"),
                            legend_pos = c(0.98, 0.98),
                            xticks = c("bottom", "all"),
                            rug_ylabels = FALSE,
                            base_theme = NULL) {
  legend_position <- match.arg(legend_position)
  xticks <- match.arg(xticks)
  .grs_check_heights(panel_heights)
  .grs_check_ylim(es_ylim)
  if (!is.numeric(linewidth) || length(linewidth) != 1L ||
        is.na(linewidth) || !is.finite(linewidth) || linewidth <= 0) {
    stop("`linewidth` must be one finite positive number.", call. = FALSE)
  }

  ranks <- .grs_ranks(x, ranks)
  sets <- .grs_sets(x, db)
  ids <- .grs_select(x, sets, pathways, top_n)
  set_labels <- .grs_labels(x, ids, labels, max_name_length)
  pal <- .grs_palette(palette, ids)

  axis_labels <- .grs_panel_levels(metric_label)
  curves <- .grs_curves(sets[ids], ranks, gsea_param)
  frames <- .grs_frames(curves, ids)

  # One x range for all three panels, pinned rather than inferred.
  # plotEnrichmentData()'s curve carries sentinel ranks at 0 and n + 1 while
  # the ticks and the metric run 1..n, so per-panel limits differ by a rank at
  # each end -- enough to visibly offset a tick from the curve feature it
  # marks. A facet shared one scale by construction; real panels have to be
  # told.
  xlim <- range(c(frames$es$rank, frames$ticks$rank, frames$stats$rank))

  # The composer is a closure over the computed curve data, so a restyle
  # re-lays-out the figure without re-running fgsea.
  compose <- function(es_ylim = c(-1, 1),
                      legend_position = "right",
                      xticks = "bottom",
                      rug_ylabels = FALSE,
                      panel_heights = c(2.4, 0.7, 0.9),
                      base_theme = NULL,
                      base_size = 14,
                      linewidth = 1.6,
                      ...) {
    legend_position <- match.arg(
      legend_position, c("inside", "right", "bottom", "none")
    )
    xticks <- match.arg(xticks, c("bottom", "all"))
    .grs_check_heights(panel_heights)
    .grs_check_ylim(es_ylim)
    show_x <- identical(xticks, "all")
    base <- base_theme %||% .grs_default_theme(base_size)

    es <- .grs_panel_es(
      frames$es, pal, set_labels, es_ylim, xlim, axis_labels[["es"]], title,
      base, show_x, legend_position, legend_pos, linewidth
    )
    ticks <- .grs_panel_ticks(
      frames$ticks, pal, length(ids), xlim, axis_labels[["ticks"]], base,
      show_x, rug_ylabels, linewidth
    )
    metric <- .grs_panel_metric(
      frames$stats, xlim, axis_labels[["stats"]], base, show_x = TRUE
    )
    .grs_compose(list(es, ticks, metric), panel_heights, legend_position)
  }

  out <- compose(
    es_ylim = es_ylim, legend_position = legend_position, xticks = xticks,
    rug_ylabels = rug_ylabels, panel_heights = panel_heights,
    base_theme = base_theme, base_size = base_size, linewidth = linewidth
  )
  # A styling contract that is only prose gets dropped in a refactor and the
  # consumer degrades silently. This attribute is asserted by the test suite.
  attr(out, "grs_restyle") <- compose
  out
}

# ---- internals (prefix `.grs_`) ---------------------------------------------

#' Default ES : ticks : metric height ratio
#'
#' The single source of truth for panel proportions.
#'
#' @return A length-3 numeric vector.
#' @keywords internal
.GRS_PANEL_HEIGHTS <- c(2.4, 0.7, 0.9)

#' Validate `panel_heights`
#'
#' @param panel_heights The renderer's `panel_heights`.
#' @return `TRUE`, invisibly.
#' @keywords internal
.grs_check_heights <- function(panel_heights) {
  if (!is.numeric(panel_heights) || length(panel_heights) != 3L ||
        any(!is.finite(panel_heights)) || any(panel_heights <= 0)) {
    stop("`panel_heights` must be three positive finite numbers ",
         "(ES : ticks : metric).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate `es_ylim`
#'
#' @param es_ylim The renderer's `es_ylim`, possibly `NULL`.
#' @return `TRUE`, invisibly.
#' @keywords internal
.grs_check_ylim <- function(es_ylim) {
  if (is.null(es_ylim)) return(invisible(TRUE))
  if (!is.numeric(es_ylim) || length(es_ylim) != 2L ||
        any(!is.finite(es_ylim)) || es_ylim[[1L]] >= es_ylim[[2L]]) {
    stop("`es_ylim` must be two finite increasing numbers, or NULL.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Reshape [fgsea::plotEnrichmentData()] output into per-panel frames
#'
#' Tick lanes span their whole lane: lane `i` covers `[n - i, n - i + 1]`, so
#' with the tick panel's limits pinned to `c(0, n)` every tick is `1/n` of that
#' panel's height regardless of how many pathways are drawn. The previous
#' version used `[n - i + 0.2, n - i + 0.8]` inside a rescaled window, which
#' left each tick at roughly 2% of the panel and made them unreadable.
#'
#' @param curves Named list of `plotEnrichmentData()` outputs.
#' @param ids Character vector of pathway ids, in plotting order.
#' @return A list of three data frames: `es`, `ticks`, `stats`.
#' @keywords internal
.grs_frames <- function(curves, ids) {
  n <- length(ids)
  es <- do.call(rbind, lapply(ids, function(id) {
    d <- curves[[id]]$curve
    data.frame(rank = d$rank, y = d$ES, pathway_id = id,
               stringsAsFactors = FALSE)
  }))
  ticks <- do.call(rbind, lapply(seq_along(ids), function(i) {
    d <- curves[[ids[i]]]$ticks
    data.frame(rank = d$rank, ymin = n - i, ymax = n - i + 1,
               pathway_id = ids[i], stringsAsFactors = FALSE)
  }))
  st <- curves[[ids[1L]]]$stats
  list(
    es = es,
    ticks = ticks,
    stats = data.frame(rank = st$rank, y = st$stat, stringsAsFactors = FALSE)
  )
}

#' The enrichment-score panel
#'
#' Owns the single colour scale and the legend. `es_ylim` is applied with
#' [ggplot2::coord_cartesian()], which zooms without dropping data.
#'
#' @param df The `es` frame from [.grs_frames()].
#' @param pal Named colour vector.
#' @param set_labels Named legend labels.
#' @param es_ylim Optional length-2 numeric, or `NULL`.
#' @param xlim Length-2 numeric x range shared by all three panels.
#' @param y_lab,title Axis label and plot title.
#' @param base Foundation theme.
#' @param show_x Logical. Draw x-axis text and ticks.
#' @param legend_position,legend_pos Legend placement.
#' @return A `ggplot`.
#' @keywords internal
.grs_panel_es <- function(df, pal, set_labels, es_ylim, xlim, y_lab, title,
                          base, show_x, legend_position, legend_pos,
                          linewidth = 1.6) {
  p <- ggplot(df, aes(x = .data$rank, y = .data$y,
                      colour = .data$pathway_id)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = linewidth) +
    scale_colour_manual(
      values = pal,
      breaks = names(pal),
      labels = unname(set_labels[names(pal)]),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = y_lab, title = title) +
    coord_cartesian(xlim = xlim, ylim = es_ylim)
  p <- p + .grs_panel_theme(base) + .grs_legend_theme(legend_position,
                                                      legend_pos)
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' The gene-tick panel
#'
#' @param df The `ticks` frame from [.grs_frames()].
#' @param pal Named colour vector.
#' @param n Number of pathways, which is the panel's y extent.
#' @param xlim Length-2 numeric x range shared by all three panels.
#' @param y_lab Axis label.
#' @param base Foundation theme.
#' @param show_x Logical. Draw x-axis text and ticks.
#' @param rug_ylabels Logical. Show the lane indices.
#' @return A `ggplot`.
#' @keywords internal
.grs_panel_ticks <- function(df, pal, n, xlim, y_lab, base, show_x,
                             rug_ylabels, linewidth = 1.6) {
  # One pathway needs no colour to tell it apart, and black reads better than
  # an arbitrary hue -- the pre-package renderer did the same.
  if (n == 1L) pal <- stats::setNames(rep("black", length(pal)), names(pal))
  p <- ggplot(df) +
    geom_segment(
      aes(x = .data$rank, xend = .data$rank,
          y = .data$ymin, yend = .data$ymax,
          colour = .data$pathway_id),
      # Ticks track the curve weight so the two read as one figure, but stay
      # thinner: at curve weight they merge into a solid block.
      linewidth = linewidth * 0.45
    ) +
    scale_colour_manual(values = pal, breaks = names(pal), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      limits = c(0, n),
      breaks = seq_len(n) - 0.5,
      labels = as.character(seq_len(n)),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = xlim) +
    labs(x = NULL, y = y_lab) +
    .grs_panel_theme(base) +
    theme(legend.position = "none")
  if (!isTRUE(rug_ylabels)) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' The ranked-metric panel
#'
#' @param df The `stats` frame from [.grs_frames()].
#' @param xlim Length-2 numeric x range shared by all three panels.
#' @param y_lab Axis label.
#' @param base Foundation theme.
#' @param show_x Logical. Draw x-axis text and ticks.
#' @return A `ggplot`.
#' @keywords internal
.grs_panel_metric <- function(df, xlim, y_lab, base, show_x) {
  p <- ggplot(df, aes(x = .data$rank)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = 0, ymax = .data$y),
                fill = "grey70", colour = "transparent") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(n.breaks = 4) +
    coord_cartesian(xlim = xlim) +
    labs(x = "Rank in ranked gene list", y = y_lab) +
    .grs_panel_theme(base) +
    theme(legend.position = "none")
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' Compose the three panels
#'
#' @param panels List of three `ggplot`s, top to bottom.
#' @param panel_heights Length-3 numeric height ratio.
#' @param legend_position Legend placement.
#' @return A `patchwork`.
#' @keywords internal
.grs_compose <- function(panels, panel_heights, legend_position) {
  pw <- patchwork::wrap_plots(panels, ncol = 1L, heights = panel_heights)
  if (legend_position %in% c("right", "bottom")) {
    # A collected guide is placed by patchwork at figure level, so the
    # justification has to be re-asserted there with `&`. Set on the ES panel
    # alone it is ignored, and the legend floats level with the tick panel.
    pw <- pw + patchwork::plot_layout(guides = "collect")
    pw <- pw & theme(
      legend.position = legend_position,
      legend.justification = if (identical(legend_position, "right")) {
        "top"
      } else {
        "left"
      },
      legend.justification.right = "top"
    )
  }
  pw
}

#' Per-panel chrome, re-applied over any foundation theme
#'
#' @param base Foundation theme.
#' @return A list of theme objects, applied in order.
#' @keywords internal
.grs_panel_theme <- function(base) {
  list(
    base,
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      # Panels sit directly above one another, so per-panel top/bottom margins
      # would read as gaps in what should be one figure.
      plot.margin = margin(2, 10, 2, 5),
      # Bold axis titles and a bold plot title, and tighter label margins:
      # this figure is read from a distance on a poster or a slide more often
      # than it is read close up. Scoped to this renderer rather than to
      # theme_bulki(), so no other figure changes weight.
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
      axis.text.y = element_text(margin = margin(r = 2)),
      axis.title.x = element_text(face = "bold", margin = margin(t = 4)),
      axis.text.x = element_text(margin = margin(t = 2)),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7)
    )
  )
}

#' Legend placement
#'
#' @param legend_position One of `"inside"`, `"right"`, `"bottom"`, `"none"`.
#' @param legend_pos Length-2 numeric inside-legend position.
#' @return A ggplot2 theme.
#' @keywords internal
.grs_legend_theme <- function(legend_position, legend_pos) {
  switch(
    legend_position,
    inside = theme(
      legend.position = "inside",
      legend.position.inside = legend_pos,
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", colour = "grey90"),
      legend.key.size = unit(0.8, "lines"),
      legend.key.spacing.y = unit(2, "pt")
    ),
    # Top-aligned, not centred: the legend belongs beside the ES panel it
    # describes. Centred over the full figure height it floats level with the
    # tick panel, which it says nothing about.
    right = theme(legend.position = "right",
                  legend.justification = "top",
                  legend.justification.right = "top",
                  legend.key.spacing.y = unit(2, "pt")),
    bottom = theme(legend.position = "bottom",
                   legend.justification = "left",
                   legend.key.spacing.y = unit(2, "pt")),
    none = theme(legend.position = "none")
  )
}

#' Hide the x axis on a panel that is not the bottom one
#'
#' @return A ggplot2 theme.
#' @keywords internal
.grs_hide_x <- function() {
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
}

#' Resolve the rank vector for [gs_plot_running()]
#'
#' @param x The renderer's `x`.
#' @param ranks The renderer's `ranks`, possibly `NULL`.
#' @return A named numeric vector, sorted decreasing.
#' @keywords internal
.grs_ranks <- function(x, ranks) {
  ranks <- ranks %||% attr(x, "ranks")
  if (is.null(ranks)) {
    stop("`ranks` is required: supply the named numeric vector from ",
         "`gs_ranks()`. A `gs_result` does not carry the ranked list.",
         call. = FALSE)
  }
  if (!is.numeric(ranks) || is.null(names(ranks))) {
    stop("`ranks` must be a *named* numeric vector; see `gs_ranks()`.",
         call. = FALSE)
  }
  sort(ranks, decreasing = TRUE)
}

#' Resolve gene-set membership for [gs_plot_running()]
#'
#' @param x The renderer's `x`.
#' @param db The renderer's `db`, possibly `NULL`.
#' @return A named list of character vectors.
#' @keywords internal
.grs_sets <- function(x, db) {
  cand <- db %||% attr(x, "gene_sets")
  if (is.null(cand) && !inherits(x, "gs_result") && .grs_is_sets(x)) {
    cand <- x
  }
  if (is.null(cand)) {
    stop("`db` is required: the running curve needs full gene-set membership, ",
         "which a `gs_result` does not carry. Pass the `gs_db` you tested ",
         "against.", call. = FALSE)
  }
  if (!.grs_is_sets(cand)) {
    stop("`db` must be a `gs_db` or a named list of character vectors.",
         call. = FALSE)
  }
  sets <- unclass(cand)
  attributes(sets) <- list(names = names(sets))
  lapply(sets, as.character)
}

#' Is `x` shaped like a named list of gene sets?
#'
#' @param x Object to test.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
.grs_is_sets <- function(x) {
  is.list(x) && length(x) > 0L && !is.null(names(x)) &&
    all(nzchar(names(x))) &&
    all(vapply(x, function(z) is.character(z) || is.factor(z), logical(1L)))
}

#' Choose which pathways to draw
#'
#' @param x The renderer's `x`.
#' @param sets Named list of gene sets.
#' @param pathways Ids, row indices, or `NULL`.
#' @param top_n Integer default count.
#' @return A character vector of pathway ids, in plotting order.
#' @keywords internal
.grs_select <- function(x, sets, pathways, top_n) {
  if (is.null(pathways)) {
    if (inherits(x, "gs_result")) {
      if (!nrow(x)) {
        stop("`x` has no rows, so there is nothing to plot.", call. = FALSE)
      }
      ord <- order(abs(x$stat), decreasing = TRUE)
      ids <- x$pathway_id[ord][seq_len(min(top_n, length(ord)))]
    } else {
      ids <- names(sets)[seq_len(min(top_n, length(sets)))]
    }
  } else if (is.numeric(pathways)) {
    if (!inherits(x, "gs_result")) {
      stop("Integer `pathways` index the rows of a `gs_result`; pass pathway ",
           "ids instead.", call. = FALSE)
    }
    bad <- pathways < 1 | pathways > nrow(x)
    if (any(bad)) {
      stop("`pathways` index outside the ", nrow(x), " rows of `x`: ",
           paste(pathways[bad], collapse = ", "), ".", call. = FALSE)
    }
    ids <- x$pathway_id[as.integer(pathways)]
  } else {
    ids <- as.character(pathways)
  }
  ids <- unique(ids)
  missing_ids <- setdiff(ids, names(sets))
  if (length(missing_ids)) {
    stop("Pathway(s) not found in `db`: ",
         paste(encodeString(missing_ids, quote = "\""), collapse = ", "),
         ".", call. = FALSE)
  }
  ids
}

#' Legend labels, keyed by pathway id
#'
#' @param x The renderer's `x`.
#' @param ids Character vector of pathway ids.
#' @param labels The renderer's `labels`.
#' @param max_name_length Wrap width.
#' @return A character vector of labels, named by pathway id.
#' @keywords internal
.grs_labels <- function(x, ids, labels, max_name_length) {
  base <- stats::setNames(ids, ids)
  if (inherits(x, "gs_result")) {
    hit <- match(ids, x$pathway_id)
    ok <- !is.na(hit)
    base[ok] <- x$pathway_name[hit[ok]]
  } else {
    pn <- attr(x, "pathway_names")
    if (!is.null(pn)) {
      hit <- pn[ids]
      base[!is.na(hit)] <- unname(hit[!is.na(hit)])
    }
  }
  if (!is.null(labels)) {
    if (!is.null(names(labels))) {
      hit <- intersect(ids, names(labels))
      base[hit] <- as.character(labels[hit])
    } else {
      labels <- rep(as.character(labels), length.out = length(ids))
      base[] <- labels
    }
  }
  # Every other gs_plot_* renderer routes labels through format_pathway_name()
  # before a reader sees them; this one did not, so a db carrying no
  # `pathway_names` (or a `gs_result` whose `pathway_name` is still the id) put
  # raw snake_case MSigDB ids straight into the legend.
  #
  # Formatted selectively, not unconditionally: format_pathway_name() implements
  # smart capitalisation for ALL_CAPS_SNAKE ids and is *not* idempotent on text
  # that is already prose -- it turns "Beta response" into "beta Response". So a
  # label an explicit `labels =` or a provider's `pathway_names` already made
  # human-readable is left exactly as given, and only labels still equal to
  # their raw id are formatted. Fixing the raw-id leak must not introduce
  # mangled capitalisation in its place.
  # A placeholder is not a name. gsdb_from_file() used to hand through the
  # literal "-" that coresh_derived_sets.gmt writes in every description field,
  # which put "-" in the legend once per curve. Falling back to the id here as
  # well as at parse time means a cache built by an older provider still
  # renders a usable legend.
  raw <- .gs_placeholder_name(base) | base == names(base)
  if (any(raw)) {
    base[raw] <- format_pathway_name(names(base)[raw])
  }
  vapply(base, .grs_wrap, character(1L), width = max_name_length,
         USE.NAMES = TRUE)
}

#' Wrap a label onto several lines rather than truncating it
#'
#' @param s Character scalar.
#' @param width Maximum line width.
#' @return A character scalar, possibly containing newlines.
#' @keywords internal
.grs_wrap <- function(s, width) {
  if (is.na(s) || !nzchar(s)) return("")
  paste(strwrap(s, width = max(width, 8L)), collapse = "\n")
}

#' Resolve a palette keyed by pathway id
#'
#' Named palettes are matched by id; unnamed ones are zipped to `ids` in the
#' caller's declared order. The result is always named by id, which is what the
#' colour aesthetic maps, so ggplot's alphabetical level order cannot permute
#' the colours.
#'
#' @param palette The renderer's `palette`.
#' @param ids Character vector of pathway ids, in plotting order.
#' @return A character vector of colours, named by pathway id.
#' @keywords internal
.grs_palette <- function(palette, ids) {
  n <- length(ids)
  if (is.null(palette) || length(palette) == 0L) {
    # Assigned in sorted-id order, so the same id set always gets the same
    # colours whatever order the caller asked for them in. Stability across
    # *different* id sets needs gs_palette() over a fixed universe.
    return(gs_palette(ids))
  }
  nm <- names(palette)
  if (!is.null(nm) && all(ids %in% nm)) {
    return(stats::setNames(as.character(palette[ids]), ids))
  }
  if (!is.null(nm) && any(nzchar(nm))) {
    warning("`palette` names do not cover every plotted pathway id; zipping ",
            "colours to `pathways` order instead.", call. = FALSE)
  }
  stats::setNames(rep(as.character(palette), length.out = n), ids)
}

#' Panel keys and their axis labels
#'
#' @param metric_label Axis label for the ranked-metric panel.
#' @return A named character vector, `panel key -> axis label`.
#' @keywords internal
.grs_panel_levels <- function(metric_label = "Ranked metric") {
  if (!is.character(metric_label) || length(metric_label) != 1L) {
    stop("`metric_label` must be a single string.", call. = FALSE)
  }
  c(es = "Enrichment score", ticks = "Genes", stats = metric_label)
}

#' Running-curve data for each pathway
#'
#' Delegates to [fgsea::plotEnrichmentData()], which returns `curve`, `ticks`
#' and `stats` -- exactly the three panels. The cumulative sum is never
#' recomputed here.
#'
#' @param sets Named list of gene sets to draw.
#' @param ranks Named numeric vector, decreasing.
#' @param gsea_param `gseaParam` for fgsea.
#' @return A named list of `plotEnrichmentData()` outputs.
#' @keywords internal
.grs_curves <- function(sets, ranks, gsea_param) {
  out <- lapply(names(sets), function(id) {
    genes <- intersect(sets[[id]], names(ranks))
    if (!length(genes)) {
      stop("Pathway ", encodeString(id, quote = "\""),
           " has no genes in `ranks`, so it has no ",
           "running curve.", call. = FALSE)
    }
    fgsea::plotEnrichmentData(pathway = genes, stats = ranks,
                              gseaParam = gsea_param)
  })
  stats::setNames(out, names(sets))
}

#' Default foundation theme
#'
#' Uses `theme_bulki()` when the package provides it and falls back to
#' [ggplot2::theme_minimal()] so this file has no copy of it.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @keywords internal
.grs_default_theme <- function(base_size) {
  tb <- get0("theme_bulki", envir = asNamespace("bulkiRNA"),
             mode = "function", ifnotfound = NULL)
  if (!is.null(tb)) tb(base_size = base_size) else theme_minimal(base_size)
}
