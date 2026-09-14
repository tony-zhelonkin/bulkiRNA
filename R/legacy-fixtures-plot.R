# Legacy plot baseline fixtures are non-exported implementations retained for
# verification against the golden baseline captured at `752481f`.

#' Convert a bare `gseaResult` (or pass through a `gs_result`) for the plot
#' shims
#'
#' The frozen plot names (`gsea_dotplot()`, `gsea_dotplot_facet()`,
#' `gsea_barplot()`, `gsea_running_sum_plot()`) historically only ever took a
#' bare `gseaResult` -- no `database`/`contrast` metadata -- so they cannot go
#' through `normalize_gsea_results()`, whose old formals require both with no
#' default. For that case, this builds just enough of a `gs_result` to drive
#' the new renderers: `pathway_id`/`pathway_name` from `ID`/`Description`,
#' `stat`/`stat_type` from `NES`/`"NES"`, `p_value`/`padj` from
#' `pvalue`/`p.adjust`, `n_genes`/`n_genes_tested` both from `setSize` (the old
#' object carries no better proxy for genes-tested), and `leading_edge` from
#' splitting `core_enrichment` on `"/"`. `database` and `contrast` are filled
#' with a constant placeholder since the frozen callers never supplied either.
#'
#' But the single most common legacy call pattern is `run_gsea()` piped
#' straight into one of these four -- and `run_gsea()` is itself a deprecated
#' shim that now returns a `gs_result`, not an S4 `gseaResult`. That object
#' already carries `attr(x, "ranks")`/`attr(x, "gene_sets")` when it came from
#' `run_gsea()`, so it is passed through untouched rather than re-derived.
#'
#' @param gsea_obj A `gseaResult` object, or a `gs_result` (e.g. from the
#'   `run_gsea()` shim).
#' @return A `gs_result`. When `gsea_obj` is already a `gs_result` it is
#'   returned unchanged (attributes and all). When it is a `gseaResult`,
#'   `attr(., "ranks")` is set to `gsea_obj@geneList` and
#'   `attr(., "gene_sets")` to `gsea_obj@geneSets`, for `gs_plot_running()`'s
#'   fallback.
#' @keywords internal
.dep_gsea_to_gs_result <- function(gsea_obj) {
  if (inherits(gsea_obj, "gs_result")) {
    return(gsea_obj)
  }
  if (!methods::is(gsea_obj, "gseaResult")) {
    stop("This deprecated function requires a `gseaResult` or a `gs_result` ",
         "object; got `", paste(class(gsea_obj), collapse = "/"), "`.",
         call. = FALSE)
  }
  res <- as.data.frame(gsea_obj@result)
  core <- res$core_enrichment
  leading_edge <- lapply(strsplit(core, "/"), function(g) g[nzchar(g)])

  df <- data.frame(
    pathway_id = as.character(res$ID),
    pathway_name = as.character(res$Description),
    n_genes = as.integer(res$setSize),
    n_genes_tested = as.integer(res$setSize),
    stat = as.numeric(res$NES),
    p_value = as.numeric(res$pvalue),
    padj = as.numeric(res$p.adjust),
    es = as.numeric(res$enrichmentScore),
    stringsAsFactors = FALSE
  )
  df$leading_edge <- leading_edge

  out <- gs_result(
    df,
    database = "gsea",
    contrast = "gsea",
    method = "fgsea",
    stat_type = "NES"
  )
  attr(out, "ranks") <- gsea_obj@geneList
  attr(out, "gene_sets") <- gsea_obj@geneSets
  out
}

#' Enhanced GSEA dotplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_dot()] instead. This shim reproduces the old
#' formals of `gsea_dotplot()` verbatim and forwards to `gs_plot_dot()`.
#'
#' @param gsea_obj GSEA result object: a `gseaResult`, or a `gs_result`
#'   (e.g. from the `run_gsea()` shim).
#' @param filterBy Method to sort/filter results: `"p.adjust"` (default),
#'   `"NES"`, `"NES_positive"`, `"NES_negative"`.
#' @param sortBy Secondary display sort (`"GeneRatio"` or `"p.adjust"`).
#' @param showCategory Number of top pathways to show.
#' @param padj_cutoff Significance threshold for highlighting.
#' @param title Plot title.
#' @param wrap_width Width for text wrapping.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits Numeric vector of length 2 for symmetric NES limits.
#' @param min.dotSize Minimum dot size.
#' @param max.dotSize Maximum dot size.
#' @param highlight_sig Whether to highlight significant points.
#' @param highlight_threshold FDR threshold for highlighting; `NULL` uses
#'   `padj_cutoff`.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#' @param use_gradient Logical, use continuous gradient for NES.
#'
#' @return A ggplot object, as returned by [gs_plot_dot()].
#' @note The old dotplot put **GeneRatio** on the x axis and NES only in the
#'   fill gradient (`GeneRatio <- count / setSize`, where `count` is the
#'   leading-edge gene count). This shim forwards `aes_x = "gene_ratio"` to
#'   reproduce that; `gs_plot_dot()`'s own default (`aes_x = "stat"`) would put
#'   NES on x, which changes what the axis *means* rather than how it looks.
#'   `gs_plot_dot()`'s `gene_ratio` is `lengths(leading_edge) / n_genes`, which
#'   is exactly `count / setSize` on both paths into this shim: a converted
#'   `gseaResult` (where `n_genes` comes from `setSize` and `leading_edge` from
#'   splitting `core_enrichment`) and a `gs_result` from the `run_gsea()` shim
#'   (`n_genes` from `lengths(sets[pathway])`, `leading_edge` from fgsea's
#'   `leadingEdge`). An exact match, not an approximation.
#' @note `sortBy` was the old *display* sort (its stage 3), not a selection
#'   tie-break. `gs_plot_dot()` orders the y axis by the chosen `sort_by` within
#'   the selection, so `sortBy` is accepted and ignored.
#' @note `use_gradient = FALSE` (the old "binary colour" mode) has no
#'   counterpart: `gs_plot_dot()` always uses the continuous diverging
#'   gradient. Accepted and ignored; the plot always renders with a gradient.
#' @keywords internal
gsea_dotplot <- function(
    gsea_obj,
    filterBy = "p.adjust",
    sortBy = "GeneRatio",
    showCategory = 10,
    padj_cutoff = 0.05,
    title = "GSEA Dotplot",
    wrap_width = 50,
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    highlight_threshold = NULL,
    strip_prefix = TRUE,
    use_gradient = TRUE,
    nes_limits = NULL) {
  .Deprecated("gs_plot_dot")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  sort_by <- "padj"
  direction <- "both"
  if (identical(filterBy, "NES")) {
    sort_by <- "stat"
  } else if (identical(filterBy, "NES_positive")) {
    sort_by <- "stat"
    direction <- "up"
  } else if (identical(filterBy, "NES_negative")) {
    sort_by <- "stat"
    direction <- "down"
  }

  highlight <- if (isTRUE(highlight_sig)) {
    highlight_threshold %||% padj_cutoff
  } else {
    NULL
  }

  gs_plot_dot(
    x,
    top_n = showCategory,
    aes_x = "gene_ratio",
    sort_by = sort_by,
    direction = direction,
    highlight = highlight,
    size_range = c(min.dotSize, max.dotSize),
    limits = nes_limits,
    palette = c(neg_color, mid_color, pos_color),
    wrap_width = wrap_width,
    strip_prefix = strip_prefix,
    title = title
  )
}

#' Enhanced GSEA faceted dotplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_dot()] with `facet = "direction"` instead. This
#' shim reproduces the old formals of `gsea_dotplot_facet()` verbatim.
#'
#' @param gsea_obj GSEA result object: a `gseaResult`, or a `gs_result`
#'   (e.g. from the `run_gsea()` shim).
#' @param showCategory Number of pathways to show per direction.
#' @param padj_cutoff Adjusted p-value cutoff used for significance
#'   highlighting.
#' @param title Plot title.
#' @param wrap_width Width for text wrapping.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits NES colour scale limits.
#' @param min.dotSize Minimum dot size.
#' @param max.dotSize Maximum dot size.
#' @param highlight_sig Whether to highlight significant points.
#' @param highlight_threshold FDR threshold for highlighting; `NULL` uses
#'   `padj_cutoff`.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#'
#' @return A ggplot object, as returned by [gs_plot_dot()].
#' @note Forwards `aes_x = "gene_ratio"`, matching the old faceted dotplot's
#'   GeneRatio x axis -- see the `gsea_dotplot()` note for the derivation.
#' @keywords internal
gsea_dotplot_facet <- function(
    gsea_obj,
    showCategory = 10,
    padj_cutoff = 0.05,
    title = "GSEA Faceted Dotplot",
    wrap_width = 50,
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    nes_limits = c(-3.5, 3.5),
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    highlight_threshold = NULL,
    strip_prefix = TRUE) {
  .Deprecated("gs_plot_dot")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  highlight <- if (isTRUE(highlight_sig)) {
    highlight_threshold %||% padj_cutoff
  } else {
    NULL
  }

  gs_plot_dot(
    x,
    top_n = showCategory,
    aes_x = "gene_ratio",
    sort_by = "padj",
    direction = "both",
    facet = "direction",
    highlight = highlight,
    size_range = c(min.dotSize, max.dotSize),
    limits = nes_limits,
    palette = c(neg_color, mid_color, pos_color),
    wrap_width = wrap_width,
    strip_prefix = strip_prefix,
    title = title
  )
}

#' Enhanced GSEA barplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_bar()] instead. This shim reproduces the old
#' formals of `gsea_barplot()` verbatim.
#'
#' @param gsea_obj GSEA result object: a `gseaResult`, or a `gs_result`
#'   (e.g. from the `run_gsea()` shim).
#' @param padj_cutoff Adjusted p-value cutoff; a hard filter, as in the old
#'   function.
#' @param top_n Number of pathways to show.
#' @param title Plot title.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits NES colour scale limits.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#'
#' @return A ggplot object, as returned by [gs_plot_bar()].
#' @note The old function pre-filtered on `padj_cutoff` and never outlined
#'   individual bars, so `highlight` is forced to `NULL` here -- outlining an
#'   already-significance-filtered set would be a no-op difference from the
#'   old figure, but is called out explicitly rather than left implicit.
#' @note After `gs_plot_bar()` returns, the plot's `$data` and `gs_source`
#'   attribute are re-sorted to ascending `stat`, matching the old function's
#'   post-selection `dplyr::arrange(NES)`. This changes only row order in the
#'   data and the table `gs_save()` writes -- bar positions come from the y-axis
#'   factor's levels and are already correct.
#' @keywords internal
gsea_barplot <- function(
    gsea_obj,
    padj_cutoff = 0.05,
    top_n = 30,
    title = "GSEA NES Barplot",
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    nes_limits = c(-3.5, 3.5),
    strip_prefix = TRUE) {
  .Deprecated("gs_plot_bar")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  p <- gs_plot_bar(
    x,
    top_n = top_n,
    sort_by = "stat",
    direction = "both",
    padj_max = padj_cutoff,
    highlight = NULL,
    limits = nes_limits,
    palette = c(neg_color, mid_color, pos_color),
    strip_prefix = strip_prefix,
    title = title
  )

  # `gs_plot_bar()` selects top-by-|stat|, matching the old
  # `arrange(desc(abs(NES))) |> head(top_n)`, and sets the y-axis factor's
  # *levels* in ascending-stat order -- levels alone decide bar position, so
  # nothing rendered depends on row order. The old function additionally
  # reordered the rows themselves (`arrange(NES)`, `gsea_barplot.R:64`), which
  # is what `gs_save()` writes to the table beside the figure. Without this the
  # table would list pathways in selection order (steepest |NES| first) instead.
  if (nrow(p$data) > 0L) {
    ord <- order(p$data$stat)
    p$data <- p$data[ord, , drop = FALSE]
    src <- attr(p, "gs_source")
    if (!is.null(src)) attr(p, "gs_source") <- src[ord, , drop = FALSE]
  }
  p
}

#' Unified GSEA running-sum plot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_running()] instead. This shim reproduces the old
#' formals of `gsea_running_sum_plot()` verbatim.
#'
#' @param gsea_obj A `gseaResult` object from clusterProfiler/fgsea, or a
#'   `gs_result` (e.g. from the `run_gsea()` shim) -- passed through
#'   untouched, so its own `ranks`/`gene_sets` attributes drive the curve.
#' @param gene_set_ids Integer vector (row indices) or character vector
#'   (pathway IDs). `NULL` (default) picks the new renderer's default top 5
#'   by \eqn{|NES|}.
#' @param palette Optional colour vector, keyed as documented in
#'   [gs_plot_running()].
#' @param labels Optional legend labels, keyed as documented in
#'   [gs_plot_running()].
#' @param legend_pos Inside-panel legend position (used when
#'   `legend_position = "inside"`).
#' @param base_size Base font size for the per-panel theme.
#' @param max_name_length Maximum character length for pathway names in the
#'   legend.
#' @param title Optional title for the ES (top) panel.
#' @param panel_heights Length-3 numeric, the ES:rug:metric height ratio.
#' @param legend_position One of `"inside"` (default), `"right"` or `"none"`.
#' @param es_ylim Optional length-2 numeric y limit for the ES panel.
#' @param xticks `"bottom"` (default) or `"all"`.
#' @param rug_ylabels Logical. Show the rug panel's lane indices.
#' @param base_theme Optional complete ggplot2 theme used as the foundation
#'   for every panel.
#'
#' @return A `patchwork`, as returned by [gs_plot_running()].
#' @note `es_ylim`, `xticks` and `rug_ylabels` are forwarded again as of
#'   1.1.0. They were absorbed and ignored while `gs_plot_running()` drew all
#'   three panels as facets of a single plot, which had no per-panel y scale to
#'   honour them with. Now that it composes three real panels, each of these
#'   formals means what it originally meant.
#' @note `gsea_obj@geneList` and `gsea_obj@geneSets` are attached to the
#'   converted `gs_result` as the `ranks` and `gene_sets` attributes, which
#'   [gs_plot_running()] falls back to when its own `ranks`/`db` arguments are
#'   `NULL` -- so this shim never duplicates that lookup logic itself.
#' @keywords internal
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids = NULL,
                                  palette = NULL,
                                  labels = NULL,
                                  legend_pos = c(.98, .98),
                                  base_size = 14,
                                  max_name_length = 40,
                                  title = NULL,
                                  panel_heights = c(2.4, 0.7, 0.9),
                                  legend_position = c("inside", "right", "none"),
                                  es_ylim = NULL,
                                  xticks = c("bottom", "all"),
                                  rug_ylabels = FALSE,
                                  base_theme = NULL) {
  .Deprecated("gs_plot_running")
  legend_position <- match.arg(legend_position)
  xticks <- match.arg(xticks)

  x <- .dep_gsea_to_gs_result(gsea_obj)

  gs_plot_running(
    x,
    pathways = gene_set_ids,
    labels = labels,
    palette = palette,
    panel_heights = panel_heights,
    es_ylim = es_ylim,
    title = title,
    base_size = base_size,
    max_name_length = max_name_length,
    legend_position = legend_position,
    legend_pos = legend_pos,
    xticks = xticks,
    rug_ylabels = rug_ylabels,
    base_theme = base_theme
  )
}

#' Custom minimal theme with grid (deprecated)
#'
#' @description
#' Deprecated: use [theme_bulki()] instead. This shim reproduces the old
#' formals of `custom_minimal_theme_with_grid()` verbatim and forwards to
#' the grid variant of the `theme_bulki()` theme.
#'
#' @param base_size Base font size for the theme.
#' @param base_family Base font family for the theme.
#'
#' @return A ggplot2 theme object.
#' @note The old theme had no base-size floor and rendered at its documented
#'   12 pt default; `theme_bulki()` raises anything below 14 pt to 14 pt on
#'   purpose. To keep this deprecated path byte-faithful, the shim routes
#'   through the internal `.theme_bulki(floor = NULL)` rather than the public
#'   `theme_bulki()`, so `base_size = 12` really does render at 12 pt. The
#'   public theme keeps its floor. No gate could have caught the difference:
#'   the golden record for this case stores only the 144 theme element *names*,
#'   never their sizes.
#' @keywords internal
custom_minimal_theme_with_grid <- function(base_size = 12, base_family = "") {
  .Deprecated("theme_bulki")
  .theme_bulki(base_size = base_size, base_family = base_family,
               grid = TRUE, floor = NULL)
}

#' Plot all GSEA results for a given analysis (deprecated)
#'
#' @description
#' Deprecated: the new golden path is to call the [gs_plot_dot()] /
#' [gs_plot_bar()] renderer you want and save it with [gs_save()]. This shim
#' reproduces the old formals of `plot_all_gsea_results()` verbatim and
#' forwards to the internal `.gs_plot_all()`, which holds the ported body.
#'
#' @param gsea_list List of GSEA results: `gseaResult` objects, `gs_result`
#'   objects (e.g. from the `run_gsea()` shim), or a mix of both.
#' @param analysis_name Name of the analysis.
#' @param out_root Output root directory.
#' @param n_pathways Number of pathways to show.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param expr_data Expression data matrix for heatmaps. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param sample_annotation Sample annotation data.frame. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param sample_order Order of samples for plots. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param ann_colors Colours for annotations. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#'
#' @return A character vector of every written path, invisibly, as returned
#'   by `.gs_plot_all()`.
#' @note `gsea_list` (a named list of `gseaResult`/data-frame objects keyed by
#'   database) is converted to a single `gs_result` via the internal
#'   `.dep_gsea_to_gs_result()`, with each element's list name substituted in
#'   as its `database` column so `.gs_plot_all()`'s per-database loop still
#'   works.
#' @note `expr_data`, `sample_annotation`, `sample_order` and `ann_colors`
#'   drove heatmap panels in the old function; `.gs_plot_all()` (the ported
#'   body behind this shim) does not render a heatmap, so these four formals
#'   are accepted and ignored. No heatmap is written where the old function
#'   would have written one.
#' @keywords internal
plot_all_gsea_results <- function(
    gsea_list,
    analysis_name,
    out_root,
    n_pathways = 20,
    padj_cutoff = 0.05,
    expr_data = NULL,
    sample_annotation = NULL,
    sample_order = NULL,
    ann_colors = NULL) {
  # `.gs_plot_all()` is internal, so naming it was unactionable advice.
  .Deprecated(msg = paste(
    "`plot_all_gsea_results()` is deprecated. Build the figures you want with",
    "gs_plot_dot() / gs_plot_bar() / gs_plot_running() and write them with",
    "gs_save(), which also writes each figure's source table."))

  if (is.null(gsea_list) || length(gsea_list) == 0) {
    return(invisible(character(0)))
  }

  parts <- lapply(names(gsea_list), function(nm) {
    part <- .dep_gsea_to_gs_result(gsea_list[[nm]])
    part[["database"]] <- nm
    part
  })
  x <- do.call(rbind, parts)

  .gs_plot_all(
    x,
    out_dir = file.path(out_root, analysis_name),
    name = analysis_name,
    top_n = n_pathways,
    padj_cutoff = padj_cutoff
  )
}

#' Save GSEA results to a text log file (deprecated)
#'
#' @description
#' Deprecated: use [gs_save()] (source-table export) or the internal
#' `.gs_write_log()` directly. This shim reproduces the old formals of
#' `save_gsea_log()` verbatim and forwards to the internal `.gs_write_log()`,
#' which holds the ported body.
#'
#' @param gsea_obj GSEA result object: a `gseaResult`, or a `gs_result`
#'   (e.g. from the `run_gsea()` shim).
#' @param filename Output filename (with or without path).
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param dir Output directory (optional).
#'
#' @return `path`, invisibly, as returned by `.gs_write_log()`.
#' @keywords internal
save_gsea_log <- function(
    gsea_obj,
    filename,
    padj_cutoff = 0.05,
    dir = NULL) {
  # `.gs_write_log()` is internal, so naming it was unactionable advice.
  .Deprecated(msg = paste(
    "`save_gsea_log()` is deprecated. gs_save() writes the source table",
    "beside each figure; there is no exported replacement for the free-text",
    "log, by design."))

  if (!is.null(dir)) {
    filename <- file.path(dir, filename)
  }
  x <- .dep_gsea_to_gs_result(gsea_obj)
  .gs_write_log(x, filename, padj_cutoff = padj_cutoff)
}

#' Create a standard volcano plot (deprecated)
#'
#' @description
#' Deprecated: use [de_volcano()] instead. This shim reproduces the old
#' formals of `create_standard_volcano()` verbatim and forwards to
#' `de_volcano()`. The frozen `color_palette` formal forwards to the new
#' `palette` formal. The new renderer's `orientation` argument remains after
#' `...`, so no positional call to this shim can shift onto it.
#'
#' @param de_results Data frame whose rownames are gene IDs and that contains
#'   at least `logFC`, `P.Value`, `adj.P.Val`.
#' @param decision_by Choose adjusted (`"fdr"`) or raw (`"p"`) p-values.
#' @param p_cutoff Numeric significance threshold.
#' @param fc_cutoff Numeric absolute log2 fold-change threshold.
#' @param top_n Integer genes labelled per side.
#' @param highlight_gene Character vector of gene IDs always labelled.
#' @param label_method Label-selection method.
#' @param x_breaks Numeric fold-change axis spacing.
#' @param title Plot title.
#' @param subtitle Optional plot subtitle.
#' @param caption Optional plot caption.
#' @param fixed_p_boundary Optional raw p-value for the FDR boundary line.
#' @param color_palette Named character vector of four category colours.
#' @param show_grid Logical. Keep the panel grid.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#' @param annotate_counts Logical. Add up/down counts to the legend.
#' @param ... Absorbs deprecated arguments such as `use_fdr`.
#'
#' @return A `ggplot2` object.
#' @keywords internal
create_standard_volcano <- function(
    de_results,
    decision_by = c("fdr", "p"),
    p_cutoff = 0.05,
    fc_cutoff = 2,
    top_n = 5,
    highlight_gene = NULL,
    label_method = "top",
    x_breaks = 1,
    title = "Volcano plot",
    subtitle = NULL,
    caption = NULL,
    fixed_p_boundary = NULL,
    color_palette = c(
      "NS"               = "#7F7F7F",
      "Log2FC"           = "#0173B2",
      "p-value"          = "#029E73",
      "p-value & Log2FC" = "#D55E00"
    ),
    show_grid = FALSE,
    max.overlaps = 10,
    annotate_counts = FALSE,
    ...) {
  .Deprecated("de_volcano")
  de_volcano(
    de_results = de_results,
    decision_by = decision_by,
    p_cutoff = p_cutoff,
    fc_cutoff = fc_cutoff,
    top_n = top_n,
    highlight_gene = highlight_gene,
    label_method = label_method,
    x_breaks = x_breaks,
    title = title,
    subtitle = subtitle,
    caption = caption,
    fixed_p_boundary = fixed_p_boundary,
    palette = color_palette,
    show_grid = show_grid,
    max.overlaps = max.overlaps,
    annotate_counts = annotate_counts,
    ...
  )
}

#' Create a mean-difference (MD) plot (deprecated)
#'
#' @description
#' Deprecated: use [de_md_plot()] instead. This shim reproduces the old
#' formals of `create_MD_plot()` verbatim and forwards to `de_md_plot()`.
#' The frozen `color_palette` formal forwards to the new `palette` formal;
#' `fdr_cutoff` is unchanged.
#'
#' @param fit An `MArrayLM` object from limma.
#' @param coef Integer index or character name of the coefficient to plot.
#' @param de_results Optional matching `topTable()` data frame.
#' @param fc_cutoff Numeric absolute log2 fold-change guide.
#' @param fdr_cutoff Numeric FDR threshold for the Up/Down call.
#' @param top_n Integer genes labelled per direction.
#' @param highlight_gene Character vector of gene IDs always labelled.
#' @param label_method Label-selection method.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#' @param title Plot title.
#' @param color_palette Named colours for `Up`, `Down`, and `NS`.
#' @param show_grid Logical. Keep the panel grid.
#' @param show_quadrant_counts Logical. Annotate significant-gene counts.
#'
#' @return A `ggplot2` object.
#' @keywords internal
create_MD_plot <- function(
    fit,
    coef,
    de_results = NULL,
    fc_cutoff = 1,
    fdr_cutoff = 0.05,
    top_n = 5,
    highlight_gene = NULL,
    label_method = "top",
    max.overlaps = 10,
    title = NULL,
    color_palette = c(
      Up   = "#D55E00",
      Down = "#0072B2",
      NS   = "#999999"
    ),
    show_grid = FALSE,
    show_quadrant_counts = TRUE) {
  .Deprecated("de_md_plot")
  de_md_plot(
    fit = fit,
    coef = coef,
    de_results = de_results,
    fc_cutoff = fc_cutoff,
    fdr_cutoff = fdr_cutoff,
    top_n = top_n,
    highlight_gene = highlight_gene,
    label_method = label_method,
    max.overlaps = max.overlaps,
    title = title,
    palette = color_palette,
    show_grid = show_grid,
    show_quadrant_counts = show_quadrant_counts
  )
}
