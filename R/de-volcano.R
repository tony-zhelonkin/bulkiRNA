#' Up/down counts for the highest-priority populated legend category
#'
#' Walks the colour categories in priority order and returns the up/down split
#' of the **first** one that has any members, so the rendered counts always
#' belong to the strongest claim the legend actually makes:
#'
#' 1. `"both"` -- combined FDR & `|log2FC|` (`sig_dec & sig_fc`).
#' 2. `"fdr"` -- FDR only (`sig_dec`), used when priority 1 is empty.
#' 3. `"none"` -- nothing crosses FDR; counts are `0 / 0`.
#'
#' `legend_key` names the colour category whose legend label should receive the
#' appended count line. Up counts `logFC > 0`, down counts `logFC < 0`.
#'
#' @param df Data frame carrying `logFC`, `sig_dec`, `sig_fc`.
#' @return A list with `category`, `n_up`, `n_down`, `legend_key`.
#' @keywords internal
.volcano_sig_counts <- function(df) {
  if (any(df$sig_dec & df$sig_fc, na.rm = TRUE)) {
    sel <- df$sig_dec & df$sig_fc
    category   <- "both"
    legend_key <- "p-value & Log2FC"
  } else if (any(df$sig_dec, na.rm = TRUE)) {
    sel <- df$sig_dec
    category   <- "fdr"
    legend_key <- "p-value"
  } else {
    return(list(category = "none", n_up = 0L, n_down = 0L,
                legend_key = NA_character_))
  }
  list(category   = category,
       n_up       = sum(sel & df$logFC > 0, na.rm = TRUE),
       n_down     = sum(sel & df$logFC < 0, na.rm = TRUE),
       legend_key = legend_key)
}

#' Volcano plot for differential-expression results
#'
#' Plots `-log10(P.Value)` against `logFC` while **deciding significance with
#' `adj.P.Val`**. Raw *p* is a per-gene test statistic and gives maximal
#' resolution; FDR is an average over the rejected set and produces stair-step
#' artefacts when many *p* values collapse onto one adjusted value. Keeping raw
#' *p* on the axis and FDR for the decision buys the multiple-testing guarantee
#' without the staircase.
#'
#' The dashed horizontal line is therefore **not** at `-log10(p_cutoff)` in FDR
#' mode. It sits at the raw *p* value that realises the FDR boundary -- the
#' largest raw *p* among the genes that pass `adj.P.Val <= p_cutoff` -- so the
#' line always coincides with the colour boundary. When no gene passes, no line
#' is drawn and an italic annotation says so instead.
#'
#' | `decision_by` | Colour decision uses | Dashed line at |
#' |---|---|---|
#' | `"fdr"` (default) | `adj.P.Val <= p_cutoff` | `-log10` of the largest raw *p* with `adj.P.Val <= p_cutoff` |
#' | `"p"` | `P.Value <= p_cutoff` | `-log10(p_cutoff)` |
#'
#' `fixed_p_boundary` pins the line to a known raw *p* instead of deriving it.
#' Use it when `de_results` has already been filtered, so the largest surviving
#' raw *p* no longer marks the true FDR boundary.
#'
#' @param de_results Data frame whose **rownames are gene IDs**, carrying at
#'   least `logFC`, `P.Value` and `adj.P.Val` (limma `topTable()` shape).
#' @param decision_by `"fdr"` (default) or `"p"`; which column makes the
#'   significance call and therefore the point colour.
#' @param p_cutoff Numeric. An FDR threshold when `decision_by = "fdr"`, a raw
#'   *p* threshold when `decision_by = "p"`.
#' @param fc_cutoff Numeric. Absolute log2 fold-change threshold.
#' @param top_n Integer. Genes labelled per side, ranked by the deciding
#'   statistic.
#' @param highlight_gene Character vector of gene IDs always labelled, in bold
#'   black, with overlap suppression disabled.
#' @param label_method One of `"top"`, `"sig"`, `"p"`, `"log2fc"`; anything else
#'   labels nothing.
#' @param x_breaks Numeric spacing between fold-change axis ticks, in both
#'   orientations -- it follows the fold-change axis, not the *x* position. The
#'   -log10(p) axis always takes ggplot2's default breaks.
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#' @param caption Optional caption. `NULL` generates one documenting the
#'   thresholds actually drawn.
#' @param fixed_p_boundary Numeric raw *p* to pin the FDR boundary line to, or
#'   `NULL` to derive it from the significant set.
#' @param palette Named character vector of four colours for `NS`,
#'   `Log2FC`, `p-value` and `p-value & Log2FC`.
#' @param show_grid Logical. Keep the panel grid.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#' @param annotate_counts Logical. Append a second legend line of up/down
#'   counts under the highest-priority populated significance category.
#' @param ... Absorbs deprecated arguments such as `use_fdr`.
#' @param orientation `"horizontal"` (default; fold change on *x*) or
#'   `"vertical"` (fold change on *y*). Vertical panels stack into a row well --
#'   see [de_volcano_grid()]. Must be supplied by name.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#' @examples
#' de <- data.frame(
#'   logFC     = c(-3, -1, 0.2, 1.5, 4),
#'   P.Value   = c(1e-6, 0.02, 0.7, 0.01, 1e-8),
#'   adj.P.Val = c(1e-4, 0.06, 0.9, 0.04, 1e-6),
#'   row.names = paste0("Gene", 1:5)
#' )
#' de_volcano(de, p_cutoff = 0.05, fc_cutoff = 1)
de_volcano <- function(
    de_results,
    decision_by      = c("fdr", "p"),
    p_cutoff         = 0.05,
    fc_cutoff        = 2,
    top_n            = 5,
    highlight_gene   = NULL,
    label_method     = "top",
    x_breaks         = 1,
    title            = "Volcano plot",
    subtitle         = NULL,
    caption          = NULL,
    fixed_p_boundary = NULL,
    palette          = c(
      "NS"               = "#7F7F7F",
      "Log2FC"           = "#0173B2",
      "p-value"          = "#029E73",
      "p-value & Log2FC" = "#D55E00"
    ),
    show_grid        = FALSE,
    max.overlaps     = 10,
    annotate_counts  = FALSE,
    ...,
    orientation      = c("horizontal", "vertical")) {
  decision_by <- match.arg(decision_by)
  orientation <- match.arg(orientation)

  req <- c("logFC", "P.Value", "adj.P.Val")
  missing <- setdiff(req, colnames(de_results))
  if (length(missing)) {
    stop("`de_results` is missing column(s): ",
         paste(sprintf("`%s`", missing), collapse = ", "),
         ". A limma `topTable()` data frame has all of them.", call. = FALSE)
  }

  if (isTRUE(list(...)[["use_fdr"]])) {
    warning("`use_fdr` is deprecated: use `decision_by = \"fdr\"` instead.",
            call. = FALSE)
  }

  # ---- 1. significance and the boundary line -------------------------------
  if (decision_by == "fdr") {
    sig_stat  <- de_results$adj.P.Val
    sig_logic <- sig_stat <= p_cutoff              # inclusive
    if (!is.null(fixed_p_boundary)) {
      sig_line  <- -log10(fixed_p_boundary)
      draw_line <- TRUE
    } else {
      # `[sig_logic]` on a logical carrying NA keeps the NA positions, so this
      # was non-empty even when nothing was significant; `max(na.rm = TRUE)`
      # then returned -Inf and the line became NaN with draw_line = TRUE,
      # captioning the figure "p <= NaN" and making the documented "no genes
      # pass FDR" branch unreachable. `which()` drops the NAs.
      sig_pvals <- de_results$P.Value[which(sig_logic)]
      sig_pvals <- sig_pvals[!is.na(sig_pvals)]
      if (length(sig_pvals) > 0) {
        # The line goes at the largest raw p among the FDR-significant genes,
        # so it lands exactly on the colour boundary.
        sig_line  <- -log10(max(sig_pvals, na.rm = TRUE))
        draw_line <- TRUE
      } else {
        sig_line  <- NA_real_
        draw_line <- FALSE
      }
    }
    legend_sig <- sprintf("FDR \u2264 %.2g", p_cutoff)
  } else {
    sig_stat   <- de_results$P.Value
    sig_logic  <- sig_stat <= p_cutoff
    sig_line   <- -log10(p_cutoff)
    draw_line  <- TRUE
    legend_sig <- sprintf("p \u2264 %.2g", p_cutoff)
  }

  df <- dplyr::mutate(
    de_results,
    sig_fc  = abs(.data$logFC) >= fc_cutoff,
    sig_dec = sig_logic,
    significance_value = sig_stat,
    cat = dplyr::case_when(
      .data$sig_fc & .data$sig_dec ~ "p-value & Log2FC",
      .data$sig_fc                 ~ "Log2FC",
      .data$sig_dec                ~ "p-value",
      TRUE                         ~ "NS"
    )
  )

  # ---- 2. label selection ---------------------------------------------------
  get_top <- function(side) {
    keep <- if (side == "up") df$logFC > 0 else df$logFC < 0
    sub <- df[which(keep), , drop = FALSE]
    sub <- sub[order(sub$significance_value), , drop = FALSE]
    utils::head(sub, top_n)
  }

  lab_df <- switch(
    label_method,
    top    = rbind(get_top("up"), get_top("down")),
    sig    = df[df$cat == "p-value & Log2FC", , drop = FALSE],
    p      = df[df$sig_dec, , drop = FALSE],
    log2fc = df[df$sig_fc, , drop = FALSE],
    df[0, , drop = FALSE]
  )

  if (!is.null(highlight_gene)) {
    extra  <- df[rownames(df) %in% highlight_gene, , drop = FALSE]
    lab_df <- lab_df[!rownames(lab_df) %in% highlight_gene, , drop = FALSE]
  } else {
    extra <- df[0, , drop = FALSE]
  }

  # ---- 3. limits and colours ------------------------------------------------
  fc_tick <- ceiling(max(abs(df$logFC), na.rm = TRUE) / x_breaks) * x_breaks
  # A p-value that underflowed to 0 made `-log10()` infinite, and that Inf went
  # straight into `coord_cartesian(ylim = c(0, p_max))`, squashing every point
  # at the bottom of a blank panel. Clamp at the smallest representable double.
  p_top   <- max(-log10(pmax(df$P.Value, .Machine$double.xmin)), na.rm = TRUE)
  p_max   <- ceiling(p_top)
  dark_pal <- .de_shade(palette)

  legend_labels <- c(
    "p-value & Log2FC" = sprintf("%s & |log2FC| \u2265 %.1f", legend_sig, fc_cutoff),
    "Log2FC"           = sprintf("|log2FC| \u2265 %.1f", fc_cutoff),
    "p-value"          = legend_sig,
    "NS"               = "NS"
  )
  if (isTRUE(annotate_counts)) {
    sc  <- .volcano_sig_counts(df)
    key <- if (is.na(sc$legend_key)) "p-value & Log2FC" else sc$legend_key
    legend_labels[key] <- paste0(
      legend_labels[key],
      sprintf("\n\u2191 %d   \u2193 %d", sc$n_up, sc$n_down)
    )
  }

  auto_caption <- if (decision_by == "fdr") {
    if (draw_line) {
      sprintf(
        if (orientation == "horizontal") {
          "Dashed lines: horiz. \u2013 FDR \u2264 %.2g (p \u2264 %.2g); vert. \u2013 |log2FC| \u2265 %.1f"
        } else {
          "Dashed lines: vert. \u2013 FDR \u2264 %.2g (p \u2264 %.2g); horiz. \u2013 |log2FC| \u2265 %.1f"
        },
        p_cutoff, signif(10^(-sig_line), 2), fc_cutoff)
    } else {
      sprintf(
        if (orientation == "horizontal") {
          "No genes pass FDR \u2264 %.2g. Dashed lines: vert. \u2013 |log2FC| \u2265 %.1f"
        } else {
          "No genes pass FDR \u2264 %.2g. Dashed lines: horiz. \u2013 |log2FC| \u2265 %.1f"
        },
        p_cutoff, fc_cutoff)
    }
  } else {
    sprintf(
      if (orientation == "horizontal") {
        "Dashed lines: horiz. \u2013 p \u2264 %.2g; vert. \u2013 |log2FC| \u2265 %.1f"
      } else {
        "Dashed lines: vert. \u2013 p \u2264 %.2g; horiz. \u2013 |log2FC| \u2265 %.1f"
      },
      p_cutoff, fc_cutoff)
  }

  # ---- 4. build -------------------------------------------------------------
  if (orientation == "horizontal") {
    g <- ggplot(df, aes(x = .data$logFC, y = -log10(.data$P.Value),
                        colour = .data$cat)) +
      geom_point(size = 2, alpha = 0.65) +
      geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
      scale_x_continuous(breaks = seq(-fc_tick, fc_tick, by = x_breaks),
                         limits = c(-fc_tick, fc_tick)) +
      coord_cartesian(ylim = c(0, p_max)) +
      labs(x = "log2(FC)", y = expression(-log[10](p - value)))
  } else {
    g <- ggplot(df, aes(x = -log10(.data$P.Value), y = .data$logFC,
                        colour = .data$cat)) +
      geom_point(size = 2, alpha = 0.65) +
      geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
      # `x_breaks` is documented as the fold-change tick spacing, but in this
      # orientation fold change is on *y*: it used to retune the -log10(p) axis
      # instead, so `x_breaks = 0.5` on a dataset reaching p = 1e-40 asked for
      # ~80 p-axis ticks while the fold-change axis it names got ggplot2's
      # defaults. Both orientations now put `x_breaks` on the fold-change axis
      # and take the breaks-aligned tick maximum as its limit.
      scale_y_continuous(breaks = seq(-fc_tick, fc_tick, by = x_breaks),
                         limits = c(-fc_tick, fc_tick)) +
      coord_cartesian(xlim = c(0, p_max)) +
      labs(y = "log2(FC)", x = expression(-log[10](p - value)))
  }

  g <- g +
    scale_colour_manual(name = NULL, values = palette,
                        breaks = names(palette),
                        labels = legend_labels) +
    labs(title = title, subtitle = subtitle,
         caption = if (is.null(caption)) auto_caption else caption) +
    .de_theme()

  if (draw_line) {
    g <- g + if (orientation == "horizontal") {
      geom_hline(yintercept = sig_line, linetype = "dashed")
    } else {
      geom_vline(xintercept = sig_line, linetype = "dashed")
    }
  } else if (decision_by == "fdr") {
    g <- g + annotate(
      "text",
      x = if (orientation == "horizontal") fc_tick * 0.5 else p_max * 0.5,
      y = if (orientation == "horizontal") p_max * 0.95 else 0,
      label = sprintf("No genes pass FDR \u2264 %.2g", p_cutoff),
      size = 4, colour = "darkred", fontface = "italic")
  }

  if (!show_grid) {
    g <- g + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  }

  # ---- 5. labels ------------------------------------------------------------
  lab_aes <- if (orientation == "horizontal") {
    aes(x = .data$logFC, y = -log10(.data$P.Value), label = .data$.gene_label)
  } else {
    aes(x = -log10(.data$P.Value), y = .data$logFC, label = .data$.gene_label)
  }

  if (nrow(lab_df)) {
    lab_df$.gene_label <- rownames(lab_df)
    g <- g + ggrepel::geom_text_repel(
      data = lab_df, mapping = lab_aes,
      colour = unname(dark_pal[lab_df$cat]), fontface = "plain",
      size = 3.5, box.padding = 0.4, point.padding = 0.3,
      max.overlaps = max.overlaps, min.segment.length = 0,
      inherit.aes = FALSE, show.legend = FALSE)
  }
  if (nrow(extra)) {
    extra$.gene_label <- rownames(extra)
    g <- g + ggrepel::geom_text_repel(
      data = extra, mapping = lab_aes,
      colour = "black", fontface = "bold",
      size = 3.5, box.padding = 0.5, point.padding = 0.3,
      max.overlaps = Inf, force = 5, min.segment.length = 0,
      inherit.aes = FALSE, show.legend = FALSE)
  }

  g
}

#' Combine volcano panels into one row with a shared scale and legend
#'
#' Puts several volcanoes side by side on common limits with one collected
#' legend. Intended for `orientation = "vertical"` panels, where the fold-change
#' axis runs top-to-bottom and reads across a row.
#'
#' @param plots A list of `ggplot` objects from [de_volcano()].
#' @param labels Character vector of panel titles; defaults to `names(plots)`.
#' @param legend_position Where the collected legend goes.
#' @param keep_first_caption Logical. Keep the first panel's caption under that
#'   panel instead of promoting it to a single figure-level caption.
#' @return A `patchwork` object.
#' @export
#' @examples
#' de <- data.frame(
#'   logFC     = c(-3, -1, 0.2, 1.5, 4),
#'   P.Value   = c(1e-6, 0.02, 0.7, 0.01, 1e-8),
#'   adj.P.Val = c(1e-4, 0.06, 0.9, 0.04, 1e-6),
#'   row.names = paste0("Gene", 1:5)
#' )
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   p <- de_volcano(de, fc_cutoff = 1, orientation = "vertical")
#'   de_volcano_grid(list(A = p, B = p))
#' }
de_volcano_grid <- function(plots,
                            labels = names(plots),
                            legend_position = "bottom",
                            keep_first_caption = FALSE) {
  .require_pkg("patchwork", "`de_volcano_grid()`")
  if (!is.list(plots) || !length(plots) ||
      !all(vapply(plots, inherits, logical(1), "ggplot"))) {
    stop("`plots` must be a non-empty list of ggplot objects from ",
         "`de_volcano()`.", call. = FALSE)
  }
  if (is.null(labels)) labels <- rep("", length(plots))

  global_y <- max(vapply(plots, function(p)
    max(abs(layer_scales(p)$y$range$range), na.rm = TRUE), numeric(1)))
  global_x <- max(vapply(plots, function(p)
    max(layer_scales(p)$x$range$range, na.rm = TRUE), numeric(1)))

  # Read the caption off the PLOT, not the build. Under ggplot2 4.0.3
  # `ggplot_build(p)$layout$plot` no longer exists, so the old read returned
  # NULL every time: the loop below then set `caption = NULL` on every panel and
  # the promotion at the end became dead code, silently dropping the caption
  # from both documented modes. The caption is the only place the *realised*
  # raw-p boundary behind the dashed lines is reported, so losing it makes the
  # lines unexplained.
  first_cap <- plots[[1]]$labels$caption
  if (is.null(first_cap) || !nzchar(first_cap)) first_cap <- NULL

  # Each panel already carries its own coord from de_volcano(). Adding a second
  # one would work but makes ggplot2 announce the replacement on every panel,
  # so the shared coord is swapped in place instead -- same result, silent.
  shared_coord <- coord_cartesian(xlim = c(0, global_x),
                                  ylim = c(-global_y, global_y), clip = "off")
  for (i in seq_along(plots)) {
    plots[[i]]$coordinates <- shared_coord
    plots[[i]] <- plots[[i]] +
      ggtitle(labels[i]) +
      labs(caption = if (keep_first_caption && i == 1) first_cap else NULL)
  }

  combined <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = legend_position)

  if (!is.null(first_cap) && !keep_first_caption) {
    combined <- combined & patchwork::plot_annotation(caption = first_cap)
  }
  combined
}
