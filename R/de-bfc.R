#' B-statistic versus log2 fold-change plot
#'
#' limma's `B` is the log-odds that a gene is differentially expressed, so
#' `B > 0` means "more likely DE than not". Plotting it against fold change
#' gives an evidence-versus-effect view that does not depend on a p-value
#' threshold. The y limits always leave room around `B_cutoff`, so the
#' threshold line stays visible even when no gene crosses it -- in that case an
#' italic annotation names the line explicitly.
#'
#' @param de_results Data frame whose rownames are gene IDs, carrying `logFC`
#'   and `B` (limma `topTable()` shape).
#' @param fc_cutoff Numeric absolute log2 fold-change threshold.
#' @param B_cutoff Numeric B-statistic threshold.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#' @param top_n Integer genes labelled per side, ranked by decreasing `B`.
#' @param highlight_gene Character vector of gene IDs always labelled, bold and
#'   black.
#' @param label_method One of `"top"`, `"sig"`, `"B"`, `"log2fc"`; anything
#'   else labels nothing.
#' @param x_breaks Numeric spacing between fold-change axis ticks.
#' @param title Plot title.
#' @param palette Named colours for `NS`, `Log2FC`, `B-statistic` and
#'   `B-statistic & Log2FC`.
#' @param y_padding Numeric head-room kept around `B_cutoff`.
#' @param show_grid Logical. Keep the panel grid.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#' @examples
#' de <- data.frame(
#'   logFC = c(-3, -1, 0.2, 1.5, 4),
#'   B     = c(4, -1, -3, 0.5, 6),
#'   row.names = paste0("Gene", 1:5)
#' )
#' de_bfc_plot(de, fc_cutoff = 1)
de_bfc_plot <- function(
    de_results,
    fc_cutoff      = 2,
    B_cutoff       = 0,
    max.overlaps   = 10,
    top_n          = 5,
    highlight_gene = NULL,
    label_method   = "top",
    x_breaks       = 2,
    title          = "B-statistic vs Log2FC",
    palette        = c(
      "NS"                   = "#7F7F7F",
      "Log2FC"               = "#0173B2",
      "B-statistic"          = "#029E73",
      "B-statistic & Log2FC" = "#D55E00"
    ),
    y_padding      = 1,
    show_grid      = FALSE) {

  missing <- setdiff(c("logFC", "B"), colnames(de_results))
  if (length(missing)) {
    stop("`de_results` is missing column(s): ",
         paste(sprintf("`%s`", missing), collapse = ", "), ".", call. = FALSE)
  }

  df <- dplyr::mutate(
    de_results,
    sig_fc = abs(.data$logFC) >= fc_cutoff,
    sig_B  = .data$B > B_cutoff,
    cat = dplyr::case_when(
      .data$sig_fc & .data$sig_B ~ "B-statistic & Log2FC",
      .data$sig_fc               ~ "Log2FC",
      .data$sig_B                ~ "B-statistic",
      TRUE                       ~ "NS")
  )

  get_top <- function(side) {
    keep <- if (side == "up") df$logFC > 0 else df$logFC < 0
    sub <- df[which(keep), , drop = FALSE]
    utils::head(sub[order(-sub$B), , drop = FALSE], top_n)
  }
  lab_df <- switch(
    label_method,
    top    = rbind(get_top("up"), get_top("down")),
    sig    = df[df$cat == "B-statistic & Log2FC", , drop = FALSE],
    B      = df[df$sig_B, , drop = FALSE],
    log2fc = df[df$sig_fc, , drop = FALSE],
    df[0, , drop = FALSE]
  )
  if (!is.null(highlight_gene)) {
    extra <- df[rownames(df) %in% highlight_gene, , drop = FALSE]
    # `rbind()` on data frames *uniquifies* colliding row names -- "Gene1"
    # becomes "Gene11" -- so appending a highlight that is already in `lab_df`
    # invented a gene that does not exist. The `!duplicated(rownames(...))` that
    # used to follow could never fire, because after the rename the names are no
    # longer duplicates, and the mangled string was then drawn on the figure as
    # a label (and, being unequal to the real id, was not bolded either).
    # Remove the overlap from `lab_df` first so `rbind()` has nothing to rename.
    lab_df <- lab_df[!rownames(lab_df) %in% rownames(extra), , drop = FALSE]
    lab_df <- rbind(lab_df, extra)
  }

  xmax <- ceiling(max(abs(df$logFC), na.rm = TRUE) / x_breaks) * x_breaks
  ymin <- min(floor(min(df$B, na.rm = TRUE)), B_cutoff - y_padding)
  ymax <- max(ceiling(max(df$B, na.rm = TRUE)), B_cutoff + y_padding)
  if (ymax <= B_cutoff) ymax <- B_cutoff + y_padding
  none_above <- max(df$B, na.rm = TRUE) < B_cutoff
  if (none_above) ymax <- B_cutoff + y_padding * 2

  dark_pal <- .de_shade(palette)

  g <- ggplot(df, aes(x = .data$logFC, y = .data$B, colour = .data$cat)) +
    geom_point(size = 2, alpha = 0.65) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = B_cutoff, linetype = "dashed") +
    scale_colour_manual(
      name = NULL, values = palette, breaks = names(palette),
      labels = c(
        "B-statistic & Log2FC" = sprintf("B > %.1f & |log2FC| \u2265 %.1f",
                                         B_cutoff, fc_cutoff),
        "Log2FC"               = sprintf("|log2FC| \u2265 %.1f", fc_cutoff),
        "B-statistic"          = sprintf("B > %.1f", B_cutoff),
        "NS"                   = "NS")) +
    scale_x_continuous(breaks = seq(-xmax, xmax, by = x_breaks),
                       limits = c(-xmax, xmax)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "log2(FC)", y = "B statistic (log-odds of DE)", title = title,
         caption = sprintf(
           "Dashed lines: horiz. \u2013 B > %.1f; vert. \u2013 |log2FC| \u2265 %.1f",
           B_cutoff, fc_cutoff))

  if (none_above) {
    g <- g + annotate("text", x = xmax * 0.85,
                      y = B_cutoff + y_padding * 0.5,
                      label = sprintf("B = %.1f threshold", B_cutoff),
                      colour = "black", size = 3.5, fontface = "italic")
  }

  g <- g + .de_theme()
  if (!show_grid) {
    g <- g + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  }

  if (nrow(lab_df)) {
    lab_df$.gene_label <- rownames(lab_df)
    g <- g + ggrepel::geom_text_repel(
      data = lab_df,
      mapping = aes(x = .data$logFC, y = .data$B, label = .data$.gene_label),
      colour = ifelse(rownames(lab_df) %in% highlight_gene, "black",
                      unname(dark_pal[lab_df$cat])),
      fontface = ifelse(rownames(lab_df) %in% highlight_gene, "bold", "plain"),
      size = 3.5, box.padding = 0.4, point.padding = 0.3,
      max.overlaps = max.overlaps, min.segment.length = 0,
      inherit.aes = FALSE, show.legend = FALSE)
  }

  g
}
