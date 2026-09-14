#' Mean-difference (MD) plot for differential-expression results
#'
#' `limma::plotMD()`'s idea with the additions that make a systematic shift
#' visible at a glance: a light grey cloud for non-significant genes so the
#' coloured significant ones stay readable, a dashed LOESS trend exposing bias
#' across the expression range, a dotted guide at median expression that
#' catches normalisation problems, and optional per-quadrant counts.
#'
#' Significance is called on `adj.P.Val` alone; `fc_cutoff` only draws guides
#' and feeds `label_method = "sig"`/`"log2fc"`.
#'
#' @param fit An `MArrayLM` object from limma.
#' @param coef Integer index or character name of the coefficient to plot.
#' @param de_results Optional `topTable()` data frame matching `fit`. When
#'   `NULL` it is extracted from `fit` with `sort.by = "none"`.
#' @param fc_cutoff Numeric absolute log2 fold-change guide.
#' @param fdr_cutoff Numeric FDR threshold for the Up/Down call.
#' @param top_n Integer genes labelled per direction.
#' @param highlight_gene Character vector of gene IDs always labelled, bold and
#'   black.
#' @param label_method One of `"top"`, `"sig"`, `"fdr"`, `"log2fc"`, `"none"`.
#' @param max.overlaps Passed to [ggrepel::geom_text_repel()].
#' @param title Plot title; defaults to `"MD plot: <coef>"`.
#' @param palette Named colours for `Up`, `Down` and `NS`.
#' @param show_grid Logical. Keep the panel grid.
#' @param show_quadrant_counts Logical. Annotate the number of significant
#'   genes in each quadrant of (median expression, zero fold change).
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   y <- matrix(rnorm(200), nrow = 50,
#'               dimnames = list(paste0("G", 1:50), paste0("S", 1:4)))
#'   design <- cbind(Intercept = 1, Group = c(0, 0, 1, 1))
#'   fit <- limma::eBayes(limma::lmFit(y, design))
#'   de_md_plot(fit, coef = "Group")
#' }
de_md_plot <- function(
    fit,
    coef,
    de_results     = NULL,
    fc_cutoff      = 1,
    fdr_cutoff     = 0.05,
    top_n          = 5,
    highlight_gene = NULL,
    label_method   = "top",
    max.overlaps   = 10,
    title          = NULL,
    palette        = c(Up = "#D55E00", Down = "#0072B2", NS = "#999999"),
    show_grid      = FALSE,
    show_quadrant_counts = TRUE) {

  coef_name <- if (is.numeric(coef)) colnames(fit)[coef] else as.character(coef)
  if (is.null(title)) title <- paste("MD plot:", coef_name)

  if (is.null(de_results)) {
    .require_pkg("limma", "`de_md_plot(de_results = NULL)`",
                'BiocManager::install("limma")')
    de_results <- limma::topTable(fit, coef = coef, number = Inf,
                                  sort.by = "none")
  }
  if (!"AveExpr" %in% colnames(de_results) && !is.null(fit$Amean)) {
    de_results$AveExpr <- fit$Amean
  }

  req <- c("logFC", "AveExpr", "adj.P.Val")
  missing <- setdiff(req, colnames(de_results))
  if (length(missing)) {
    stop("`de_results` is missing column(s): ",
         paste(sprintf("`%s`", missing), collapse = ", "), ".", call. = FALSE)
  }

  df <- dplyr::mutate(
    de_results,
    sig_fc  = abs(.data$logFC) >= fc_cutoff,
    sig_fdr = .data$adj.P.Val  <= fdr_cutoff,
    status  = dplyr::case_when(
      .data$sig_fdr & .data$logFC > 0 ~ "Up",
      .data$sig_fdr & .data$logFC < 0 ~ "Down",
      TRUE                            ~ "NS"),
    gene_id = rownames(de_results)
  )

  get_top <- function(direction) {
    keep <- if (direction == "up") df$logFC > 0 else df$logFC < 0
    sub <- df[which(keep), , drop = FALSE]
    utils::head(sub[order(sub$adj.P.Val), , drop = FALSE], top_n)
  }
  lab_df <- switch(
    label_method,
    top    = rbind(get_top("up"), get_top("down")),
    sig    = df[df$sig_fc & df$sig_fdr, , drop = FALSE],
    fdr    = df[df$sig_fdr, , drop = FALSE],
    log2fc = df[df$sig_fc, , drop = FALSE],
    df[0, , drop = FALSE]
  )
  if (!is.null(highlight_gene)) {
    lab_df <- rbind(lab_df, df[df$gene_id %in% highlight_gene, , drop = FALSE])
    lab_df <- lab_df[!duplicated(lab_df$gene_id), , drop = FALSE]
  }

  df_ns  <- df[df$status == "NS", , drop = FALSE]
  df_sig <- df[df$status != "NS", , drop = FALSE]

  x_pad <- diff(range(df$AveExpr)) * 0.05
  y_pad <- diff(range(df$logFC))  * 0.05
  dark_pal <- .de_shade(palette)

  g <- ggplot() +
    geom_point(data = df_ns, aes(x = .data$AveExpr, y = .data$logFC),
               colour = palette[["NS"]], size = 1.3, alpha = 0.25,
               shape = 16) +
    geom_point(data = df_sig,
               aes(x = .data$AveExpr, y = .data$logFC, colour = .data$status),
               size = 2.2, alpha = 0.9, shape = 16) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_smooth(data = df, aes(x = .data$AveExpr, y = .data$logFC),
                method = "loess", se = FALSE, colour = "red",
                linetype = 2, linewidth = 0.6) +
    geom_vline(xintercept = stats::median(df$AveExpr), colour = "grey60",
               linetype = 3) +
    scale_colour_manual(values = palette, name = "Significance",
      breaks = c("Up", "Down"),
      labels = c(sprintf("Up  (FDR \u2264 %.2g)", fdr_cutoff),
                 sprintf("Down (FDR \u2264 %.2g)", fdr_cutoff))) +
    coord_cartesian(xlim = range(df$AveExpr) + c(-x_pad, x_pad),
                    ylim = range(df$logFC)  + c(-y_pad, y_pad)) +
    labs(x = "Average Expression", y = "log2 Fold Change", title = title,
         caption = sprintf("Dashed: |logFC| \u2265 %.1f", fc_cutoff)) +
    .de_theme()

  if (!show_grid) {
    g <- g + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  }

  if (show_quadrant_counts && nrow(df_sig)) {
    qx <- stats::median(df$AveExpr)
    quad <- data.frame(
      x     = c(min(df$AveExpr), min(df$AveExpr),
                max(df$AveExpr), max(df$AveExpr)),
      y     = c(max(df$logFC), min(df$logFC), max(df$logFC), min(df$logFC)),
      hjust = c(0, 0, 1, 1),
      vjust = c(1, 0, 1, 0),
      lbl   = c(sum(df_sig$logFC > 0 & df_sig$AveExpr < qx),
                sum(df_sig$logFC < 0 & df_sig$AveExpr < qx),
                sum(df_sig$logFC > 0 & df_sig$AveExpr > qx),
                sum(df_sig$logFC < 0 & df_sig$AveExpr > qx))
    )
    g <- g + geom_text(
      data = quad,
      mapping = aes(x = .data$x, y = .data$y, label = .data$lbl,
                    hjust = .data$hjust, vjust = .data$vjust),
      colour = "grey30", size = 3)
  }

  if (nrow(lab_df)) {
    g <- g + ggrepel::geom_text_repel(
      data = lab_df,
      mapping = aes(x = .data$AveExpr, y = .data$logFC, label = .data$gene_id),
      colour = ifelse(lab_df$gene_id %in% highlight_gene, "black",
                      unname(dark_pal[lab_df$status])),
      fontface = ifelse(lab_df$gene_id %in% highlight_gene, "bold", "plain"),
      size = 3.3, box.padding = 0.35, point.padding = 0.25,
      max.overlaps = max.overlaps, inherit.aes = FALSE, show.legend = FALSE)
  }

  g
}
