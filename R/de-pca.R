#' Principal components of a DGEList
#'
#' Shared by [de_pca()] and [de_pca_3d()]. Computes logCPM, drops zero-variance
#' rows and columns (PCA is undefined on them), runs [stats::prcomp()] and
#' returns the scores joined to the sample metadata.
#'
#' @param dge A `DGEList`.
#' @param min_dim Minimum number of surviving rows/columns required.
#' @return A list with `scores` (data frame of PCs plus `dge$samples`) and
#'   `percent_var` (numeric, percent variance per PC).
#' @keywords internal
.de_pca_scores <- function(dge, min_dim = 2L) {
  if (!inherits(dge, "DGEList")) {
    stop("`dge` must be a `DGEList`; build one with `build_dge()`.",
         call. = FALSE)
  }
  .require_pkg("edgeR", "PCA of a `DGEList`", 'BiocManager::install("edgeR")')

  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  keep_r <- apply(logcpm, 1, stats::var) > 0
  keep_c <- apply(logcpm, 2, stats::var) > 0
  logcpm <- logcpm[keep_r, keep_c, drop = FALSE]
  if (nrow(logcpm) < min_dim || ncol(logcpm) < min_dim) {
    stop(sprintf(
      "`dge` has too little variation for PCA: %d gene(s) and %d sample(s) survive the ",
      nrow(logcpm), ncol(logcpm)),
      "zero-variance filter, need at least ", min_dim, " of each.",
      call. = FALSE)
  }

  pca <- stats::prcomp(t(logcpm))
  n_pc <- min(ncol(pca$x), 3L)
  scores <- data.frame(pca$x[, seq_len(n_pc), drop = FALSE],
                       dge$samples[colnames(logcpm), , drop = FALSE],
                       check.names = FALSE)
  scores$sample_id <- rownames(scores)
  list(scores = scores,
       percent_var = pca$sdev^2 / sum(pca$sdev^2) * 100)
}

#' Two-dimensional PCA of a DGEList
#'
#' PCA on logCPM, plotted on a fixed 1:1 aspect ratio so distances read
#' correctly. Grouping is by metadata column name -- nothing about organs,
#' tissues or study design is hard-coded.
#'
#' @param dge A `DGEList`, e.g. from [build_dge()].
#' @param colour_by Name of a `dge$samples` column mapped to point colour, or
#'   `NULL` for one colour.
#' @param shape_by Name of a `dge$samples` column mapped to point shape, or
#'   `NULL`.
#' @param label Logical. Draw the `colour_by` value above each point.
#' @param title Plot title.
#' @param xlim_abs Optional symmetric x-axis limits, applied as a zoom
#'   (`coord_fixed()`), so a sample outside them is clipped from view rather
#'   than dropped from the figure. `NULL` derives them from the data with 10%
#'   head-room.
#' @param ylim_abs Optional symmetric y-axis limits, with the same zoom
#'   behaviour as `xlim_abs`.
#' @param point_size Numeric point size.
#' @param label_size Numeric label size.
#' @param palette Character vector of colours recycled across the
#'   `colour_by` levels.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom rlang .data
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(rpois(400, 50), nrow = 100,
#'     dimnames = list(paste0("G", 1:100), paste0("S", 1:4)))
#'   samples <- data.frame(group = c("WT", "WT", "KO", "KO"),
#'                         row.names = colnames(counts))
#'   dge <- build_dge(counts, samples, data.frame(gene = rownames(counts)))
#'   de_pca(dge, colour_by = "group")
#' }
de_pca <- function(dge,
                   colour_by  = "group",
                   shape_by   = NULL,
                   label      = TRUE,
                   title      = "PCA Plot",
                   xlim_abs   = NULL,
                   ylim_abs   = NULL,
                   point_size = 5,
                   label_size = 4,
                   palette = c("#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00",
                               "#CC79A7", "#000000")) {
  pca <- .de_pca_scores(dge, min_dim = 2L)
  d   <- pca$scores
  pv  <- pca$percent_var

  for (col in c(colour_by, shape_by)) {
    if (!is.null(col) && !col %in% colnames(d)) {
      stop(sprintf("`%s` is not a column of `dge$samples`. Available: %s.",
                   col, paste(colnames(dge$samples), collapse = ", ")),
           call. = FALSE)
    }
  }

  mapping <- aes(x = .data$PC1, y = .data$PC2)
  if (!is.null(colour_by)) mapping$colour <- rlang::sym(colour_by)
  if (!is.null(shape_by))  mapping$shape  <- rlang::sym(shape_by)

  xlim_val <- xlim_abs %||% (max(abs(d$PC1)) * 1.1)
  ylim_val <- ylim_abs %||% (max(abs(d$PC2)) * 1.1)

  g <- ggplot(d, mapping) +
    geom_point(size = point_size) +
    labs(colour = colour_by, shape = shape_by,
         x = sprintf("PC1: %.1f%% variance", pv[1]),
         y = sprintf("PC2: %.1f%% variance", pv[2]),
         title = title)

  if (label && !is.null(colour_by)) {
    d$.label <- as.character(d[[colour_by]])
    g <- g + geom_text(data = d, aes(x = .data$PC1, y = .data$PC2,
                                     label = .data$.label),
                       vjust = -1.5, hjust = 0.5, size = label_size,
                       fontface = "bold", inherit.aes = FALSE,
                       show.legend = FALSE)
  }
  if (!is.null(colour_by)) {
    n_lev <- length(unique(d[[colour_by]]))
    g <- g + scale_colour_manual(
      values = rep_len(palette, max(n_lev, 1L)))
  }

  # `xlim()`/`ylim()` set *scale* limits, which convert out-of-range samples to
  # NA and drop them -- the only signal being a "Removed n rows" warning at
  # print time. Pinning `xlim_abs` to make panels comparable across organs
  # therefore deleted the outlier sample instead of zooming past it. Every
  # sibling renderer clips with coord_*, and this one already called
  # `coord_fixed()`, so the limits belong there.
  g +
    .de_theme() +
    theme(legend.position = "right",
          legend.box = "vertical",
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10)) +
    coord_fixed(xlim = c(-xlim_val, xlim_val),
                ylim = c(-ylim_val, ylim_val))
}

#' Interactive three-dimensional PCA of a DGEList
#'
#' The 3D counterpart of [de_pca()], rendered with plotly. Needs at least three
#' principal components, so at least three samples and three variable genes.
#'
#' @param dge A `DGEList`, e.g. from [build_dge()].
#' @param colour_by Name of a `dge$samples` column mapped to marker colour, or
#'   `NULL`.
#' @param symbol_by Name of a `dge$samples` column mapped to marker symbol, or
#'   `NULL`.
#' @param title Plot title.
#' @param point_size Marker size.
#'
#' @return A `plotly` object.
#' @export
#' @examples
#' if (requireNamespace("edgeR", quietly = TRUE) &&
#'     requireNamespace("plotly", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(rpois(600, 50), nrow = 100,
#'     dimnames = list(paste0("G", 1:100), paste0("S", 1:6)))
#'   samples <- data.frame(group = rep(c("WT", "KO"), each = 3),
#'                         row.names = colnames(counts))
#'   dge <- build_dge(counts, samples, data.frame(gene = rownames(counts)))
#'   de_pca_3d(dge, colour_by = "group")
#' }
de_pca_3d <- function(dge,
                      colour_by  = "group",
                      symbol_by  = NULL,
                      title      = "3D PCA Plot",
                      point_size = 8) {
  .require_pkg("plotly", "`de_pca_3d()`")
  pca <- .de_pca_scores(dge, min_dim = 3L)
  d   <- pca$scores
  if (!"PC3" %in% colnames(d)) {
    stop("PCA produced fewer than 3 components; `de_pca_3d()` needs 3. ",
         "Use `de_pca()` instead.", call. = FALSE)
  }
  pv <- round(pca$percent_var, 1)

  for (col in c(colour_by, symbol_by)) {
    if (!is.null(col) && !col %in% colnames(d)) {
      stop(sprintf("`%s` is not a column of `dge$samples`. Available: %s.",
                   col, paste(colnames(dge$samples), collapse = ", ")),
           call. = FALSE)
    }
  }

  hover <- paste0("Sample: ", d$sample_id)
  if (!is.null(colour_by)) {
    hover <- paste0(hover, "<br>", colour_by, ": ", as.character(d[[colour_by]]))
  }

  args <- list(
    data = d, x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers",
    text = hover, hoverinfo = "text",
    marker = list(size = point_size)
  )
  if (!is.null(colour_by)) args$color  <- stats::as.formula(paste0("~", colour_by))
  if (!is.null(symbol_by)) args$symbol <- stats::as.formula(paste0("~", symbol_by))

  p <- do.call(plotly::plot_ly, args)
  plotly::layout(
    p,
    title = list(text = title, x = 0.5, xanchor = "center"),
    scene = list(
      xaxis = list(title = sprintf("PC1: %.1f%% variance", pv[1])),
      yaxis = list(title = sprintf("PC2: %.1f%% variance", pv[2])),
      zaxis = list(title = sprintf("PC3: %.1f%% variance", pv[3]))
    ),
    legend = list(orientation = "v", x = 1.05, y = 0.9,
                  title = list(text = colour_by %||% ""))
  )
}
