# Builders for the de_* / IO tests.

# A small limma-shaped topTable. `n_sig` genes are FDR-significant, and their
# raw p-values are strictly smaller than every non-significant one, so the
# decision-by-FDR boundary is unambiguous and checkable by hand.
fake_de_table <- function(n = 20, n_sig = 5, seed = 1) {
  set.seed(seed)
  p <- c(10^seq(-8, -4, length.out = n_sig),
         seq(0.2, 0.9, length.out = n - n_sig))
  adj <- c(rep(0.01, n_sig), rep(0.5, n - n_sig))
  data.frame(
    logFC     = seq(-4, 4, length.out = n),
    AveExpr   = seq(1, 10, length.out = n),
    t         = seq(-6, 6, length.out = n),
    P.Value   = p,
    adj.P.Val = adj,
    B         = seq(6, -6, length.out = n),
    row.names = paste0("Gene", seq_len(n))
  )
}

# The y-intercept of every dashed horizontal line in a built plot.
hline_yintercepts <- function(p) {
  b <- ggplot2::ggplot_build(p)
  out <- numeric(0)
  for (i in seq_along(b$plot$layers)) {
    l <- b$plot$layers[[i]]
    if (inherits(l$geom, "GeomHline")) out <- c(out, b$data[[i]]$yintercept)
  }
  unique(out)
}

vline_xintercepts <- function(p) {
  b <- ggplot2::ggplot_build(p)
  out <- numeric(0)
  for (i in seq_along(b$plot$layers)) {
    l <- b$plot$layers[[i]]
    if (inherits(l$geom, "GeomVline")) out <- c(out, b$data[[i]]$xintercept)
  }
  unique(out)
}

fake_counts <- function(n_gene = 40, n_sample = 6, seed = 2) {
  set.seed(seed)
  m <- matrix(rpois(n_gene * n_sample, lambda = 60), nrow = n_gene,
              dimnames = list(paste0("G", seq_len(n_gene)),
                              paste0("S", seq_len(n_sample))))
  m
}

fake_dge <- function(...) {
  testthat::skip_if_not_installed("edgeR")
  counts <- fake_counts(...)
  samples <- data.frame(
    group = rep(c("WT", "KO"), length.out = ncol(counts)),
    organ = rep(c("Lung", "Spleen"), length.out = ncol(counts)),
    row.names = colnames(counts)
  )
  build_dge(counts, samples, data.frame(gene = rownames(counts)))
}
