# Builders for the renderer tests. Richer than helper-gs.R's fake_gs_result():
# these carry leading edges, several databases and several contrasts, which is
# what the facet / compare / heatmap paths need.

fake_plot_result <- function(n = 6L,
                             databases = "msigdb_H",
                             contrasts = "KO-WT",
                             stat_type = "NES",
                             leading_edge = TRUE) {
  grid <- expand.grid(
    i = seq_len(n), database = databases, contrast = contrasts,
    stringsAsFactors = FALSE
  )
  df <- data.frame(
    pathway_id = paste0("HALLMARK_SET_", grid$i),
    pathway_name = paste0("HALLMARK_SET_", grid$i),
    database = grid$database,
    contrast = grid$contrast,
    method = "fgsea",
    n_genes = 20L + grid$i,
    n_genes_tested = 15L + grid$i,
    # Indexed by `i`, so every database and contrast spans both directions.
    stat = seq(-2.5, 2.5, length.out = n)[grid$i],
    stat_type = stat_type,
    p_value = seq(1e-6, 0.4, length.out = n)[grid$i],
    padj = seq(1e-4, 0.6, length.out = n)[grid$i],
    stringsAsFactors = FALSE
  )
  if (leading_edge) {
    df$leading_edge <- lapply(seq_len(nrow(df)), function(i) {
      paste0("G", seq_len(1 + (i %% 5)))
    })
  }
  res <- bulkiRNA:::gs_result(df)
  attr(res, "database_label") <- stats::setNames(
    paste("Label for", databases), databases
  )
  res
}

fake_plot_matrix <- function(n_path = 5L, n_samp = 6L) {
  set.seed(42)
  m <- matrix(
    stats::rnorm(n_path * n_samp),
    nrow = n_path,
    dimnames = list(paste0("HALLMARK_SET_", seq_len(n_path)),
                    paste0("s", seq_len(n_samp)))
  )
  bulkiRNA:::gs_matrix(
    m,
    database = "msigdb_H", method = "gsva", score_type = "gsva",
    pathway_names = stats::setNames(rownames(m), rownames(m)),
    sample_data = data.frame(
      group = rep(c("WT", "KO"), length.out = n_samp),
      row.names = colnames(m)
    )
  )
}

layer_data_for <- function(p, i = 1L) ggplot2::ggplot_build(p)$data[[i]]

# A fresh temp directory per call, so file-writing tests never collide.
tmp_dir <- function() {
  d <- tempfile("bulkirna-test-")
  dir.create(d, recursive = TRUE)
  d
}
