test_that("de_pca returns a ggplot with variance-labelled axes", {
  dge <- fake_dge()
  p <- de_pca(dge, colour_by = "group", shape_by = "organ")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "^PC1: [0-9.]+% variance$")
  expect_match(p$labels$y, "^PC2: [0-9.]+% variance$")
})

test_that("an unknown grouping column errors and lists the real ones", {
  dge <- fake_dge()
  expect_error(de_pca(dge, colour_by = "tissue"), "not a column")
  expect_error(de_pca(dge, colour_by = "tissue"), "group")
})

test_that("de_pca rejects anything that is not a DGEList", {
  expect_error(de_pca(matrix(1:4, 2)), "must be a `DGEList`")
})

test_that("de_pca refuses data with no variation", {
  skip_if_not_installed("edgeR")
  counts <- matrix(10L, nrow = 5, ncol = 4,
                   dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))
  samples <- data.frame(group = c("a", "a", "b", "b"),
                        row.names = colnames(counts))
  dge <- build_dge(counts, samples, data.frame(gene = rownames(counts)))
  expect_error(de_pca(dge), "`dge` has too little variation for PCA", fixed = TRUE)
})

test_that("de_pca_3d returns a plotly object", {
  skip_if_not_installed("plotly")
  dge <- fake_dge(n_sample = 6)
  p <- de_pca_3d(dge, colour_by = "group")
  expect_s3_class(p, "plotly")
})

# --- regression from the DE review -------------------------------------------

test_that("xlim_abs zooms instead of deleting samples", {
  skip_if_not_installed("edgeR")
  # `xlim()`/`ylim()` set scale limits, so a sample outside them became NA and
  # was dropped -- the only signal a "Removed n rows" warning at print time.
  dge <- fake_dge()
  full <- de_pca(dge, colour_by = "group")
  n <- nrow(ggplot2::ggplot_build(full)$data[[1]])

  tight <- de_pca(dge, colour_by = "group", xlim_abs = 1e-3, ylim_abs = 1e-3)
  pts <- ggplot2::ggplot_build(tight)$data[[1]]
  expect_identical(nrow(pts), n)
  expect_false(any(is.na(pts$x)))
})
