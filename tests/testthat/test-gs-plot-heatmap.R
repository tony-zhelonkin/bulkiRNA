test_that("a gs_result heatmap tiles pathways by contrast", {
  res <- fake_plot_result(n = 4L, contrasts = c("A-B", "C-D"))
  p <- gs_plot_heatmap(res, top_n = 3)
  expect_s3_class(p, "ggplot")
  src <- attr(p, "gs_source")
  expect_setequal(unique(src$column), c("A-B", "C-D"))
  expect_equal(length(unique(src$pathway_id)), 3L)
})

test_that("by = 'database' shows labels, not keys", {
  res <- fake_plot_result(n = 3L, databases = c("msigdb_H", "mitopathways"))
  src <- attr(
    gs_plot_heatmap(res, top_n = 4, by = "database"), "gs_source"
  )
  expect_true(all(grepl("^Label for ", src$column)))
})

test_that("the fill legend is named from stat_type", {
  p <- gs_plot_heatmap(fake_plot_result(stat_type = "t"))
  expect_equal(p$scales$get_scales("fill")$name, "t statistic")
})

test_that("significant tiles get an asterisk layer", {
  res <- fake_plot_result(n = 6L)
  p <- gs_plot_heatmap(res, top_n = 6, highlight = 0.05)
  expect_length(ggplot2::ggplot_build(p)$data, 2L)
  expect_equal(nrow(layer_data_for(p, 2L)), sum(res$padj < 0.05))

  p2 <- gs_plot_heatmap(res, top_n = 6, highlight = NULL)
  expect_length(ggplot2::ggplot_build(p2)$data, 1L)
})

test_that("a gs_matrix heatmap tiles pathways by sample", {
  m <- fake_plot_matrix(n_path = 5L, n_samp = 6L)
  p <- gs_plot_heatmap(m, top_n = 3)
  src <- attr(p, "gs_source")
  expect_equal(length(unique(src$sample)), 6L)
  expect_equal(length(unique(src$pathway_id)), 3L)
  expect_equal(p$scales$get_scales("fill")$name, "gsva score")
})

test_that("samples select and order the columns", {
  m <- fake_plot_matrix(n_samp = 6L)
  src <- attr(gs_plot_heatmap(m, samples = c("s3", "s1")), "gs_source")
  expect_equal(levels(src$sample), c("s3", "s1"))
  expect_error(gs_plot_heatmap(m, samples = "nope"), "`samples`")
})

test_that("group facets the columns by sample metadata", {
  m <- fake_plot_matrix()
  src <- attr(gs_plot_heatmap(m, group = "group"), "gs_source")
  expect_setequal(unique(src$.facet_col), c("WT", "KO"))
  expect_error(gs_plot_heatmap(m, group = "missing"), "sample_data")
})

test_that("an unsupported class is rejected", {
  expect_error(gs_plot_heatmap(mtcars), "gs_result or a gs_matrix")
})

test_that("an empty result gives an empty plot", {
  expect_equal(
    gs_plot_heatmap(fake_plot_result(n = 3L)[0, ])$labels$subtitle,
    "No pathways to plot"
  )
})

# --- regression: the literal "contrast" placeholder ----------------------------
# gs_test()'s `contrast` formal defaults to the string "contrast", so a
# single-contrast result printed an axis tick reading the word "contrast".
test_that("the contrast axis blanks gs_test()'s placeholder but keeps real names", {
  r <- fake_plot_result(n = 3L)
  r$contrast <- "contrast"
  p <- gs_plot_heatmap(r)
  expect_true(all(as.character(p$data$column) == ""))

  r2 <- fake_plot_result(n = 3L, contrasts = "KO-WT")
  p2 <- gs_plot_heatmap(r2)
  expect_true(all(as.character(p2$data$column) == "KO-WT"))
})
