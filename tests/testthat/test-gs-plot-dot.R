test_that("gs_plot_dot returns a ggplot with one point per selected pathway", {
  res <- fake_plot_result(n = 8L)
  p <- gs_plot_dot(res, top_n = 5)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(layer_data_for(p)), 5L)
})

test_that("the x axis label comes from stat_type, never a literal NES", {
  expect_equal(
    gs_plot_dot(fake_plot_result())$labels$x, "NES"
  )
  expect_equal(
    gs_plot_dot(fake_plot_result(stat_type = "t"))$labels$x, "t statistic"
  )
  expect_equal(
    gs_plot_dot(
      fake_plot_result(stat_type = "log2_fold_enrichment")
    )$labels$x,
    "log2 fold enrichment"
  )
})

test_that("selection and highlighting are independent", {
  res <- fake_plot_result(n = 8L)
  # Every selected pathway is drawn ...
  p <- gs_plot_dot(res, top_n = 8, highlight = 0.05)
  expect_equal(nrow(layer_data_for(p, 1L)), 8L)
  # ... but only the significant ones get the outline layer.
  n_sig <- sum(res$padj < 0.05)
  expect_equal(nrow(layer_data_for(p, 2L)), n_sig)

  # highlight = NULL drops the outline layer entirely
  p2 <- gs_plot_dot(res, top_n = 8, highlight = NULL)
  expect_length(ggplot2::ggplot_build(p2)$data, 1L)
})

test_that("base dots use a transparent outline, not NA", {
  d <- layer_data_for(gs_plot_dot(fake_plot_result()))
  expect_identical(unique(d$colour), "transparent")
  expect_false(anyNA(d$colour))
  expect_equal(unique(d$shape), 21)
})

test_that("direction = 'up' / 'down' filter the rows drawn", {
  res <- fake_plot_result(n = 8L)
  up <- gs_plot_dot(res, top_n = 20, direction = "up")
  expect_true(all(attr(up, "gs_source")$stat > 0))
  down <- gs_plot_dot(res, top_n = 20, direction = "down")
  expect_true(all(attr(down, "gs_source")$stat < 0))
})

test_that("facet = 'direction' selects top N within each facet", {
  res <- fake_plot_result(n = 10L)
  src <- attr(
    gs_plot_dot(res, top_n = 3, facet = "direction"), "gs_source"
  )
  expect_setequal(as.character(unique(src$.facet_row)), c("Up", "Down"))
  expect_true(all(table(src$.facet_row) <= 3L))
})

test_that("facet strips show database_label, never the snake_case key", {
  res <- fake_plot_result(n = 4L, databases = c("msigdb_H", "mitopathways"))
  src <- attr(
    gs_plot_dot(res, top_n = 4, facet = "database"), "gs_source"
  )
  expect_setequal(
    as.character(unique(src$.facet_row)),
    c("Label for msigdb_H", "Label for mitopathways")
  )
  expect_false(any(src$.facet_row %in% c("msigdb_H", "mitopathways")))
})

test_that("database_labels overrides the attribute", {
  res <- fake_plot_result(n = 3L)
  src <- attr(
    gs_plot_dot(res, facet = "database",
                database_labels = c(msigdb_H = "Hallmark")),
    "gs_source"
  )
  expect_identical(unique(as.character(src$.facet_row)), "Hallmark")
})

test_that("an unmapped database key falls back to the key itself", {
  res <- fake_plot_result(n = 2L)
  attr(res, "database_label") <- c(other = "Other")
  src <- attr(gs_plot_dot(res, facet = "database"), "gs_source")
  expect_identical(unique(as.character(src$.facet_row)), "msigdb_H")
})

test_that("compare puts every panel on the same pathways", {
  res <- fake_plot_result(n = 4L, contrasts = c("A-B", "C-D"))
  src <- attr(
    gs_plot_dot(res, top_n = 3, compare = "contrast"), "gs_source"
  )
  by_panel <- split(src$pathway_id, src$.facet_col)
  expect_length(by_panel, 2L)
  expect_setequal(by_panel[[1L]], by_panel[[2L]])
})

test_that("compare rejects an unsupported column", {
  expect_error(gs_plot_dot(fake_plot_result(), compare = "method"),
               "`compare`")
})

test_that("aes_x = 'gene_ratio' uses the leading edge, and errors without it", {
  res <- fake_plot_result(n = 4L)
  src <- attr(gs_plot_dot(res, aes_x = "gene_ratio"), "gs_source")
  expect_equal(src$gene_ratio, lengths(res$leading_edge[
    match(src$pathway_id, res$pathway_id)
  ]) / res$n_genes[match(src$pathway_id, res$pathway_id)])

  bare <- fake_plot_result(n = 4L, leading_edge = FALSE)
  expect_error(gs_plot_dot(bare, aes_x = "gene_ratio"), "gene_ratio")
})

test_that("fill limits are symmetric and out-of-range values are squished", {
  res <- fake_plot_result(n = 6L)
  d <- layer_data_for(gs_plot_dot(res, top_n = 6, limits = c(-1, 1)))
  expect_false(anyNA(d$fill))
})

test_that("label = TRUE adds a ggrepel layer", {
  p <- gs_plot_dot(fake_plot_result(n = 3L), label = TRUE)
  classes <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomTextRepel" %in% classes)
})

test_that("a zero-row gs_result yields an empty plot, not an error", {
  res <- fake_plot_result(n = 3L)[0, ]
  p <- gs_plot_dot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$subtitle, "No pathways to plot")
})

test_that("non-gs_result input is rejected by name", {
  expect_error(gs_plot_dot(mtcars), "`x` must be a gs_result")
})

test_that("top_n must be a positive number", {
  expect_error(gs_plot_dot(fake_plot_result(), top_n = 0), "`top_n`")
})

test_that("the source table travels with the plot", {
  p <- gs_plot_dot(fake_plot_result(n = 5L), top_n = 4)
  src <- attr(p, "gs_source")
  expect_s3_class(src, "data.frame")
  expect_equal(nrow(src), 4L)
  expect_true(all(c("pathway_id", "label", "stat", "padj", "significant")
                  %in% names(src)))
})

# --- regression: padj == 0 -----------------------------------------------------
# `-log10(0)` is `Inf`, which no size scale can map, so ggplot2 dropped the point
# with no warning -- silently deleting the strongest hit in the analysis from the
# figure. Before the clamp in `.gs_plot_frame()` the built layer's `size` column
# came back `Inf, 10, 2` and only two of three points were drawn.
test_that("a pathway with padj == 0 is still drawn, at the largest size", {
  r <- fake_plot_result(n = 3L)
  r$padj <- c(0, 1e-3, 0.2)
  p <- gs_plot_dot(r)
  d <- layer_data_for(p)
  expect_equal(nrow(d), 3L)
  expect_true(all(is.finite(d$size)))
  # The clamped point must remain the most significant, not merely finite.
  src <- attr(p, "gs_source")
  expect_equal(src$pathway_id[which.max(src$neg_log_padj)],
               r$pathway_id[which(r$padj == 0)])
})
