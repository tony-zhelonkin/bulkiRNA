test_that("gs_plot_bar draws one bar per selected pathway", {
  res <- fake_plot_result(n = 8L)
  p <- gs_plot_bar(res, top_n = 4)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(layer_data_for(p)), 4L)
})

test_that("the value axis is labelled from stat_type", {
  expect_equal(gs_plot_bar(fake_plot_result())$labels$y, "NES")
  expect_equal(
    gs_plot_bar(fake_plot_result(stat_type = "t"))$labels$y, "t statistic"
  )
})

test_that("non-significant pathways are shown, not dropped", {
  res <- fake_plot_result(n = 8L)
  src <- attr(gs_plot_bar(res, top_n = 8), "gs_source")
  expect_equal(nrow(src), 8L)
  expect_true(any(!src$significant))
})

test_that("significant bars are outlined and others are transparent", {
  src <- attr(gs_plot_bar(fake_plot_result(n = 8L), top_n = 8), "gs_source")
  expect_setequal(unique(src$outline), c("black", "transparent"))
  expect_true(all(src$outline[src$significant] == "black"))
})

test_that("padj_max applies a hard filter before selection", {
  res <- fake_plot_result(n = 8L)
  src <- attr(gs_plot_bar(res, top_n = 8, padj_max = 0.05), "gs_source")
  expect_true(all(src$padj < 0.05))
  expect_lt(nrow(src), 8L)
})

test_that("selection by absolute statistic keeps both tails", {
  res <- fake_plot_result(n = 10L)
  src <- attr(gs_plot_bar(res, top_n = 4, sort_by = "stat"), "gs_source")
  expect_true(any(src$stat > 0) && any(src$stat < 0))
})

test_that("facets read database_label", {
  res <- fake_plot_result(n = 3L, databases = c("msigdb_H", "mitopathways"))
  src <- attr(gs_plot_bar(res, top_n = 2, facet = "database"), "gs_source")
  expect_true(all(grepl("^Label for ", as.character(src$.facet_row))))
})

test_that("an empty result gives an empty plot", {
  p <- gs_plot_bar(fake_plot_result(n = 3L)[0, ])
  expect_equal(p$labels$subtitle, "No pathways to plot")
})

test_that("bad input is rejected", {
  expect_error(gs_plot_bar(list()), "`x` must be a gs_result")
  expect_error(gs_plot_bar(fake_plot_result(), padj_max = "low"),
               "`padj_max`")
})

# --- regression: colliding display labels --------------------------------------
# `label` is the y-axis factor, so two pathways formatting to the same string
# collapsed onto ONE row with their bars stacked -- the figure showed n-1 of n
# results, with no warning. Confirmed before the fix: 2 rows of data, 1 distinct
# y level, 2 bars rendered.
test_that("two pathways with the same formatted name get separate axis rows", {
  r <- fake_plot_result(n = 2L)
  r$pathway_id   <- c("HALLMARK_MYC_TARGETS_V1", "KEGG_MYC_TARGETS_V1")
  r$pathway_name <- c("MYC_TARGETS_V1", "MYC_TARGETS_V1")
  p <- gs_plot_bar(r)
  expect_equal(nrow(p$data), 2L)
  expect_equal(length(unique(as.character(p$data$label))), 2L)
  # Disambiguated by the machine id, the one field guaranteed unique.
  expect_true(all(grepl("HALLMARK_MYC_TARGETS_V1|KEGG_MYC_TARGETS_V1",
                        as.character(p$data$label))))
})

# --- regression: the empty plot must answer nrow() -----------------------------
# `ggplot()` with no data leaves `$data` a `waiver`, so `nrow()` returned NULL
# and `if (nrow(p$data) > 0L)` threw "argument is of length zero" in callers.
test_that("an all-filtered-out selection yields a zero-row, not dataless, plot", {
  r <- fake_plot_result(n = 3L)
  p <- gs_plot_bar(r, padj_max = 1e-12)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 0L)
  expect_false(is.null(nrow(p$data)))
})
