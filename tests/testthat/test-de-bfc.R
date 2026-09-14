test_that("de_bfc_plot returns a ggplot and needs logFC + B", {
  de <- fake_de_table()
  expect_s3_class(de_bfc_plot(de), "ggplot")
  expect_error(de_bfc_plot(de[, "logFC", drop = FALSE]), "`B`")
})

test_that("the B threshold line is always inside the y limits", {
  de <- fake_de_table()
  de$B <- seq(-20, -10, length.out = nrow(de))   # nothing crosses B = 0
  p <- de_bfc_plot(de, B_cutoff = 0, y_padding = 1)

  ylim <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
  expect_gt(ylim[2], 0)

  labels <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) d$label))
  expect_true(any(grepl("threshold", labels)))
})

test_that("categories combine B and fold change", {
  de <- data.frame(logFC = c(3, 0.1, 3, 0.1), B = c(2, 2, -2, -2),
                   row.names = paste0("g", 1:4))
  p <- de_bfc_plot(de, fc_cutoff = 1, B_cutoff = 0, label_method = "none")
  cats <- ggplot2::ggplot_build(p)$plot$data$cat
  expect_equal(cats, c("B-statistic & Log2FC", "B-statistic", "Log2FC", "NS"))
})

test_that("top labels are ranked by decreasing B on each side", {
  de <- fake_de_table(n = 20)
  p  <- de_bfc_plot(de, top_n = 2)
  repel <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(repel, 1L)
  expect_equal(nrow(repel[[1]]$data), 4L)
})

# --- regression: highlight overlapping the top-N ------------------------------
# `rbind()` uniquifies colliding row names, so a highlighted gene that was also
# top-N was appended as "Gene11"/"Gene51" -- names of genes that do not exist --
# and the following `!duplicated(rownames(...))` could never fire. The invented
# string was drawn on the figure and, being unequal to the real id, was not
# bolded either.
test_that("de_bfc_plot never invents a gene name when highlight overlaps top-N", {
  de <- data.frame(logFC = c(-3, -1, 0.2, 1.5, 4, 2.5, -2.2),
                   B = c(4, -1, -3, 0.5, 6, 3, 2),
                   AveExpr = 1:7, P.Value = 0.01, adj.P.Val = 0.02,
                   row.names = paste0("Gene", 1:7))
  p <- de_bfc_plot(de, fc_cutoff = 1, top_n = 3,
                   highlight_gene = c("Gene5", "Gene1"))
  drawn <- unlist(lapply(p$layers, function(l) {
    d <- l$data
    if (is.data.frame(d) && ".gene_label" %in% names(d)) d$.gene_label else NULL
  }))
  expect_true(all(drawn %in% rownames(de)))
  expect_false(any(c("Gene11", "Gene51") %in% drawn))
  # Both highlights must still be labelled exactly once.
  expect_equal(sum(drawn == "Gene1"), 1L)
  expect_equal(sum(drawn == "Gene5"), 1L)
})
