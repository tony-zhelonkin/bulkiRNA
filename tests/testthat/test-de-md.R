test_that("de_md_plot returns a ggplot with a loess trend", {
  de <- fake_de_table()
  fit <- list(Amean = de$AveExpr)
  p <- de_md_plot(fit, coef = "KO_vs_WT", de_results = de)

  expect_s3_class(p, "ggplot")
  # The geom_smooth trend is load-bearing: it is what exposes intensity bias.
  expect_true(any(vapply(p$layers, function(l) inherits(l$stat, "StatSmooth"),
                         logical(1))))
})

test_that("the title defaults to the coefficient name", {
  de <- fake_de_table()
  p <- de_md_plot(list(Amean = de$AveExpr), coef = "KO_vs_WT", de_results = de)
  expect_equal(p$labels$title, "MD plot: KO_vs_WT")
})

test_that("significance is called on FDR alone", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_md_plot(list(Amean = de$AveExpr), coef = 1, de_results = de,
                   fdr_cutoff = 0.05)
  built <- ggplot2::ggplot_build(p)
  # layer 2 is the coloured significant cloud
  expect_equal(nrow(built$data[[2]]), sum(de$adj.P.Val <= 0.05))
})

test_that("AveExpr is filled from fit$Amean when absent", {
  de <- fake_de_table()
  amean <- de$AveExpr
  de$AveExpr <- NULL
  p <- de_md_plot(list(Amean = amean), coef = 1, de_results = de)
  expect_s3_class(p, "ggplot")
})

test_that("a missing required column is an error naming the column", {
  de <- fake_de_table()
  de$adj.P.Val <- NULL
  expect_error(de_md_plot(list(Amean = de$AveExpr), coef = 1, de_results = de),
               "adj.P.Val")
})

test_that("quadrant counts sum to the number of significant genes", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_md_plot(list(Amean = de$AveExpr), coef = 1, de_results = de,
                   show_quadrant_counts = TRUE)
  txt <- Filter(function(l) inherits(l$geom, "GeomText") &&
                  !inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(txt, 1L)
  expect_equal(sum(txt[[1]]$data$lbl), sum(de$adj.P.Val <= 0.05))
})
