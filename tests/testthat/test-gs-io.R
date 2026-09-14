io_res <- function() {
  bulkiRNA:::gs_result(
    data.frame(
      pathway_id = c("A", "B", "C"),
      pathway_name = c("Alpha set", "Beta set", "Gamma set"),
      contrast = c("KO-WT", "KO-WT", "X vs Y"),
      n_genes = c(20L, 30L, 40L),
      n_genes_tested = c(18L, 25L, 35L),
      stat = c(2.1, -1.8, 0.4),
      p_value = c(0.001, 0.01, 0.4),
      padj = c(0.01, 0.04, 0.6),
      leading_edge = I(list(c("g1", "g2"), "g3", character(0))),
      stringsAsFactors = FALSE
    ),
    database = "demo", method = "fgsea", stat_type = "NES"
  )
}

test_that("gs_write lays out by_contrast/ and _overview/", {
  d <- file.path(tempdir(), "gsw1")
  unlink(d, recursive = TRUE)
  out <- gs_write(io_res(), d)
  expect_true(file.exists(file.path(d, "_overview", "gsea_all.tsv")))
  expect_true(file.exists(file.path(d, "_overview", "gsea_summary.tsv")))
  expect_true(file.exists(file.path(d, "by_contrast", "KO-WT", "gsea_demo.tsv")))
  # the contrast label is slugged into a safe directory name
  expect_true(file.exists(file.path(d, "by_contrast", "X_vs_Y", "gsea_demo.tsv")))
  expect_length(attr(out, "files"), 4L)
})

test_that("the write/read round trip is lossless", {
  d <- file.path(tempdir(), "gsw2")
  unlink(d, recursive = TRUE)
  r <- io_res()
  gs_write(r, d)
  back <- gs_read(d)
  expect_s3_class(back, "gs_result")
  expect_equal(nrow(back), nrow(r))
  expect_equal(back$pathway_id, r$pathway_id)
  expect_equal(back$stat, r$stat)
  expect_equal(back$padj, r$padj)
  expect_equal(back$direction, r$direction)
  expect_equal(back$leading_edge, unclass(r$leading_edge))
})

test_that("gs_read falls back to the by_contrast files", {
  d <- file.path(tempdir(), "gsw3")
  unlink(d, recursive = TRUE)
  gs_write(io_res(), d, overview = FALSE)
  back <- gs_read(d)
  expect_equal(nrow(back), 3L)
})

test_that("gs_read errors when there is nothing to read", {
  expect_error(gs_read(file.path(tempdir(), "nope-not-here")), "No gs_result")
})

test_that("gs_write validates its arguments", {
  expect_error(gs_write(data.frame(a = 1), tempdir()), "gs_result")
  expect_error(gs_write(io_res(), ""), "non-empty")
})

# --- regressions from the IO review ------------------------------------------

test_that("gs_read warns about by_contrast files the last write did not make", {
  # gs_write() only ever added files, and gs_read()'s fallback globs
  # by_contrast/*/<name>_*.tsv, so a contrast dropped from the analysis was
  # silently row-bound back in from the previous run.
  d <- file.path(tempdir(), "gsw-stale")
  unlink(d, recursive = TRUE)
  gs_write(io_res(), d, overview = FALSE)

  dropped <- dplyr::filter(io_res(), .data$contrast == "KO-WT")
  gs_write(dropped, d, overview = FALSE)

  expect_warning(back <- gs_read(d), "not written by the last")
  expect_setequal(unique(back$contrast), c("KO-WT", "X vs Y"))
})

test_that("gs_write(prune = TRUE) clears the stale contrast", {
  d <- file.path(tempdir(), "gsw-prune")
  unlink(d, recursive = TRUE)
  gs_write(io_res(), d, overview = FALSE)

  dropped <- dplyr::filter(io_res(), .data$contrast == "KO-WT")
  gs_write(dropped, d, overview = FALSE, prune = TRUE)

  expect_false(dir.exists(file.path(d, "by_contrast", "X_vs_Y")))
  expect_silent(back <- gs_read(d))
  expect_identical(unique(back$contrast), "KO-WT")
})

test_that("a matching write/read pair does not warn", {
  d <- file.path(tempdir(), "gsw-clean")
  unlink(d, recursive = TRUE)
  gs_write(io_res(), d, overview = FALSE)
  expect_silent(gs_read(d))
})
