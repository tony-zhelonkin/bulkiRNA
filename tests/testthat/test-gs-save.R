test_that("gs_save writes pdf, png and the same-stem source table", {
  p <- gs_plot_dot(fake_plot_result(n = 5L), top_n = 4)
  stem <- file.path(tmp_dir(), "figures", "demo_dot")
  out <- gs_save(p, stem, width = 6, height = 4)

  expect_true(file.exists(paste0(stem, ".pdf")))
  expect_true(file.exists(paste0(stem, ".png")))
  expect_true(file.exists(paste0(stem, ".tsv")))
  expect_equal(out, paste0(stem, c(".pdf", ".png", ".tsv")))
})

test_that("a known extension on the path is stripped, not doubled", {
  p <- gs_plot_bar(fake_plot_result(n = 3L))
  stem <- file.path(tmp_dir(), "demo")
  gs_save(p, paste0(stem, ".pdf"))
  expect_true(file.exists(paste0(stem, ".pdf")))
  expect_false(file.exists(paste0(stem, ".pdf.pdf")))
})

test_that("the written table is the plot's source frame", {
  p <- gs_plot_dot(fake_plot_result(n = 6L), top_n = 3)
  stem <- file.path(tmp_dir(), "demo")
  gs_save(p, stem, formats = "png")
  tab <- utils::read.delim(paste0(stem, ".tsv"), check.names = FALSE)
  expect_equal(nrow(tab), 3L)
  expect_true("pathway_id" %in% names(tab))
})

test_that("list columns are collapsed with / on the way out", {
  df <- data.frame(a = 1:2)
  df$leading_edge <- list(c("G1", "G2"), "G3")
  f <- file.path(tmp_dir(), "t.tsv")
  bulkiRNA:::.gs_save_tsv(df, f)
  expect_equal(utils::read.delim(f)$leading_edge, c("G1/G2", "G3"))
})

test_that("table = FALSE writes images only", {
  p <- gs_plot_bar(fake_plot_result(n = 3L))
  stem <- file.path(tmp_dir(), "demo")
  out <- gs_save(p, stem, formats = "png", table = FALSE)
  expect_equal(out, paste0(stem, ".png"))
  expect_false(file.exists(paste0(stem, ".tsv")))
})

test_that("an explicit data argument overrides the plot's source", {
  p <- gs_plot_bar(fake_plot_result(n = 5L), top_n = 5)
  stem <- file.path(tmp_dir(), "demo")
  gs_save(p, stem, formats = "png", data = data.frame(x = 1:2))
  expect_equal(names(utils::read.delim(paste0(stem, ".tsv"))), "x")
})

test_that("bad arguments error by name", {
  expect_error(gs_save("not a plot", "x"), "`plot`")
  expect_error(gs_save(gs_plot_bar(fake_plot_result()), ""), "`path`")
  expect_error(
    gs_save(gs_plot_bar(fake_plot_result()), "x", formats = "gif"),
    "Unsupported `formats`"
  )
})

test_that(".gs_write_log reports counts and both directions", {
  res <- fake_plot_result(n = 8L)
  f <- file.path(tmp_dir(), "log.txt")
  bulkiRNA:::.gs_write_log(res, f, padj_cutoff = 0.5)
  txt <- readLines(f)
  expect_true(any(grepl("^Total pathways: 8", txt)))
  expect_true(any(grepl("UPREGULATED PATHWAYS", txt)))
  expect_true(any(grepl("DOWNREGULATED PATHWAYS", txt)))
  expect_true(any(grepl("End of log", txt)))
})

test_that(".gs_write_log says so when nothing is significant", {
  res <- fake_plot_result(n = 4L)
  f <- file.path(tmp_dir(), "log.txt")
  bulkiRNA:::.gs_write_log(res, f, padj_cutoff = 1e-12)
  expect_true(any(grepl("No significant pathways", readLines(f))))
})

test_that(".gs_plot_all writes the standard set per database", {
  res <- fake_plot_result(n = 6L,
                          databases = c("msigdb_H", "mitopathways"))
  d <- tmp_dir()
  out <- bulkiRNA:::.gs_plot_all(res, d, top_n = 3)
  expect_true(all(file.exists(out)))
  expect_true(dir.exists(file.path(d, "msigdb_H")))
  expect_true(dir.exists(file.path(d, "mitopathways")))
  # 4 figures x (pdf + png + tsv) + 1 log, per database
  expect_equal(length(out), 2L * (4L * 3L + 1L))
})
