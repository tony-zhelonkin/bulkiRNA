test_that("gs_result fills scalar metadata and derives direction", {
  res <- fake_gs_result(3L)
  expect_s3_class(res, "gs_result")
  expect_true(all(c("database", "contrast", "method", "stat_type") %in% names(res)))
  expect_identical(unique(res$database), "testdb")
  expect_identical(res$direction, c("down", "ns", "up"))
})

test_that("core columns are present, typed, and canonically ordered", {
  res <- fake_gs_result(2L)
  core <- bulkiRNA:::.gs_core_cols
  expect_identical(names(res)[seq_along(core)], core)
  expect_type(res$pathway_id, "character")
  expect_type(res$n_genes, "integer")
  expect_type(res$stat, "double")
})

test_that("a missing core column is an error naming the column", {
  df <- fake_gs_df()
  df$padj <- NULL
  expect_error(
    bulkiRNA:::gs_result(df, database = "d", contrast = "c",
                         method = "fgsea", stat_type = "NES"),
    "padj"
  )
})

test_that("stat_type is a closed vocabulary", {
  expect_error(
    bulkiRNA:::gs_result(fake_gs_df(), database = "d", contrast = "c",
                         method = "fgsea", stat_type = "odds_ratio"),
    "Unknown stat_type"
  )
  expect_identical(
    names(gs_stat_types()),
    c("NES", "t", "log2_fold_enrichment", "signed_log10p", "pct_var")
  )
  expect_identical(gs_stat_types()[["pct_var"]], "% variance explained")
})

test_that("validate_gs_result rejects a bad direction", {
  res <- fake_gs_result()
  res$direction <- "UP"
  expect_error(bulkiRNA:::validate_gs_result(res), "direction")
})

test_that("optional columns survive and stay after the core block", {
  df <- fake_gs_df(2L)
  df$es <- c(0.4, -0.4)
  df$log2err <- c(0.3, 0.8)
  res <- bulkiRNA:::gs_result(df, database = "d", contrast = "c",
                              method = "fgsea", stat_type = "NES")
  res$leading_edge <- list(c("A", "B"), "C")
  res <- bulkiRNA:::validate_gs_result(res)
  expect_true(all(c("es", "log2err") %in% names(res)))
  expect_identical(res$log2err, c(0.3, 0.8))
  expect_gt(match("es", names(res)), length(bulkiRNA:::.gs_core_cols))
})

test_that("an infinite log2err is valid and is not missing", {
  df <- fake_gs_df(2L)
  df$log2err <- c(0.3, Inf)

  expect_silent(
    res <- bulkiRNA:::gs_result(
      df, database = "d", contrast = "c",
      method = "fgsea", stat_type = "NES"
    )
  )
  expect_silent(bulkiRNA:::validate_gs_result(res))
  expect_true(is.infinite(res$log2err[2L]))
  expect_false(is.na(res$log2err[2L]))
  expect_s3_class(dplyr::arrange(res, .data$log2err), "gs_result")
})

test_that("gs_to_master deliberately drops log2err", {
  df <- fake_gs_df(2L)
  df$log2err <- c(0.3, Inf)
  res <- bulkiRNA:::gs_result(
    df, database = "d", contrast = "c",
    method = "fgsea", stat_type = "NES"
  )

  master <- gs_to_master(res)
  expect_false("log2err" %in% names(master))
  expect_identical(ncol(master), 14L)
})

test_that("rbind pools contrasts and keeps the class", {
  pooled <- rbind(fake_gs_result(2L, "A-B"), fake_gs_result(2L, "C-D"))
  expect_s3_class(pooled, "gs_result")
  expect_identical(nrow(pooled), 4L)
  expect_setequal(unique(pooled$contrast), c("A-B", "C-D"))
})

test_that("rbind fills optional columns absent from one input", {
  a <- fake_gs_result(1L, "A-B")
  b <- fake_gs_result(1L, "C-D")
  b$es <- 0.5
  pooled <- rbind(a, b)
  expect_true("es" %in% names(pooled))
  expect_true(is.na(pooled$es[1]))
})

test_that("dplyr verbs keep the class, and drop it when a core column goes", {
  res <- fake_gs_result(3L)
  expect_s3_class(dplyr::filter(res, stat > 0), "gs_result")
  expect_s3_class(dplyr::arrange(res, padj), "gs_result")
  thin <- dplyr::select(res, "pathway_id", "stat")
  expect_false(inherits(thin, "gs_result"))
  expect_s3_class(thin, "tbl_df")
})

test_that("as_tibble drops the subclass without losing data", {
  res <- fake_gs_result(3L)
  tb <- tibble::as_tibble(res)
  expect_false(inherits(tb, "gs_result"))
  expect_identical(tb$stat, res$stat)
})

test_that("summary counts pathways and significant hits per database x contrast", {
  s <- summary(rbind(fake_gs_result(3L, "A-B"), fake_gs_result(3L, "C-D")),
               padj_cutoff = 0.05)
  expect_identical(nrow(s), 2L)
  expect_true(all(s$n_pathways == 3L))
  expect_true(all(s$n_sig >= s$n_up))
})

test_that("gs_stat_label reads the axis label off stat_type", {
  expect_identical(bulkiRNA:::gs_stat_label(fake_gs_result()), "NES")
  mixed <- fake_gs_result(2L)
  mixed$stat_type <- c("NES", "t")
  expect_identical(bulkiRNA:::gs_stat_label(mixed), "statistic")
})

test_that("print shows the gs_result header", {
  expect_output(print(fake_gs_result()), "gs_result")
})

# --- regressions from the compute-layer review -------------------------------

test_that("a dplyr verb that breaks the contract downgrades to a tibble", {
  # `dplyr_reconstruct` re-attached the class after checking column *presence*
  # only, so validate_gs_result() was unreachable from every verb: a character
  # `padj` or a restored "Up"/"Down" capitalisation still claimed to be a
  # gs_result, and gs_filter(direction = "up") then returned zero rows silently.
  res <- fake_gs_result(3L)
  expect_warning(bad <- dplyr::mutate(res, padj = "oops"), "gs_result contract")
  expect_false(inherits(bad, "gs_result"))
  expect_s3_class(bad, "tbl_df")

  expect_warning(dir_bad <- dplyr::mutate(res, direction = "Up"), "gs_result contract")
  expect_false(inherits(dir_bad, "gs_result"))

  # A verb that keeps the contract is untouched, and silent.
  expect_silent(ok <- dplyr::arrange(res, .data$padj))
  expect_s3_class(ok, "gs_result")
})
