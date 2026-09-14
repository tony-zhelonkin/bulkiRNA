test_that("gs_coregulation returns the GESECA gs_result contract", {
  expect_identical(formals(gs_coregulation)$center, TRUE)
  expect_identical(formals(gs_coregulation)$scale, FALSE)
  expect_identical(formals(gs_coregulation)$sample_size, 101L)
  expect_identical(formals(gs_coregulation)$seed, 123L)

  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 91L
  )

  expect_s3_class(res, "gs_result")
  expect_gt(nrow(res), 0L)
  expect_true(all(res$method == "geseca"))
  expect_true(all(res$stat_type == "pct_var"))
  expect_true(all(res$contrast == "coregulation"))
  expect_true(all(res$database == "coreg_demo"))
  expect_true(all(is.na(res$direction)))
  expect_true(all(res$stat >= 0))
  expect_type(res$log2err, "double")
  expect_length(res$log2err, nrow(res))
  expect_true(all(grepl("^Pretty ", res$pathway_name)))
  expect_silent(bulkiRNA:::validate_gs_result(res))
})

test_that("the planted coregulated set outranks a scrambled set", {
  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 92L
  )

  score <- stats::setNames(res$stat, res$pathway_id)
  expect_gt(score[["COREGULATED"]], score[["SCRAMBLED"]])
  expect_identical(
    gs_top(res, n = 1L, by = "stat", per = character(0))$pathway_id,
    "COREGULATED"
  )
})

test_that("default centering and no scaling reproduce the GESECA score", {
  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 93L
  )

  centered <- fixture$expr - rowMeans(fixture$expr)
  genes <- fixture$coregulated
  expected <- sum(colSums(centered[genes, , drop = FALSE])^2) /
    length(genes) / sum(centered^2) * 100
  observed <- res$stat[res$pathway_id == "COREGULATED"]
  expect_equal(observed, expected, tolerance = 1e-12)

  shifted <- fixture$expr + seq_len(nrow(fixture$expr))
  shifted_res <- gs_coregulation(
    shifted, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 93L
  )
  expect_equal(res$stat, shifted_res$stat, tolerance = 1e-12)
})

test_that("unsigned direction stays unavailable throughout result operations", {
  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 94L
  )

  expect_true(all(is.na(res$direction)))
  expect_equal(nrow(gs_filter(res, direction = "up")), 0L)
  expect_equal(nrow(gs_filter(res, direction = "down")), 0L)
  expect_equal(nrow(gs_filter(res, direction = "ns")), 0L)
  expect_equal(nrow(gs_top(res, n = 2L, by_direction = TRUE,
                           per = character(0))), 2L)

  overview <- summary(res, padj_cutoff = Inf)
  expect_true(all(is.na(overview$n_up)))
  expect_true(all(is.na(overview$n_down)))
})

test_that("existing renderers and persistence accept pct_var unchanged", {
  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 95L
  )

  dot <- gs_plot_dot(res, top_n = 4L)
  bar <- gs_plot_bar(res, top_n = 4L)
  expect_s3_class(dot, "ggplot")
  expect_s3_class(bar, "ggplot")
  expect_identical(dot$labels$x, "% variance explained")
  expect_identical(bar$labels$y, "% variance explained")
  expect_equal(nrow(attr(dot, "gs_source")), 4L)
  expect_equal(nrow(attr(bar, "gs_source")), 4L)

  path <- tempfile("coregulation-write-")
  expect_silent(out <- gs_write(res, path))
  expect_true(all(file.exists(attr(out, "files"))))
})

test_that("pct_var requires the deliberate master-table override", {
  fixture <- fake_coregulation_input()
  res <- gs_coregulation(
    fixture$expr, fixture$db,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 96L
  )

  expect_error(gs_to_master(res), "got.*pct_var")
  expect_error(gs_to_master(res), "stat_as_nes = TRUE", fixed = TRUE)

  master <- gs_to_master(res, stat_as_nes = TRUE)
  expect_identical(
    names(master),
    c(
      "pathway_id", "pathway_name", "database", "contrast", "nes",
      "pvalue", "padj", "set_size", "leading_edge_size", "gene_ratio",
      "core_enrichment", "genes_full_set", "direction", "neg_log_padj"
    )
  )
  expect_identical(ncol(master), 14L)
  expect_true(all(is.na(master$direction)))
  expect_silent(gs_validate_master(master))
})

test_that("multiple databases remain distinct and empty answers stay typed", {
  fixture <- fake_coregulation_input()
  dbs <- list(
    first = fixture$db,
    second = gsdb_register(
      list(ONLY = fixture$coregulated), database = "ignored",
      species = "Homo sapiens"
    )
  )
  res <- gs_coregulation(
    fixture$expr, dbs,
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 97L
  )
  expect_setequal(unique(res$database), c("first", "second"))

  empty <- gs_coregulation(
    fixture$expr,
    gsdb_register(
      list(TOO_SMALL = fixture$coregulated[1:2]), database = "empty",
      species = "Homo sapiens"
    ),
    min_size = 10L, max_size = 50L,
    sample_size = 21L, seed = 97L
  )
  expect_s3_class(empty, "gs_result")
  expect_equal(nrow(empty), 0L)
  expect_type(empty$direction, "character")
  expect_type(empty$log2err, "double")
})

test_that("gs_coregulation rejects ambiguous matrices and bad controls", {
  fixture <- fake_coregulation_input()
  expr <- fixture$expr

  expect_error(gs_coregulation(as.data.frame(expr), fixture$db),
               "numeric matrix")
  expect_error(gs_coregulation(unname(expr), fixture$db), "row names")
  duplicated <- expr
  rownames(duplicated)[2L] <- rownames(duplicated)[1L]
  expect_error(gs_coregulation(duplicated, fixture$db), "duplicated")
  nonfinite <- expr
  nonfinite[1L] <- NA_real_
  expect_error(gs_coregulation(nonfinite, fixture$db), "finite")
  expect_error(gs_coregulation(expr, fixture$db, center = NA), "`center`")
  expect_error(gs_coregulation(expr, fixture$db, scale = NA), "`scale`")
  expect_error(gs_coregulation(expr, fixture$db, min_size = 0), "`min_size`")
  expect_error(gs_coregulation(expr, fixture$db,
                               min_size = 50L, max_size = 10L),
               "must not exceed")
  expect_error(gs_coregulation(expr, fixture$db, sample_size = 1.5),
               "`sample_size`")
  expect_error(gs_coregulation(expr, fixture$db, eps = -1), "`eps`")
  expect_error(gs_coregulation(expr, fixture$db, seed = NA), "`seed`")
  expect_error(gs_coregulation(expr, fixture$db, verbose = NA), "`verbose`")

  expect_error(gs_coregulation(fake_gs_matrix(), fixture$db),
               "pathways in its rows")
})

test_that("gs_coregulation pins RNG inside a BiocParallel task", {
  skip_if_not_installed("BiocParallel")
  fixture <- fake_coregulation_input()
  call <- function() {
    gs_coregulation(
      fixture$expr, fixture$db,
      min_size = 10L, max_size = 50L,
      sample_size = 21L, seed = 98L
    )
  }

  parent <- call()
  nested <- BiocParallel::bplapply(
    1L, function(i) call(), BPPARAM = BiocParallel::SerialParam()
  )[[1L]]

  expect_gt(nrow(parent), 0L)
  expect_identical(parent$p_value, nested$p_value)
  expect_identical(parent$log2err, nested$log2err)
})
