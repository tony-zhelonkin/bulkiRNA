test_that("a named numeric vector dispatches to fgsea and returns a gs_result", {
  res <- gs_test(fake_ranks(), fake_gs_db(), contrast = "KO-WT",
                 min_size = 5, max_size = 50)
  expect_s3_class(res, "gs_result")
  expect_gt(nrow(res), 0L)
  expect_true(all(res$stat_type == "NES"))
  expect_true(all(res$method == "fgsea"))
  expect_true(all(res$database == "testdb"))
  expect_true(all(res$contrast == "KO-WT"))
  expect_true("log2err" %in% names(res))
  expect_type(res$log2err, "double")
  expect_length(res$log2err, nrow(res))
  expect_true(all(res$log2err >= 0))
  expect_true("leading_edge" %in% names(res))
  expect_true(is.list(res$leading_edge))
  # names come from the gs_db pathway_names attribute
  expect_true(all(grepl("^Pretty ", res$pathway_name)))
})

test_that("fgsea recovers the planted direction", {
  res <- gs_test(fake_ranks(), fake_gs_db(), min_size = 5, max_size = 50)
  expect_equal(res$direction[res$pathway_id == "SET_UP"], "up")
  expect_equal(res$direction[res$pathway_id == "SET_DOWN"], "down")
})

test_that("a character vector dispatches to ORA", {
  db <- fake_gs_db()
  res <- gs_test(paste0("G", 1:20), db,
                 universe = paste0("G", 1:100),
                 min_size = 5, max_size = 50)
  expect_s3_class(res, "gs_result")
  expect_true(all(res$method == "ora"))
  expect_true(all(res$stat_type == "log2_fold_enrichment"))
  expect_true(all(c("fold_enrichment", "overlap") %in% names(res)))
  expect_false("log2err" %in% names(res))
  up <- res[res$pathway_id == "SET_UP", ]
  expect_equal(up$overlap, 20L)
  expect_gt(up$stat, 0)
})

test_that("ORA rejects query genes outside the universe", {
  expect_error(
    gs_test("ZZZ9", fake_gs_db(), universe = paste0("G", 1:100),
            min_size = 1),
    "universe"
  )
})

test_that("multiple databases become multiple `database` values", {
  db <- list(
    alpha = fake_gs_db(database = "alpha"),
    beta = fake_gs_db(list(BETA_1 = paste0("G", 1:15)), database = "beta")
  )
  res <- gs_test(fake_ranks(), db, min_size = 5, max_size = 50)
  expect_setequal(unique(res$database), c("alpha", "beta"))
})

test_that("pooling contrasts is rbind", {
  db <- fake_gs_db()
  a <- gs_test(fake_ranks(), db, contrast = "A", min_size = 5, max_size = 50)
  b <- gs_test(rev(fake_ranks()), db, contrast = "B",
               min_size = 5, max_size = 50)
  pooled <- rbind(a, b)
  expect_s3_class(pooled, "gs_result")
  expect_equal(nrow(pooled), nrow(a) + nrow(b))
  expect_setequal(unique(pooled$contrast), c("A", "B"))
})

test_that("gs_test errors informatively on an unsupported input", {
  expect_error(gs_test(list(1), fake_gs_db()), "no method")
  expect_error(gs_test(unname(fake_ranks()), fake_gs_db()), "named")
  expect_error(gs_test(fake_ranks(), list(a = 1)), "gs_db")
  expect_error(gs_test(fake_ranks(), fake_gs_db(), method = "ora"), "fgsea")
})

test_that("an empty result is a zero-row gs_result, not NULL", {
  db <- fake_gs_db(list(TOO_BIG = paste0("G", 1:100)))
  res <- gs_test(fake_ranks(), db, min_size = 5, max_size = 10)
  expect_s3_class(res, "gs_result")
  expect_equal(nrow(res), 0L)
})

test_that("gs_test on a gs_matrix runs limma and reports t statistics", {
  skip_if_not_installed("limma")
  m <- matrix(
    c(rnorm(4, 2), rnorm(4, -2), rnorm(8, 0)), nrow = 2, byrow = TRUE,
    dimnames = list(c("SET_UP", "SET_MID"), paste0("s", 1:8))
  )
  gm <- bulkiRNA:::gs_matrix(
    m, database = "testdb", method = "gsva",
    sample_data = fake_sample_data()
  )
  res <- gs_test(gm, design = ~ 0 + group, contrast = "groupKO-groupWT")
  expect_s3_class(res, "gs_result")
  expect_true(all(res$stat_type == "t"))
  expect_true(all(res$method == "gsva"))
  expect_equal(res$contrast[1], "groupKO-groupWT")
  expect_true(all(is.na(res$n_genes)))
})

test_that("gs_test on a gs_matrix needs a design", {
  gm <- bulkiRNA:::gs_matrix(
    matrix(1:8, nrow = 2, dimnames = list(c("a", "b"), paste0("s", 1:4))),
    database = "d", method = "gsva"
  )
  expect_error(gs_test(gm), "design")
})

# --- regressions from the compute-layer review -------------------------------

test_that("gs_test leaves the caller's RNG stream untouched", {
  # `set.seed(seed)` was called on the global stream and never restored, once
  # per database, so anything random *after* a gs_test() call was drawn from the
  # seed-123 stream rather than the script's own.
  db <- fake_gs_db()
  set.seed(999)
  before <- runif(1)
  set.seed(999)
  invisible(gs_test(fake_ranks(), db, min_size = 5, max_size = 50))
  after <- runif(1)
  expect_identical(before, after)
})

test_that("gs_test pins its RNG inside a BiocParallel task", {
  skip_if_not_installed("BiocParallel")
  ranks <- fake_ranks()
  db <- fake_gs_db()

  parent <- gs_test(ranks, db, min_size = 5, max_size = 50, seed = 71L)
  nested <- BiocParallel::bplapply(
    1L,
    function(i) {
      gs_test(ranks, db, min_size = 5, max_size = 50, seed = 71L)
    },
    BPPARAM = BiocParallel::SerialParam()
  )[[1L]]

  expect_gt(nrow(parent), 0L)
  expect_true(all(is.finite(parent$p_value)))
  expect_identical(parent$p_value, nested$p_value)
})

test_that("a partially named db list still labels every database", {
  # `names(db) %||% ...` never fired for a partially named list, because names()
  # returns "" (not NULL) for the unnamed slots, so half the rows carried
  # `database = ""` and grouped under a blank label.
  db <- list(fake_gs_db(database = "alpha"),
             beta = fake_gs_db(database = "betaDB"))
  res <- gs_test(fake_ranks(), db, min_size = 5, max_size = 50)
  expect_setequal(unique(res$database), c("alpha", "beta"))
  expect_false(any(!nzchar(res$database)))
})

test_that("gs_test on a gs_matrix records the scoring method that made it", {
  skip_if_not_installed("limma")
  m <- matrix(
    c(rnorm(4, 2), rnorm(4, -2), rnorm(8, 0)), nrow = 2, byrow = TRUE,
    dimnames = list(c("SET_UP", "SET_MID"), paste0("s", 1:8))
  )
  gm <- bulkiRNA:::gs_matrix(
    m, database = "testdb", method = "ssgsea",
    sample_data = fake_sample_data()
  )
  res <- gs_test(gm, design = ~ 0 + group, contrast = "groupKO-groupWT")
  expect_gt(nrow(res), 0L)
  expect_true(all(res$method == "ssgsea"))
})

test_that("an empty result keeps the method's optional columns", {
  # A contrast that yields no pathways used to lose `leading_edge`, so a
  # per-contrast gs_leading_edge() reported a usage error for an empty answer.
  res <- gs_test(fake_ranks(), fake_gs_db(list(TOO_BIG = paste0("G", 1:100))),
                 min_size = 5, max_size = 10)
  expect_true(all(c("log2err", "leading_edge") %in% names(res)))
  expect_type(res$log2err, "double")
  expect_true(is.list(res$leading_edge))
  expect_length(gs_leading_edge(res), 0L)
})
