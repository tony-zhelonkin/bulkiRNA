test_that("gs_score returns a gs_matrix carrying its metadata", {
  skip_if_not_installed("GSVA")
  sc <- gs_score(fake_expr(), fake_gs_db(), min_size = 5, max_size = 50,
                 sample_data = fake_sample_data())
  expect_s3_class(sc, "gs_matrix")
  expect_equal(ncol(sc), 8L)
  expect_true(all(rownames(sc) %in% c("SET_UP", "SET_DOWN", "SET_MID")))
  expect_equal(attr(sc, "database"), "testdb")
  expect_equal(attr(sc, "score_type"), "gsva")
  expect_true(all(grepl("^Pretty ", attr(sc, "pathway_names"))))
  expect_equal(nrow(attr(sc, "sample_data")), 8L)
})

test_that("gs_score -> gs_test is the GSVA to limma pipeline", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("limma")
  sc <- gs_score(fake_expr(), fake_gs_db(), min_size = 5, max_size = 50,
                 sample_data = fake_sample_data())
  res <- gs_test(sc, design = ~ 0 + group, contrast = "groupKO-groupWT")
  expect_s3_class(res, "gs_result")
  expect_true(all(res$stat_type == "t"))
  # the planted effect is in SET_UP, high in the KO samples
  expect_gt(res$stat[res$pathway_id == "SET_UP"], 0)
  expect_lt(res$padj[res$pathway_id == "SET_UP"], 0.05)
})

test_that("gs_score validates its input", {
  expect_error(gs_score(1:10, fake_gs_db()), "numeric matrix")
  m <- matrix(1, 2, 2)
  expect_error(gs_score(m, fake_gs_db()), "row names")
  expect_error(
    gs_score(fake_expr(), list(a = fake_gs_db(), b = fake_gs_db())),
    "one database at a time"
  )
})
