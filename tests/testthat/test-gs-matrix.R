test_that("gs_matrix carries its metadata", {
  m <- fake_gs_matrix()
  expect_s3_class(m, "gs_matrix")
  expect_identical(bulkiRNA:::gs_database(m), "testdb")
  expect_identical(bulkiRNA:::gs_score_type(m), "gsva")
  expect_identical(names(bulkiRNA:::gs_pathway_names(m)), rownames(m))
  expect_identical(rownames(bulkiRNA:::gs_sample_data(m)), colnames(m))
})

test_that("dimnames are required and must be unique", {
  m <- matrix(1:4, nrow = 2)
  expect_error(bulkiRNA:::gs_matrix(m, "d", "gsva"), "row names")
  dimnames(m) <- list(c("P1", "P1"), c("s1", "s2"))
  expect_error(bulkiRNA:::gs_matrix(m, "d", "gsva"), "duplicated pathway")
})

test_that("sample_data is aligned to column order by row name", {
  m <- fake_gs_matrix(2L, 3L)
  shuffled <- data.frame(group = c("c", "a", "b"),
                         row.names = c("s3", "s1", "s2"))
  out <- bulkiRNA:::gs_matrix(unclass(m)[, , drop = FALSE], "d", "gsva",
                              sample_data = shuffled)
  expect_identical(bulkiRNA:::gs_sample_data(out)$group, c("a", "b", "c"))
})

test_that("mismatched sample_data is an error", {
  m <- fake_gs_matrix(2L, 3L)
  expect_error(
    bulkiRNA:::gs_matrix(unclass(m), "d", "gsva",
                         sample_data = data.frame(group = c("a", "b"))),
    "sample_data"
  )
})

test_that("subsetting keeps the class and subsets the metadata", {
  m <- fake_gs_matrix(4L, 4L)
  sub <- m[1:2, c("s1", "s3")]
  expect_s3_class(sub, "gs_matrix")
  expect_identical(dim(sub), c(2L, 2L))
  expect_identical(names(bulkiRNA:::gs_pathway_names(sub)), c("SET_1", "SET_2"))
  expect_identical(rownames(bulkiRNA:::gs_sample_data(sub)), c("s1", "s3"))
})

test_that("subsetting to one row drops to a plain vector", {
  m <- fake_gs_matrix(3L, 3L)
  v <- m[1, ]
  expect_false(inherits(v, "gs_matrix"))
  expect_type(v, "double")
})

test_that("as_tibble gives long form joined to sample metadata", {
  m <- fake_gs_matrix(2L, 4L)
  tb <- tibble::as_tibble(m)
  expect_identical(nrow(tb), 8L)
  expect_true(all(
    c("pathway_id", "pathway_name", "sample", "score", "score_type", "group")
      %in% names(tb)
  ))
  expect_identical(tb$score[1], unclass(m)[1, 1])
})

test_that("summary reports dimensions and score range", {
  s <- summary(fake_gs_matrix(3L, 4L))
  expect_identical(s$n_pathways, 3L)
  expect_identical(s$n_samples, 4L)
  expect_identical(s$n_na, 0L)
})

test_that("print shows the gs_matrix header", {
  expect_output(print(fake_gs_matrix()), "gs_matrix")
})

# --- regressions from the compute-layer review -------------------------------

test_that("subsetting to zero pathways or samples yields an empty gs_matrix", {
  # `[` funnelled through the full gs_matrix() constructor, whose "must have row
  # names" check rejects character(0) dimnames, so a pathway filter that matched
  # nothing aborted the script instead of returning an empty matrix.
  m <- fake_gs_matrix(3L, 4L)
  none <- m[rownames(m) %in% "not_a_pathway", , drop = FALSE]
  expect_s3_class(none, "gs_matrix")
  expect_identical(dim(none), c(0L, 4L))
  expect_length(bulkiRNA:::gs_pathway_names(none), 0L)

  no_samp <- m[, character(0), drop = FALSE]
  expect_s3_class(no_samp, "gs_matrix")
  expect_identical(dim(no_samp), c(3L, 0L))
  expect_identical(nrow(bulkiRNA:::gs_sample_data(no_samp)), 0L)
})

test_that("a repeated subscript is allowed and an unmatched name is named", {
  m <- fake_gs_matrix(3L, 4L)
  dup <- m[c(1L, 1L), , drop = FALSE]
  expect_s3_class(dup, "gs_matrix")
  expect_identical(dim(dup), c(2L, 4L))

  # Was a raw "subscript out of bounds".
  expect_error(m[c("SET_1", "nope"), ], "nope")
  expect_error(m[, c("s1", "nosuch")], "nosuch")
})
