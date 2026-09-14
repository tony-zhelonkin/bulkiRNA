res3 <- function() {
  bulkiRNA:::gs_result(
    data.frame(
      pathway_id = c("A", "B", "C"),
      pathway_name = c("Alpha set", "Beta set", "Gamma set"),
      n_genes = c(20L, 30L, 40L),
      n_genes_tested = c(18L, 25L, 35L),
      stat = c(2.1, -1.8, 0.4),
      p_value = c(0.001, 0.01, 0.4),
      padj = c(0.01, 0.04, 0.6),
      leading_edge = I(list(c("g1", "g2"), "g3", character(0))),
      stringsAsFactors = FALSE
    ),
    database = "demo", contrast = "KO-WT", method = "fgsea",
    stat_type = "NES"
  )
}

test_that("gs_filter applies each constraint and keeps the class", {
  r <- res3()
  expect_equal(nrow(gs_filter(r, padj = 0.05)), 2L)
  expect_equal(gs_filter(r, direction = "up")$pathway_id, c("A", "C"))
  expect_equal(nrow(gs_filter(r, stat = 2)), 1L)
  expect_equal(nrow(gs_filter(r, min_genes = 20, max_genes = 30)), 1L)
  expect_equal(gs_filter(r, pattern = "beta")$pathway_id, "B")
  expect_s3_class(gs_filter(r, padj = 0.05), "gs_result")
  expect_equal(nrow(gs_filter(r)), 3L)
  expect_equal(nrow(gs_filter(r, padj = 1e-9)), 0L)
})

test_that("gs_filter rejects a bad direction", {
  expect_error(gs_filter(res3(), direction = "UP"), "\"up\"")
})

test_that("gs_top selects, ordered by the ranking column", {
  r <- res3()
  expect_equal(gs_top(r, n = 2)$pathway_id, c("A", "B"))
  expect_equal(gs_top(r, n = 1, by = "stat")$pathway_id, "A")
  both <- gs_top(r, n = 1, by_direction = TRUE)
  expect_setequal(both$pathway_id, c("A", "B"))
})

test_that("gs_top applies n per group", {
  r <- rbind(res3(), res3()[1:2, ])
  r$contrast[4:5] <- "X-Y"
  r <- bulkiRNA:::gs_result(as.data.frame(r))
  expect_equal(nrow(gs_top(r, n = 1)), 2L)
  expect_equal(nrow(gs_top(r, n = 1, per = character(0))), 1L)
})

test_that("gs_top keeps an allowed missing direction as its own group", {
  r <- res3()
  r$direction <- NA_character_
  r <- bulkiRNA:::gs_result(as.data.frame(r))

  top <- gs_top(r, n = 2L, by_direction = TRUE, per = character(0))
  expect_equal(nrow(top), 2L)
  expect_true(all(is.na(top$direction)))
})

test_that("gs_split is the inverse of the rbind that pools", {
  a <- res3()
  b <- res3()
  b$contrast <- "X-Y"
  pooled <- rbind(a, b)
  parts <- gs_split(pooled)
  expect_length(parts, 2L)
  expect_true(all(vapply(parts, inherits, logical(1), "gs_result")))
  expect_equal(nrow(do.call(rbind, unname(parts))), nrow(pooled))
})

test_that("gs_leading_edge extracts, filters and pools genes", {
  r <- res3()
  le <- gs_leading_edge(r)
  expect_named(le, c("A", "B", "C"))
  expect_equal(le$A, c("g1", "g2"))
  expect_named(gs_leading_edge(r, padj = 0.05), c("A", "B"))
  expect_named(gs_leading_edge(r, top_n = 1), "A")
  expect_setequal(gs_leading_edge(r, unique_genes = TRUE), c("g1", "g2", "g3"))
})

test_that("gs_leading_edge errors when there is no such column", {
  r <- res3()
  r$leading_edge <- NULL
  expect_error(gs_leading_edge(bulkiRNA:::gs_result(as.data.frame(r))),
               "no `leading_edge`")
})

test_that("the result ops reject non-gs_result input", {
  expect_error(gs_filter(data.frame(a = 1)), "gs_result")
  expect_error(gs_top(data.frame(a = 1)), "gs_result")
  expect_error(gs_split(data.frame(a = 1)), "gs_result")
})
