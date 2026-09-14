test_that("gs_ranks builds a named, decreasing vector from a DE frame", {
  de <- data.frame(t = c(1, 3, 2), row.names = c("A", "B", "C"))
  r <- gs_ranks(de)
  expect_named(r, c("B", "C", "A"))
  expect_equal(unname(r), c(3, 2, 1))
})

test_that("gs_ranks takes genes from a column or a vector", {
  de <- data.frame(t = c(1, 2), sym = c("A", "B"))
  expect_named(gs_ranks(de, genes = "sym"), c("B", "A"))
  expect_named(gs_ranks(de, genes = c("X", "Y")), c("Y", "X"))
})

test_that("gs_ranks drops NA and infinite values with a warning", {
  de <- data.frame(t = c(1, NA, Inf, 2), row.names = c("A", "B", "C", "D"))
  expect_warning(expect_warning(r <- gs_ranks(de)))
  expect_named(r, c("D", "A"))
})

test_that("gs_ranks passes a named numeric vector through", {
  v <- c(B = 2, A = 1)
  expect_equal(gs_ranks(v), c(B = 2, A = 1))
  expect_error(gs_ranks(c(1, 2)), "named")
})

test_that("gs_ranks collapses duplicated identifiers on request", {
  de <- data.frame(t = c(1, -5, 2), row.names = NULL)
  rownames(de) <- NULL
  r <- gs_ranks(de, genes = c("A", "A", "B"), collapse = "max_abs")
  expect_equal(r[["A"]], -5)
  r2 <- gs_ranks(de, genes = c("A", "A", "B"), collapse = "mean")
  expect_equal(r2[["A"]], -2)
})

test_that("gs_ranks errors clearly on a missing metric column", {
  expect_error(gs_ranks(data.frame(x = 1), metric = "t"), "\"t\"")
})
