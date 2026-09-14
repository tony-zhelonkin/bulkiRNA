test_that("pinned seeding restores an absent RNG seed", {
  local_pinned_rng()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  kind_before <- RNGkind()

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  value <- .with_pinned_seed(42L, rnorm(2L))

  expect_length(value, 2L)
  expect_identical(RNGkind(), kind_before)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("pinned seeding restores an existing RNG state byte-for-byte", {
  local_pinned_rng()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(919L)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  .with_pinned_seed(42L, runif(3L))

  expect_identical(RNGkind(), kind_before)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

test_that("pinned seeding restores RNG state when evaluation fails", {
  local_pinned_rng()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(818L)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  expect_error(
    .with_pinned_seed(42L, {
      runif(1L)
      stop("seeded failure", call. = FALSE)
    }),
    "seeded failure"
  )

  expect_identical(RNGkind(), kind_before)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

test_that("a NULL pinned seed evaluates without touching RNG state", {
  local_pinned_rng()
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(717L)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L

  value <- .with_pinned_seed(NULL, {
    counter$n <- counter$n + 1L
    "evaluated"
  })

  expect_identical(value, "evaluated")
  expect_identical(counter$n, 1L)
  expect_identical(RNGkind(), kind_before)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

test_that("a pinned expression is evaluated exactly once", {
  local_pinned_rng()
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L

  value <- .with_pinned_seed(42L, {
    counter$n <- counter$n + 1L
    counter$n
  })

  expect_identical(value, 1L)
  expect_identical(counter$n, 1L)
})

test_that("pinned seeding preserves default-generator draws", {
  local_pinned_rng()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  expected <- {
    set.seed(42)
    rnorm(5L)
  }

  # This equivalence guarantees that existing GATOM and GSEA results do not
  # move.
  actual <- .with_pinned_seed(42, rnorm(5L))

  expect_identical(actual, expected)
})
