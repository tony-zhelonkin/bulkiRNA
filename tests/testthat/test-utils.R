test_that("ensure_dir creates directories recursively and is idempotent", {
  d <- file.path(tempdir(), "bulkirna-utils", "a", "b")
  on.exit(unlink(file.path(tempdir(), "bulkirna-utils"), recursive = TRUE))
  expect_identical(ensure_dir(d), d)
  expect_true(dir.exists(d))
  expect_silent(ensure_dir(d))
})

test_that("ensure_dir is vectorized over directory paths", {
  root <- file.path(tempdir(), "bulkirna-utils-vec")
  on.exit(unlink(root, recursive = TRUE))
  ds <- file.path(root, c("x", "y"))
  ensure_dir(ds)
  expect_true(all(dir.exists(ds)))
})

test_that("ensure_dir rejects empty, NULL and NA paths", {
  expect_error(ensure_dir(NULL), "non-empty")
  expect_error(ensure_dir(character(0)), "non-empty")
  expect_error(ensure_dir(NA_character_), "non-missing")
  expect_error(ensure_dir(""), "non-missing")
})

test_that("ensure_dir errors (rather than silently returning) when creation fails", {
  skip_on_os("windows")
  parent <- file.path(tempdir(), "bulkirna-utils-readonly")
  unlink(parent, recursive = TRUE)
  dir.create(parent)
  on.exit({
    Sys.chmod(parent, "0755")
    unlink(parent, recursive = TRUE)
  })
  Sys.chmod(parent, "0500")
  target <- file.path(parent, "child")
  expect_error(ensure_dir(target), "Failed to create directory")
  expect_false(dir.exists(target))
})

test_that("ensure_parent_dir takes a file path and creates its parent only", {
  root <- file.path(tempdir(), "bulkirna-utils-parent")
  on.exit(unlink(root, recursive = TRUE))
  f <- file.path(root, "figures", "plot.pdf")
  bulkiRNA:::ensure_parent_dir(f)
  expect_true(dir.exists(dirname(f)))
  expect_false(file.exists(f))
  expect_false(dir.exists(f))
})

test_that("%||% returns the default only for NULL", {
  expect_identical(bulkiRNA:::`%||%`(NULL, 1), 1)
  expect_identical(bulkiRNA:::`%||%`(2, 1), 2)
  expect_identical(bulkiRNA:::`%||%`(NA, 1), NA)
})
