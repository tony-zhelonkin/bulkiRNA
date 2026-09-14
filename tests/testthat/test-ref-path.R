test_that(".ref_path resolves current and records its snapshot", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot <- file.path(source_dir, "coresh_snapshot_test")
  dir.create(snapshot, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  withr::local_envvar(c(REFCACHE_ROOT = root))

  expect_no_warning(out <- .ref_path("coresh"))

  expect_identical(as.vector(out), file.path(source_dir, "current"))
  expect_identical(attr(out, "snapshot"), basename(snapshot))
  expect_identical(attr(out, "source"), "coresh")
  expect_identical(attr(out, "root"), root)
  expect_false(attr(out, "caller_supplied"))
})

test_that(".ref_path resolves a relative current target", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot_tag <- "syn66227307_20260721"
  dir.create(file.path(source_dir, snapshot_tag), recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot_tag, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }

  expect_no_warning(out <- .ref_path("coresh", root = root))

  expect_identical(attr(out, "snapshot"), snapshot_tag)
})

test_that(".ref_path accepts a nested snapshot inside its source", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot_tag <- "syn66227307_20260721"
  snapshot <- file.path(source_dir, "snapshots", snapshot_tag)
  dir.create(snapshot, recursive = TRUE)
  linked <- suppressWarnings(file.symlink(
    file.path("snapshots", snapshot_tag),
    file.path(source_dir, "current")
  ))
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }

  expect_no_warning(out <- .ref_path("coresh", root = root))

  expect_identical(attr(out, "snapshot"), snapshot_tag)
})

test_that(".ref_path warns once when current resolves outside its source", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  external <- file.path(root, "mounted", "syn66227307_20260721")
  dir.create(source_dir, recursive = TRUE)
  dir.create(external, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(external, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }

  warnings <- character()
  first <- withCallingHandlers(
    .ref_path("coresh", root = root),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(attr(first, "snapshot"), basename(external))
  second <- withCallingHandlers(
    .ref_path("coresh", root = root),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(
    warnings[[1L]],
    normalizePath(external, winslash = "/", mustWork = TRUE),
    fixed = TRUE
  )
  expect_match(
    warnings[[1L]],
    "may not identify a refcache snapshot",
    fixed = TRUE
  )
  expect_identical(attr(second, "snapshot"), basename(external))
})

test_that("clearing reference resolutions re-arms the layout warning", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  external <- file.path(root, "mounted", "syn66227307_20260721")
  dir.create(source_dir, recursive = TRUE)
  dir.create(external, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(external, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }

  expect_warning(.ref_path("coresh", root = root), "outside")
  expect_no_warning(.ref_path("coresh", root = root))
  .clear_ref_resolutions()
  expect_warning(.ref_path("coresh", root = root), "outside")
})

test_that(".ref_path appends path components below current", {
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  snapshot <- file.path(source_dir, "coresh_snapshot_test")
  leaf <- file.path(snapshot, "chunks", "part-01")
  dir.create(leaf, recursive = TRUE)
  linked <- suppressWarnings(
    file.symlink(snapshot, file.path(source_dir, "current"))
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  withr::local_envvar(c(REFCACHE_ROOT = root))

  out <- .ref_path("coresh", "chunks", "part-01")

  expect_identical(
    as.vector(out),
    file.path(source_dir, "current", "chunks", "part-01")
  )
  expect_identical(attr(out, "snapshot"), basename(snapshot))
})

test_that("an explicit path wins and is marked caller-supplied", {
  root <- withr::local_tempdir()
  explicit <- file.path(root, "pinned-copy")
  dir.create(explicit)
  withr::local_envvar(c(REFCACHE_ROOT = file.path(root, "unused")))

  out <- .ref_path("coresh", "ignored-component", path = explicit)

  expect_identical(as.vector(out), explicit)
  expect_identical(attr(out, "snapshot"), NA_character_)
  expect_identical(attr(out, "root"), NA_character_)
  expect_true(attr(out, "caller_supplied"))
})

test_that("an unset REFCACHE_ROOT names both remedies", {
  withr::local_envvar(c(REFCACHE_ROOT = NA_character_))

  expect_error(
    .ref_path("coresh"),
    "REFCACHE_ROOT.*path="
  )
})

test_that("a missing source and missing current have different errors", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(REFCACHE_ROOT = root))

  expect_error(.ref_path("coresh"), "source directory is missing")

  source_dir <- file.path(root, "coresh")
  dir.create(source_dir)
  expect_error(.ref_path("coresh"), "current.*missing or dangling")
})

test_that("a dangling current identifies a failed refresh", {
  root <- withr::local_tempdir()
  source_dir <- file.path(root, "coresh")
  dir.create(source_dir)
  withr::local_envvar(c(REFCACHE_ROOT = root))

  linked <- suppressWarnings(
    file.symlink(
      file.path(source_dir, "absent-snapshot"),
      file.path(source_dir, "current")
    )
  )
  if (!isTRUE(linked)) {
    skip("This filesystem does not support symbolic links")
  }
  expect_error(.ref_path("coresh"), "refcache refresh may have failed")
})

test_that("a real-directory current warns once", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  current <- file.path(root, "coresh", "current")
  dir.create(current, recursive = TRUE)
  withr::local_envvar(c(REFCACHE_ROOT = root))

  warnings <- character()
  out <- withCallingHandlers(
    .ref_path("coresh"),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  withCallingHandlers(
    .ref_path("coresh"),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings, "is not a symlink", fixed = TRUE)
  expect_match(warnings, "snapshot cannot be identified", fixed = TRUE)
  expect_identical(as.vector(out), current)
  expect_identical(attr(out, "snapshot"), basename(current))
})

test_that("source must be a single safe path segment", {
  expect_error(.ref_path("../etc", path = tempdir()), "single path segment")
})

test_that("the session record keeps one row per source and clears", {
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  root <- withr::local_tempdir()
  first <- file.path(root, "first")
  second <- file.path(root, "second")
  other <- file.path(root, "other")
  dir.create(first)
  dir.create(second)
  dir.create(other)

  .ref_path("coresh", path = first)
  .ref_path("coresh", path = second)
  .ref_path("cistarget", path = other)
  record <- .ref_resolutions()

  expect_identical(record$source, c("cistarget", "coresh"))
  expect_equal(nrow(record), 2L)
  expect_identical(record$path[record$source == "coresh"], second)

  .clear_ref_resolutions()
  expect_equal(nrow(.ref_resolutions()), 0L)
})
