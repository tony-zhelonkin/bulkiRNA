# Structural checks on the package namespace itself.
#
# These exist because of a bug that reached the merged branch undetected: an
# exported deprecation shim named `filter_by_size` and an internal helper of the
# same name were both defined at top level. A namespace holds one binding per
# name, R resolves it by collation order, and the internal silently won -- so
# the exported function emitted no deprecation warning and used the wrong
# defaults. `document()`, `test()`, `R CMD check --as-cran` and the golden gate
# were all indifferent to it. Only a human reading a test skip found it.

pkg_r_dir <- function() {
  # Works under both load_all() (source tree) and an installed-package test run.
  for (p in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(p) && length(list.files(p, "[.]R$"))) return(p)
  }
  NULL
}

top_level_functions <- function(dir) {
  out <- list()
  for (f in list.files(dir, "[.]R$", full.names = TRUE)) {
    for (q in parse(f)) {
      is_fn <- is.call(q) && length(q) == 3L &&
        as.character(q[[1]]) %in% c("<-", "=") &&
        is.call(q[[3]]) && identical(as.character(q[[3]][[1]]), "function")
      if (!is_fn) next
      nm <- as.character(q[[2]])
      if (length(nm) == 1L) out[[nm]] <- c(out[[nm]], basename(f))
    }
  }
  out
}

test_that("no function name is defined twice at top level in R/", {
  dir <- pkg_r_dir()
  skip_if(is.null(dir), "package R/ sources not reachable from the test dir")

  fns <- top_level_functions(dir)
  dups <- fns[vapply(fns, length, integer(1L)) > 1L]

  # Reported with the offending files, because "length 0" tells the next
  # person nothing about which two files to look at.
  expect_identical(
    vapply(names(dups), function(n) paste(n, "in", paste(dups[[n]], collapse = " + ")),
           character(1L), USE.NAMES = FALSE),
    character(0L)
  )
})

test_that("no deprecation message sends users to an internal function", {
  dir <- pkg_r_dir()
  skip_if(is.null(dir), "package R/ sources not reachable from the test dir")

  # Renamed to legacy-fixtures-* at 1.0.0. Matching both spellings keeps this
  # from going quietly inert the next time the files move.
  files <- list.files(dir, "^(deprecated|legacy-fixtures)-.*[.]R$", full.names = TRUE)
  expect_true(length(files) > 0L, info = "legacy fixture sources must be present")
  txt <- unlist(lapply(files, readLines, warn = FALSE))
  calls <- grep("\\.Deprecated\\(", txt, value = TRUE)
  expect_true(length(calls) > 0L)

  # Six shims used to name their private successor -- ".gs_plot_all",
  # ".gsdb_human_to_mouse()", "the internal filter_by_size() in R/gs-db.R" --
  # which a user cannot call, and one named a function that no longer existed
  # under that name at all. A deprecation warning's only job is to say what to
  # do instead, so pointing at an unreachable name is worse than silence.
  offenders <- grep('the internal|\\("\\.[a-zA-Z]', calls, value = TRUE)
  expect_identical(offenders, character(0L))
})

test_that("every exported name is actually defined in R/", {
  dir <- pkg_r_dir()
  skip_if(is.null(dir), "package R/ sources not reachable from the test dir")
  ns_file <- file.path(dirname(dir), "NAMESPACE")
  skip_if_not(file.exists(ns_file), "NAMESPACE not reachable")

  ns <- readLines(ns_file, warn = FALSE)
  exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
  defined <- names(top_level_functions(dir))

  # An export naming nothing is silent until a user calls it. S3 methods and
  # data objects are registered by other directives, so only export() counts.
  expect_identical(setdiff(exported, defined), character(0L))
})

# The two checks below widen the guard beyond the exact shape of the bug that
# prompted this file. The duplicate-definition scan only sees
# `name <- function(...)` at top level, and the deprecation scan only sees the
# one phrasing that was found, so a deprecation target that simply does not
# exist, or an S3 method registered for a generic nothing provides, both passed
# silently.

test_that("every identifier a deprecation message names is reachable", {
  dir <- pkg_r_dir()
  skip_if(is.null(dir), "package R/ sources not reachable from the test dir")
  # Renamed to legacy-fixtures-* at 1.0.0. Matching both spellings keeps this
  # from going quietly inert the next time the files move.
  files <- list.files(dir, "^(deprecated|legacy-fixtures)-.*[.]R$", full.names = TRUE)
  expect_true(length(files) > 0L, info = "legacy fixture sources must be present")

  txt <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  # Function-call shapes inside a .Deprecated() message: `foo()`, `pkg::foo()`,
  # `foo(x = )`. Anything without parentheses is prose, not a target.
  msgs <- regmatches(txt, gregexpr("\\.Deprecated\\((?:[^()]|\\([^()]*\\))*\\)", txt))[[1L]]
  named <- unique(unlist(lapply(msgs, function(m) {
    hits <- regmatches(m, gregexpr("[a-zA-Z._][a-zA-Z0-9._]*(?=\\()", m, perl = TRUE))[[1L]]
    setdiff(hits, c(".Deprecated", "sprintf", "paste", "paste0", "c"))
  })))

  defined <- names(top_level_functions(dir))
  # A target may legitimately live in another package (BiocManager::install).
  unreachable <- setdiff(named, c(defined, ls(asNamespace("bulkiRNA"))))
  expect_identical(unreachable, character(0L))
})

test_that("every S3 method registered in NAMESPACE has a reachable generic", {
  dir <- pkg_r_dir()
  skip_if(is.null(dir), "package R/ sources not reachable from the test dir")
  ns_file <- file.path(dirname(dir), "NAMESPACE")
  skip_if_not(file.exists(ns_file), "NAMESPACE not reachable")

  ns <- readLines(ns_file, warn = FALSE)
  entries <- sub("^S3method\\((.*)\\)$", "\\1",
                 grep("^S3method\\(", ns, value = TRUE))
  skip_if(!length(entries), "no S3 methods registered")

  generics <- trimws(vapply(strsplit(entries, ","), `[`, character(1L), 1L))
  generics <- gsub('"', "", generics, fixed = TRUE)
  # A generic may come from this package, from base, or from an import.
  reachable <- vapply(unique(generics), function(g) {
    parts <- strsplit(g, "::", fixed = TRUE)[[1L]]
    if (length(parts) == 2L) {
      # A generic the package re-exports (tibble::tbl_sum comes from pillar)
      # lives in that namespace's imports, not its own frame, so the lookup
      # must not be `inherits = FALSE`.
      requireNamespace(parts[1L], quietly = TRUE) &&
        exists(parts[2L], envir = asNamespace(parts[1L]))
    } else {
      exists(g, envir = asNamespace("bulkiRNA"))
    }
  }, logical(1L))
  expect_identical(names(reachable)[!reachable], character(0L))
})
