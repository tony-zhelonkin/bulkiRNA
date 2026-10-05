# Argument validation only — nothing here touches the network.

test_that("download_gatom_references keeps its frozen formals", {
  f <- formals(download_gatom_references)
  expect_identical(names(f),
                   c("dest_dir", "species", "networks", "overwrite"))
  expect_identical(f$dest_dir, "00_data/references/gatom")
  expect_identical(f$species, "Mus_musculus")
  expect_identical(eval(f$networks), c("kegg", "combined"))
  expect_false(eval(f$overwrite))
})

test_that("download_gatom_references is a deprecated compatibility shim", {
  deprecation <- NULL
  try(
    withCallingHandlers(
      download_gatom_references(networks = character()),
      warning = function(w) {
        if (inherits(w, "deprecatedWarning")) deprecation <<- w
        invokeRestart("muffleWarning")
      }
    ),
    silent = TRUE
  )
  expect_s3_class(deprecation, "deprecatedWarning")
  expect_match(conditionMessage(deprecation), "gatom_download_refs", fixed = TRUE)
  expect_false(grepl("bulkiRNA-deprecated", conditionMessage(deprecation),
                     fixed = TRUE))
  expect_identical(
    names(formals(gatom_download_refs)),
    c("dir", "species", "networks", "overwrite")
  )
})

test_that("download_gatom_references forwards its destination", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  dest <- tempfile("gatom-shim-")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)
  destinations <- character()

  testthat::local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      destinations <<- c(destinations, destfile)
      writeLines("payload", destfile)
      0L
    },
    .package = "utils"
  )

  expect_warning(
    suppressMessages(
      paths <- download_gatom_references(
        dest_dir = dest,
        species = "Homo_sapiens",
        networks = "kegg"
      )
    ),
    class = "deprecatedWarning"
  )
  expect_true(length(destinations) > 0L)
  expect_true(all(dirname(destinations) == dest))
  expect_true(all(dirname(paths) == dest))
})

test_that("an unsupported species fails before anything is written", {
  dest <- file.path(tempdir(), "gatom-unsupported")
  expect_error(
    gatom_download_refs(dir = dest, species = "Xenopus"),
    "`species` must be one of"
  )
  expect_false(dir.exists(dest))
})

test_that("networks must be among the four the server provides", {
  dest <- file.path(tempdir(), "gatom-badnetworks")
  for (bad in list(character(), NA_character_, "reactome")) {
    expect_error(
      gatom_download_refs(dir = dest, networks = bad),
      "`networks` must be a subset of \"kegg\", \"rhea\", \"combined\", \"lipids\"",
      fixed = TRUE
    )
  }
  expect_false(dir.exists(dest))
})

test_that("the downloader fetches the table's files and writes SHA256SUMS", {
  # It used to request gene2reaction.kegg.*.tsv, which the server does not
  # have, so every KEGG download ended in a [FAIL] warning.
  dest <- withr::local_tempdir()
  requested <- character()
  testthat::local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      requested <<- c(requested, basename(url))
      writeLines(basename(url), destfile)
      0L
    },
    .package = "utils"
  )
  expect_no_warning(suppressMessages(
    gatom_download_refs(dir = dest, species = "Mus musculus",
                        networks = c("kegg", "lipids"))
  ))
  expect_setequal(requested, c(
    "network.kegg.rds", "met.kegg.db.rds", "network.rhea.lipids.rds",
    "met.lipids.db.rds", "gene2reaction.rhea.mmu.eg.tsv",
    "org.Mm.eg.gatom.anno.rds"
  ))
  skip_if_not(exists("sha256sum", asNamespace("tools")), "needs R >= 4.5")
  sums <- readLines(file.path(dest, "SHA256SUMS"))
  expect_length(sums, length(requested))
  expect_identical(bulkiRNA:::.gatom_verify_sums(
    dest, c(network = "network.kegg.rds"),
    unname(tools::sha256sum(file.path(dest, "network.kegg.rds")))
  ), "verified")
})

test_that("gsdb_load points GATOM users at the downloader", {
  expect_error(gsdb_load("gatom"),
               "gatom_download_refs", fixed = TRUE)
  expect_error(gsdb_load("GATOM"), "not bundled")
})

# --- regressions from the IO review ------------------------------------------

# The downloader has no URL seam, so the transfer itself is mocked; these tests
# are about the skip/rename logic around it, not about the network.

test_that("an empty leftover file is refetched rather than skipped forever", {
  # An interrupted transfer left a partial file that every later run reported as
  # "[skip] ... (exists)", so gatom_refs() then died in readRDS() naming neither
  # the file nor the fix.
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  dest <- file.path(tempdir(), "gatom-empty")
  unlink(dest, recursive = TRUE)
  dir.create(dest, recursive = TRUE)
  file.create(file.path(dest, "network.kegg.rds"))

  testthat::local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      writeLines("payload", destfile)
      0L
    },
    .package = "utils"
  )
  msgs <- character()
  withCallingHandlers(
    gatom_download_refs(dir = dest, species = "Homo_sapiens",
                        networks = "kegg"),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("present but empty", msgs)))
  expect_false(any(grepl("[skip] network.kegg.rds", msgs, fixed = TRUE)))
  expect_gt(file.info(file.path(dest, "network.kegg.rds"))$size, 0)
})

test_that("a download that produces nothing leaves no file behind", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  dest <- file.path(tempdir(), "gatom-fail")
  unlink(dest, recursive = TRUE)

  # Simulates an interrupted transfer: the .part file appears but is empty.
  testthat::local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      file.create(destfile)
      0L
    },
    .package = "utils"
  )
  suppressWarnings(suppressMessages(
    out <- gatom_download_refs(dir = dest, species = "Homo_sapiens",
                               networks = "kegg")
  ))
  expect_length(out, 0L)
  expect_length(list.files(dest), 0L)
})
