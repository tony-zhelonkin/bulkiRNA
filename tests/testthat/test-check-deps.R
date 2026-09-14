test_that("requested features return one row per package", {
  de <- bulkirna_check_deps("de", quiet = TRUE)

  expect_s3_class(de, "tbl_df")
  expect_identical(de$package, c("edgeR", "limma"))
  expect_true(all(de$feature == "de"))
  expect_identical(
    names(de),
    c("package", "feature", "repository", "installed", "version", "install")
  )
})

test_that("all covers the optional union and excludes development packages", {
  all <- bulkirna_check_deps("all", quiet = TRUE)
  expected <- c(
    "edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
    "babelgene", "biomaRt", "homologene", "GSVA", "gatom", "mwcsr",
    # patchwork became a hard Import at 1.1.0, so it must NOT appear here.
    "igraph", "qs2", "BiocParallel", "plotly", "readxl", "yaml"
  )

  expect_identical(all$package, expected)
  expect_equal(nrow(all), length(unique(all$package)))
  expect_false("testthat" %in% all$package)
})

test_that("multiple features select their union", {
  deps <- bulkirna_check_deps(c("scoring", "io"), quiet = TRUE)
  expect_setequal(deps$package, c("GSVA", "readxl", "yaml"))
  expect_equal(nrow(deps), 3L)
})

test_that("unknown feature names error", {
  expect_error(
    bulkirna_check_deps("not-a-feature", quiet = TRUE),
    "should be one of"
  )
})

test_that("dependency status agrees with installed package state", {
  status <- bulkiRNA:::.bulkirna_dependency_status(c("dplyr", "base"))

  expect_identical(
    status$installed,
    c(requireNamespace("dplyr", quietly = TRUE),
      requireNamespace("base", quietly = TRUE))
  )
  expect_identical(
    status$version,
    c(as.character(utils::packageVersion("dplyr")),
      as.character(utils::packageVersion("base")))
  )
})

test_that("an absent dependency has no version", {
  status <- bulkiRNA:::.bulkirna_dependency_status(
    "definitelyNotAPackageBulkiRNA"
  )

  expect_false(status$installed)
  expect_true(is.na(status$version))
})

test_that("error mode stops when a requested dependency is missing", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  testthat::local_mocked_bindings(
    .bulkirna_optional_deps = function() {
      data.frame(
        package = "definitelyNotAPackageBulkiRNA",
        feature = "de",
        repository = "CRAN",
        install = 'install.packages("definitelyNotAPackageBulkiRNA")',
        stringsAsFactors = FALSE
      )
    },
    .package = "bulkiRNA"
  )

  expect_error(
    bulkirna_check_deps("de", quiet = TRUE, error = TRUE),
    "Missing optional package.*definitelyNotAPackageBulkiRNA"
  )
})

test_that("error mode does not stop when requested dependencies are present", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  testthat::local_mocked_bindings(
    .bulkirna_optional_deps = function() {
      data.frame(
        package = "dplyr",
        feature = "de",
        repository = "CRAN",
        install = 'install.packages("dplyr")',
        stringsAsFactors = FALSE
      )
    },
    .package = "bulkiRNA"
  )

  expect_s3_class(
    bulkirna_check_deps("de", quiet = TRUE, error = TRUE),
    "tbl_df"
  )
})

test_that("Bioconductor install commands name BiocManager", {
  deps <- bulkirna_check_deps("annotation", quiet = TRUE)
  bioc <- deps$repository == "Bioconductor"

  expect_true(any(bioc))
  expect_true(all(grepl("BiocManager", deps$install[bioc], fixed = TRUE)))
})

test_that("network repositories match their standard sources", {
  deps <- bulkirna_check_deps("network", quiet = TRUE)

  expect_identical(
    deps$repository,
    c("Bioconductor", "CRAN", "CRAN")
  )
  expect_match(deps$install[deps$package == "gatom"], "BiocManager")
  expect_match(deps$install[deps$package == "mwcsr"], "install.packages")
})

test_that("printing makes the returned report invisible", {
  expect_invisible(bulkirna_check_deps("de", quiet = FALSE))
})

test_that("the registry cannot silently drift from DESCRIPTION's Suggests", {
  # ADR-003's premise is that Suggests is the machine-readable manifest. A
  # hardcoded registry that diverges from it defeats exactly that, and the
  # divergence would be invisible: a new Suggests entry simply never appears
  # in the report. This test is the seam.
  suggests <- utils::packageDescription("bulkiRNA", fields = "Suggests")
  skip_if(is.na(suggests))

  declared <- trimws(strsplit(suggests, ",", fixed = TRUE)[[1]])
  declared <- sub("\\s*\\(.*\\)$", "", declared)   # drop version constraints
  declared <- declared[nzchar(declared)]

  # Never user-facing: needed to test or to build this package's own
  # documentation, not to use any feature. `bulkirna_check_deps()` answers
  # "what must I install to use this?", and the answer never includes knitr.
  # colorspace simulates colour-vision deficiency in test-colors.R only:
  # no runtime path reaches it, so it is dev tooling, not an optional feature.
  dev_only <- c("testthat", "knitr", "rmarkdown", "withr", "colorspace")
  expected <- setdiff(declared, dev_only)
  reported <- bulkirna_check_deps("all", quiet = TRUE)$package

  expect_setequal(reported, expected)
})

test_that("the coresh feature names the packages the chunk reader needs", {
  deps <- bulkirna_check_deps("coresh", quiet = TRUE)
  expect_identical(deps$package, c("qs2", "BiocParallel"))
  expect_identical(deps$repository, c("CRAN", "Bioconductor"))
})
