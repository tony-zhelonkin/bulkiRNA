test_that("gsdb_list reports the bundled databases", {
  skip_if_not_installed("yaml")
  l <- gsdb_list()
  expect_true(all(c("database", "name", "bundled", "species",
                    "description") %in% names(l)))
  expect_true(all(c("mitopathways", "mitoxplorer", "mito_unified",
                    "transportdb") %in% l$database))
  expect_true(grepl("Mus_musculus",
                    l$species[l$database == "mito_unified"]))
  # gatom is registered but not bundled
  expect_false(l$bundled[l$database == "gatom"])
  expect_identical(l$species[l$database == "gatom"], "")
})

test_that("gsdb_load returns a gs_db for a bundled database", {
  db <- gsdb_load("mito_unified")
  expect_s3_class(db, "gs_db")
  expect_gt(length(db), 100)
  expect_identical(attr(db, "species"), "Mus musculus")
  expect_identical(attr(db, "database"), "mito_unified")  # join key
  expect_identical(attr(db, "database_label"),
                   "Unified Mitochondrial Pathways")  # display only
  expect_true(all(vapply(db, is.character, logical(1))))
  expect_true(all(vapply(db, length, integer(1)) > 0L))
  expect_identical(names(attr(db, "pathway_names")), names(db))
  expect_true(any(grepl("^MITOPATHWAYS_", names(db))))
})

test_that("gsdb_load honours species, size bounds and multiple databases", {
  hs <- gsdb_load("mitopathways", species = "Homo sapiens")
  expect_identical(attr(hs, "species"), "Homo sapiens")
  expect_true(all(toupper(hs[[1]]) == hs[[1]]))  # human symbols

  # underscore form accepted
  expect_identical(attr(gsdb_load("mitopathways", "Mus_musculus"), "species"),
                   "Mus musculus")

  small <- gsdb_load("mitopathways", min_size = 10, max_size = 50)
  sizes <- vapply(small, length, integer(1))
  expect_true(all(sizes >= 10 & sizes <= 50))

  both <- gsdb_load(c("mitopathways", "transportdb"))
  expect_identical(names(both), c("mitopathways", "transportdb"))
  # the list name and the database key agree, so rbind()/gs_filter() join
  expect_identical(vapply(both, attr, character(1), "database"),
                   c(mitopathways = "mitopathways",
                     transportdb = "transportdb"))
  expect_true(all(vapply(both, inherits, logical(1), "gs_db")))
})

test_that("gsdb_load errors are explicit", {
  expect_error(gsdb_load("nope"), "must be one of")
  expect_error(gsdb_load(character()), "non-empty character vector")
  expect_error(gsdb_load("mitoxplorer", species = "Homo sapiens"),
               "No processed")
  expect_error(gsdb_load("mitopathways", rebuild = TRUE),
               "raw reference files are not shipped")
})

test_that("gsdb_info returns metadata and citations", {
  skip_if_not_installed("yaml")
  info <- gsdb_info("mitopathways")
  expect_identical(info$database, "mitopathways")
  expect_identical(info$name, "MitoPathways 3.0")
  expect_true(nzchar(info$citations_path))
  expect_true(length(info$citations_text) > 0)

  object_info <- gsdb_info(gsdb_load("mitopathways"))
  expect_identical(
    object_info[c("database", "name")],
    info[c("database", "name")]
  )

  expect_error(gsdb_info("nope"), "must be one of")
  expect_error(gsdb_info(c("a", "b")), "single database name")
})

test_that("the rebuild driver refuses to run without a source checkout", {
  expect_error(
    bulkiRNA:::.gsdb_rebuild(refs_dir = file.path(tempdir(), "no-such-refs")),
    "raw reference files are not shipped"
  )
})
