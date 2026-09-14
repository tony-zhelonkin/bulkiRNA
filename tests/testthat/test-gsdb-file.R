gmt_lines <- c(
  "SET_A\tfirst set\tActb\tGapdh\tSdha",
  "SET_B\tna\tSdha\tNdufa1",
  "SET_C\t\tOnly1"
)

gmx_lines <- c(
  "first set\tsecond set",
  "SET_A\tSET_B",
  "Actb\tSdha",
  "Gapdh\tNdufa1",
  "Sdha\t"
)

test_that("gsdb_from_file reads GMT and derives labels", {
  f <- tempfile(fileext = ".gmt")
  writeLines(gmt_lines, f)
  db <- gsdb_from_file(f, species = "Mus musculus")
  expect_s3_class(db, "gs_db")
  expect_identical(names(db), c("SET_A", "SET_B", "SET_C"))
  expect_identical(db$SET_A, c("Actb", "Gapdh", "Sdha"))
  labs <- attr(db, "pathway_names")
  expect_identical(labs[["SET_A"]], "first set")
  expect_identical(labs[["SET_B"]], "SET_B")  # "na" -> id
  expect_identical(labs[["SET_C"]], "SET_C")  # empty -> id
  expect_identical(attr(db, "database"), bulkiRNA:::.gsdb_key_from_path(f))
  expect_false(grepl("[^A-Za-z0-9_]", attr(db, "database")))  # typeable
  expect_identical(attr(db, "database_label"), basename(f))
})

test_that("gsdb_from_file reads GMX with ragged columns", {
  f <- tempfile(fileext = ".gmx")
  writeLines(gmx_lines, f)
  db <- gsdb_from_file(f, species = "Mus musculus", database = "mine")
  expect_identical(names(db), c("SET_A", "SET_B"))
  expect_identical(db$SET_A, c("Actb", "Gapdh", "Sdha"))
  expect_identical(db$SET_B, c("Sdha", "Ndufa1"))
  expect_identical(attr(db, "pathway_names")[["SET_B"]], "second set")
  expect_identical(attr(db, "database"), "mine")
  # `database_label` must not silently become the machine key: it defaults
  # to the file's base name, independent of `database` (see gs-db.R:9-19).
  expect_identical(attr(db, "database_label"), basename(f))
})

test_that("gsdb_from_file's database_label is independent of database", {
  f <- tempfile(fileext = ".gmt")
  writeLines(gmt_lines, f)
  db <- gsdb_from_file(f, database = "mine", database_label = "My Sets")
  expect_identical(attr(db, "database"), "mine")
  expect_identical(attr(db, "database_label"), "My Sets")
})

test_that("the format is sniffed when the extension does not say", {
  gmt <- tempfile(fileext = ".txt")
  writeLines(gmt_lines, gmt)
  expect_identical(bulkiRNA:::.gsdb_sniff_format(gmt, gmt_lines), "gmt")
  expect_identical(names(gsdb_from_file(gmt)), c("SET_A", "SET_B", "SET_C"))

  gmx <- tempfile(fileext = ".txt")
  writeLines(gmx_lines, gmx)
  expect_identical(bulkiRNA:::.gsdb_sniff_format(gmx, gmx_lines), "gmx")
  expect_identical(names(gsdb_from_file(gmx)), c("SET_A", "SET_B"))
})

test_that("a uniform-width GMT is not mistaken for GMX just by extension-less sniffing", {
  # Every "row" here has the same number of tab fields (name, desc, gene,
  # gene), which used to satisfy the GMX width heuristic on its own and get
  # parsed as GMX -- turning ids into gene symbols and the shared "desc"
  # description into a bogus 4th set.
  uniform <- c("S1\tdesc\tg1\tg2", "S2\tdesc\tg3\tg4", "S3\tdesc\tg5\tg6")
  f <- tempfile(fileext = ".txt")
  writeLines(uniform, f)
  expect_identical(bulkiRNA:::.gsdb_sniff_format(f, uniform), "gmt")
  expect_identical(names(gsdb_from_file(f)), c("S1", "S2", "S3"))
})

test_that("gsdb_from_file applies prefix and size bounds", {
  f <- tempfile(fileext = ".gmt")
  writeLines(gmt_lines, f)
  db <- gsdb_from_file(f, prefix = "MY", min_size = 2)
  expect_identical(names(db), c("MY_SET_A", "MY_SET_B"))
  expect_identical(attr(db, "pathway_names")[["MY_SET_A"]], "first set")
})

test_that("gsdb_from_file errors on bad paths and empty files", {
  expect_error(gsdb_from_file(file.path(tempdir(), "nope.gmt")),
               "does not exist")
  empty <- tempfile(fileext = ".gmt")
  writeLines(character(), empty)
  expect_error(gsdb_from_file(empty), "is empty")
  short <- tempfile(fileext = ".gmx")
  writeLines(c("a", "b"), short)
  expect_error(gsdb_from_file(short), "at least three rows")
  nogmt <- tempfile(fileext = ".gmt")
  writeLines(c("A\tB", "C\tD"), nogmt)
  expect_error(gsdb_from_file(nogmt), "No GMT records")
})

test_that("gsdb_register promotes a plain list to a gs_db", {
  db <- gsdb_register(
    list(MY_SET = c("Actb", "Gapdh"), OTHER = c("Sdha", "Ndufa1", "Sdhb")),
    database = "my_signatures", species = "Mus musculus",
    pathway_names = c(MY_SET = "My favourite set"),
    database_label = "My signatures"
  )
  expect_s3_class(db, "gs_db")
  expect_identical(attr(db, "database"), "my_signatures")
  expect_identical(attr(db, "database_label"), "My signatures")
  expect_identical(attr(db, "pathway_names")[["MY_SET"]], "My favourite set")
  expect_identical(names(gsdb_register(list(A = "a", B = c("b", "c")),
                                       "d", min_size = 2)), "B")
  expect_error(gsdb_register("not a list", "d"), "named list")
  expect_error(gsdb_register(list(A = 1:3), "d"), "character vectors")
})

test_that(".gsdb_key_from_path yields a typeable key", {
  k <- bulkiRNA:::.gsdb_key_from_path
  expect_identical(k("/tmp/MitoPathways3.0.gmx"), "MitoPathways3_0")
  expect_identical(k("my sets (v2).gmt"), "my_sets_v2")
  expect_identical(k("/tmp/---.gmt"), "genesets")
})
