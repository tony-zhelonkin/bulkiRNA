test_that("gs_db carries the frozen contract", {
  db <- bulkiRNA:::gs_db(
    list(A = c("Actb", "Actb", "Gapdh"), B = c("Sdha")),
    database = "testdb", species = "Mus_musculus",
    pathway_names = c(A = "Set A")
  )
  expect_s3_class(db, "gs_db")
  expect_true(is.list(db))
  expect_identical(names(db), c("A", "B"))
  expect_identical(db$A, c("Actb", "Gapdh"))  # de-duplicated
  expect_identical(attr(db, "database"), "testdb")
  expect_identical(attr(db, "species"), "Mus musculus")  # normalised
  expect_identical(attr(db, "gene_id_type"), "symbol")
  expect_identical(attr(db, "pathway_names"),
                   c(A = "Set A", B = "B"))
  expect_identical(names(attr(db, "pathway_names")), names(db))
})

test_that("database is the join key and database_label is display-only", {
  db <- bulkiRNA:::gs_db(list(A = "a"), database = "mito_unified",
                         species = "Mus musculus",
                         database_label = "Unified Mitochondrial Pathways")
  expect_identical(attr(db, "database"), "mito_unified")
  expect_identical(attr(db, "database_label"),
                   "Unified Mitochondrial Pathways")

  # label defaults to the key, so the attribute is never absent
  plain <- bulkiRNA:::gs_db(list(A = "a"), database = "k",
                            species = "Mus musculus")
  expect_identical(attr(plain, "database_label"), "k")
  expect_output(print(plain), "<gs_db> k", fixed = TRUE)
  expect_output(print(db), "[mito_unified]", fixed = TRUE)

  # both survive subsetting
  expect_identical(attr(db["A"], "database"), "mito_unified")
  expect_identical(attr(db["A"], "database_label"),
                   "Unified Mitochondrial Pathways")

  expect_error(bulkiRNA:::gs_db(list(A = "a"), database = "k",
                                species = "Mus musculus",
                                database_label = c("a", "b")),
               "`database_label` must be")
})

test_that("a filter_by_size(result, min_size, max_size) shim is writable", {
  # Step C must be able to reproduce the old frozen signature, which took and
  # returned the {T2G, T2N} list, via from_t2g -> filter -> as_t2g.
  legacy <- list(
    T2G = data.frame(
      gs_name = c(rep("SMALL", 3), rep("OK", 6), rep("BIG", 12)),
      gene_symbol = paste0("g", seq_len(21)),
      stringsAsFactors = FALSE
    ),
    T2N = data.frame(gs_name = c("SMALL", "OK", "BIG"),
                     description = c("s", "o", "b"),
                     stringsAsFactors = FALSE)
  )
  shim <- function(result, min_size = 5, max_size = 500) {
    db <- bulkiRNA:::.gsdb_from_t2g(result, "legacy", "Mus musculus")
    bulkiRNA:::.gsdb_as_t2g(bulkiRNA:::.gs_filter_size(db, min_size, max_size))
  }
  out <- shim(legacy)                      # old default min_size = 5
  expect_identical(sort(unique(out$T2G$gs_name)), c("BIG", "OK"))
  expect_identical(sort(out$T2N$gs_name), c("BIG", "OK"))
  expect_identical(names(out$T2G), c("gs_name", "gene_symbol"))
  expect_identical(names(out$T2N), c("gs_name", "description"))
  expect_identical(sort(unique(shim(legacy, 1, 5)$T2G$gs_name)), "SMALL")
})

test_that("gs_db drops empty sets and rejects bad input", {
  db <- bulkiRNA:::gs_db(list(A = "Actb", B = character(), C = c(NA, "")),
                         database = "d", species = "Mus musculus")
  expect_identical(names(db), "A")

  expect_error(bulkiRNA:::gs_db(list(1, 2), database = "d",
                                species = "Mus musculus"),
               "named list")
  expect_error(bulkiRNA:::gs_db(list(A = "x", A = "y"), database = "d",
                                species = "Mus musculus"),
               "unique")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "",
                                species = "Mus musculus"),
               "non-empty string")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "d", species = NA),
               "`species` must be")
  expect_error(bulkiRNA:::gs_db(list(A = "x"), database = "d",
                                species = "Mus musculus",
                                gene_id_type = "entrez"),
               "gene_id_type")
})

test_that("subsetting and summary keep the contract", {
  db <- bulkiRNA:::gs_db(list(A = "a", B = c("b", "c")),
                         database = "d", species = "Mus musculus")
  sub <- db["B"]
  expect_s3_class(sub, "gs_db")
  expect_identical(names(sub), "B")
  expect_identical(names(attr(sub, "pathway_names")), "B")
  expect_identical(attr(sub, "database"), "d")

  s <- summary(db)
  expect_identical(names(s), c("pathway_id", "pathway_name", "n_genes"))
  expect_identical(s$n_genes, c(1L, 2L))
  expect_output(print(db), "<gs_db>")
})

test_that("pathway_descriptions survives subsetting instead of vanishing", {
  db <- bulkiRNA:::gs_db(
    list(A = "a", B = c("b", "c")),
    database = "d", species = "Mus musculus",
    pathway_descriptions = c(A = "desc A", B = "desc B")
  )
  expect_identical(attr(db, "pathway_descriptions"), c(A = "desc A", B = "desc B"))

  sub <- db["A"]
  expect_identical(attr(sub, "pathway_descriptions"), c(A = "desc A"))

  # ids no longer present are dropped, not left dangling
  expect_null(attr(bulkiRNA:::gs_db(list(A = "a"), database = "d",
                                    species = "Mus musculus",
                                    pathway_descriptions = c(Z = "unused")),
                   "pathway_descriptions"))
})

test_that("database and set provenance survive subsetting", {
  set_provenance <- tibble::tibble(
    set_name = c("A", "B", "C"),
    source_row = c(11L, 22L, 33L)
  )
  provenance <- list(
    source = "test-snapshot",
    snapshot = "snapshot-20260818"
  )
  db <- bulkiRNA:::gs_db(
    list(A = "a", B = "b", C = "c"),
    database = "d",
    species = "Mus musculus",
    set_provenance = set_provenance,
    provenance = provenance
  )

  sub <- db["B"]
  expect_identical(attr(sub, "provenance"), provenance)
  expect_identical(attr(sub, "set_provenance")$set_name, "B")

  info <- gsdb_info(sub)
  expect_identical(info$name, "d")
  expect_identical(info$provenance, provenance)
  expect_identical(info$set_provenance$set_name, "B")
  printed <- capture.output(print(sub))
  expect_identical(sum(grepl("^Provenance:", printed)), 1L)
  expect_true(any(grepl("snapshot=snapshot-20260818", printed, fixed = TRUE)))
})

test_that("long database provenance wraps at the configured display width", {
  db <- bulkiRNA:::gs_db(
    list(A = "a"),
    database = "d",
    species = "Mus musculus",
    provenance = list(
      source = "coresh",
      snapshot = "syn66227307_20260721",
      queries = "iron_uptake(5), ferroptosis(4), heme(4)",
      top_hits = 5L,
      top_n = 50L
    )
  )
  withr::local_options(width = 60L)

  printed <- capture.output(print(db))
  provenance_lines <- grep("^(Provenance: | {12})", printed, value = TRUE)

  expect_gt(length(provenance_lines), 1L)
  expect_true(all(nchar(provenance_lines, type = "width") <= 60L))
  expect_true(grepl("^Provenance: source=coresh", provenance_lines[[1L]]))
  expect_true(all(grepl("^ {12}", provenance_lines[-1L])))
})

test_that("set provenance is restricted and reordered to retained sets", {
  set_provenance <- tibble::tibble(
    set_name = c("A", "B", "C"),
    source_row = c(11L, 22L, 33L)
  )
  db <- bulkiRNA:::gs_db(
    list(A = "a", B = "b", C = "c"),
    database = "d",
    species = "Mus musculus",
    set_provenance = set_provenance,
    provenance = list(snapshot = "snapshot-20260818")
  )

  sub <- db[c("C", "A")]
  expect_identical(
    attr(sub, "set_provenance"),
    set_provenance[c(3L, 1L), ]
  )
  expect_identical(attr(sub, "set_provenance")$set_name, c("C", "A"))
})

test_that("gs_db rejects provenance that cannot be read or keyed", {
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a"), database = "d", species = "Mus musculus",
      provenance = c("snapshot", "current")
    ),
    "one-row data frame or a named list"
  )
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a"), database = "d", species = "Mus musculus",
      provenance = list(source = c("one", "two"))
    ),
    "scalar values"
  )
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a"), database = "d", species = "Mus musculus",
      provenance = data.frame(source = I(list("nested")))
    ),
    "atomic scalar value"
  )
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a"), database = "d", species = "Mus musculus",
      provenance = data.frame(
        source = I(matrix(c("one", "two"), nrow = 1L))
      )
    ),
    "atomic scalar value"
  )
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a", B = "b"), database = "d", species = "Mus musculus",
      set_provenance = data.frame(set_name = "A")
    ),
    "must match the retained set names"
  )
  expect_error(
    bulkiRNA:::gs_db(
      list(A = "a"), database = "d", species = "Mus musculus",
      set_provenance = data.frame(set_name = c("A", "A"))
    ),
    "must contain unique"
  )
})

test_that("subsetting a gs_db by an unknown name errors about the index, not `sets`", {
  db <- bulkiRNA:::gs_db(list(A = "a"), database = "d", species = "Mus musculus")
  expect_error(db["nope"], "not in this database")
  expect_error(db["nope"], "nope")
  expect_error(db[c("A", "nope")], "nope")
  expect_identical(names(db["A"]), "A")  # valid names still work
})

test_that(".gs_filter_size honours open bounds", {
  db <- bulkiRNA:::gs_db(
    list(A = "a", B = c("b", "c"), C = c("d", "e", "f")),
    database = "d", species = "Mus musculus"
  )
  expect_identical(names(bulkiRNA:::.gs_filter_size(db)), c("A", "B", "C"))
  expect_identical(names(bulkiRNA:::.gs_filter_size(db, min_size = 2)),
                   c("B", "C"))
  expect_identical(names(bulkiRNA:::.gs_filter_size(db, max_size = 2)),
                   c("A", "B"))
  expect_identical(names(bulkiRNA:::.gs_filter_size(db, 2, 2)), "B")
  expect_error(bulkiRNA:::.gs_filter_size(db, 5, 2), "must not exceed")
  expect_error(bulkiRNA:::.gs_filter_size(list(A = "a"), 1), "must be a `gs_db`")
})

test_that(".gsdb_as_t2g round-trips the legacy shape", {
  db <- bulkiRNA:::gs_db(list(A = c("a", "b"), B = "c"),
                         database = "d", species = "Mus musculus",
                         pathway_names = c(A = "Set A", B = "Set B"))
  legacy <- bulkiRNA:::.gsdb_as_t2g(db)
  expect_identical(names(legacy), c("T2G", "T2N"))
  expect_identical(names(legacy$T2G), c("gs_name", "gene_symbol"))
  expect_identical(names(legacy$T2N), c("gs_name", "description"))
  expect_identical(nrow(legacy$T2G), 3L)
  expect_identical(legacy$T2N$description, c("Set A", "Set B"))

  back <- bulkiRNA:::.gsdb_from_t2g(legacy, database = "d",
                                    species = "Mus musculus")
  expect_identical(unclass(back)[names(db)], unclass(db)[names(db)])
  expect_identical(attr(back, "pathway_names"), attr(db, "pathway_names"))
})

test_that(".gsdb_from_t2g validates its input", {
  expect_error(bulkiRNA:::.gsdb_from_t2g(list(), "d", "Mus musculus"), "T2G")
  expect_error(
    bulkiRNA:::.gsdb_from_t2g(list(T2G = data.frame(a = 1)), "d",
                              "Mus musculus"),
    "gene_symbol"
  )
})
