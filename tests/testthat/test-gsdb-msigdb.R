test_that("gsdb_msigdb validates arguments before touching the network", {
  expect_error(gsdb_msigdb(collection = character()), "`collection` must be")
  expect_error(gsdb_msigdb(species = NA), "`species` must be")
  expect_error(gsdb_msigdb(db_species = "rat"), "'arg' should be one of")
})

test_that("gsdb_msigdb returns a gs_db for the hallmarks", {
  skip_if_not_installed("msigdbr")
  db <- tryCatch(
    gsdb_msigdb("Mus musculus", collection = "H"),
    error = function(e) skip(paste("MSigDB unavailable:", conditionMessage(e)))
  )
  expect_s3_class(db, "gs_db")
  expect_identical(attr(db, "database"), "msigdb_H")  # join key
  expect_identical(attr(db, "database_label"), "MSigDB H")
  expect_identical(attr(db, "species"), "Mus musculus")
  expect_identical(attr(db, "gene_id_type"), "symbol")
  expect_equal(length(db), 50)
  expect_true(all(grepl("^HALLMARK_", names(db))))
  expect_identical(names(attr(db, "pathway_names")), names(db))
  expect_true(all(vapply(db, is.character, logical(1))))

  sub <- gsdb_msigdb("Mus musculus", collection = "H", min_size = 200)
  expect_true(all(vapply(sub, length, integer(1)) >= 200))
  expect_lt(length(sub), length(db))
})

test_that("pathway_descriptions is attached and survives subsetting", {
  skip_if_not_installed("msigdbr")
  db <- tryCatch(
    gsdb_msigdb("Mus musculus", collection = "H"),
    error = function(e) skip(paste("MSigDB unavailable:", conditionMessage(e)))
  )
  desc <- attr(db, "pathway_descriptions")
  expect_true(!is.null(desc))
  expect_identical(sort(names(desc)), sort(names(db)))
  expect_true(all(nzchar(desc)))

  one <- names(db)[1]
  sub <- db[one]
  expect_identical(attr(sub, "pathway_descriptions"), desc[one])
})

test_that("a collection is independent of what was queried before it", {
  # Regression test for an upstream msigdbr bug (confirmed 26.1.0): the
  # ortholog table is built from the current collection's genes but cached
  # under a collection-independent key, so without the workaround in
  # `.msigdbr_drop_ortholog_cache()` a Hallmark query first shrinks Reactome
  # from 1839 sets / 10762 symbols to 1817 / 3688. It is not an error -- just
  # quietly under-tested sets and a wrong BH family -- so only a test that
  # compares the two call orders can catch it.
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("babelgene")
  # `gsdb_msigdb()` now clears the cache on every call, so poisoning has to go
  # through msigdbr directly -- and the baseline has to be taken from a cleared
  # cache, because any earlier test in this file may already have poisoned it.
  # (An earlier version of this test called `gsdb_msigdb()` for both arms and
  # was vacuous: both were equally truncated, so it passed with the fix
  # disabled.)
  raw <- function(coll, sub = NULL) {
    a <- list(db_species = "HS", species = "Mus musculus", collection = coll)
    if (!is.null(sub)) a$subcollection <- sub
    tryCatch(do.call(msigdbr::msigdbr, a),
             error = function(e) skip(paste("MSigDB unavailable:",
                                            conditionMessage(e))))
  }
  symbols <- function(x) sort(unique(as.character(x)))

  .msigdbr_drop_ortholog_cache()
  truth <- symbols(raw("C2", "CP:REACTOME")$gene_symbol)

  # Poisoning only bites from a *cleared* cache: msigdbr keeps whatever table is
  # already cached, so a Hallmark call on top of Reactome's (superset) table
  # reuses it and changes nothing. The small collection has to go first.
  poison <- function() {
    .msigdbr_drop_ortholog_cache()
    invisible(raw("H"))
  }

  poison()                                  # bypasses the workaround
  poisoned <- symbols(raw("C2", "CP:REACTOME")$gene_symbol)
  # If upstream ever fixes this, the workaround becomes a harmless no-op and
  # there is nothing left to regress against -- skip rather than fail.
  skip_if(identical(poisoned, truth),
          "msigdbr no longer leaks its ortholog cache across collections")

  poison()                                  # poison again, then go through us
  db <- gsdb_msigdb("Mus musculus", collection = "C2",
                    subcollection = "CP:REACTOME", db_species = "HS")
  expect_identical(symbols(unlist(db, use.names = FALSE)), truth)
  expect_gt(length(truth), length(poisoned))
})

test_that("dropping the ortholog cache is safe when there is nothing to drop", {
  skip_if_not_installed("msigdbr")
  # Idempotent, and never an error even on a cold cache or a changed upstream.
  expect_type(bulkiRNA:::.msigdbr_drop_ortholog_cache(), "logical")
  expect_false(bulkiRNA:::.msigdbr_drop_ortholog_cache())
})

test_that("an unknown collection errors with guidance", {
  skip_if_not_installed("msigdbr")
  err <- tryCatch(gsdb_msigdb(collection = "ZZZ"), error = conditionMessage)
  skip_if(grepl("download|connect|internet|cache", err, ignore.case = TRUE),
          "MSigDB unavailable")
  expect_match(err, "msigdbr_collections|no gene sets|collection")
})

test_that("truncated ortholog mapping is an error, not a quiet result", {
  # The seatbelt: even if the cache-clearing workaround stops working, a
  # truncated collection must never be returned. Simulated by defeating the
  # workaround the same way the bug does -- poison the cache, then call the
  # coverage check directly on the truncated table.
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("babelgene")
  raw <- function(coll, sub = NULL) {
    a <- list(db_species = "HS", species = "Mus musculus", collection = coll)
    if (!is.null(sub)) a$subcollection <- sub
    tryCatch(do.call(msigdbr::msigdbr, a),
             error = function(e) skip(paste("MSigDB unavailable:",
                                            conditionMessage(e))))
  }
  .msigdbr_drop_ortholog_cache()
  invisible(raw("H"))                                    # poison with Hallmark
  args <- list(db_species = "HS", species = "Mus musculus",
               collection = "C2", subcollection = "CP:REACTOME")
  truncated <- raw("C2", "CP:REACTOME")
  skip_if(length(unique(truncated$db_gene_symbol)) > 8000,
          "msigdbr no longer leaks its ortholog cache across collections")

  expect_error(.msigdbr_assert_ortholog_coverage(truncated, args),
               "ortholog mapping dropped")
  # And the honest result passes the same check unchanged.
  clean <- gsdb_msigdb("Mus musculus", collection = "C2",
                       subcollection = "CP:REACTOME", db_species = "HS")
  expect_s3_class(clean, "gs_db")
})

test_that("the coverage check stands down when no ortholog mapping happens", {
  # Human target and mouse-native sets never hit the buggy code path, so the
  # check must not fire -- and must not cost a reference fetch either.
  args_hs <- list(db_species = "HS", species = "Homo sapiens", collection = "H")
  args_mm <- list(db_species = "MM", species = "Mus musculus", collection = "H")
  bad <- data.frame(gene_symbol = "A", db_gene_symbol = "A")
  expect_silent(.msigdbr_assert_ortholog_coverage(bad, args_hs))
  expect_silent(.msigdbr_assert_ortholog_coverage(bad, args_mm))
})
