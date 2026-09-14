fake_coresh_loading_object <- function(gse = "GSE_LOAD", n_genes = 30L) {
  E <- outer(
    seq_len(n_genes),
    seq_len(5L),
    function(i, j) sin(i * 0.37 + j * 0.61) + i / 100
  )
  list(
    gseId = gse,
    gplId = "GPL_LOAD",
    E1024 = round(E * 1024),
    rownames = seq_len(n_genes),
    totalVar = sum(E^2)
  )
}

test_that("CoReSh query provenance is compact and uses unique sizes", {
  queries <- list(
    iron_uptake = c(7037L, 4891L, 7037L, 55240L),
    ferroptosis = c(2879L, 3764L, 112398L, 79071L)
  )

  expect_identical(
    bulkiRNA:::.coresh_query_summary(queries),
    "iron_uptake(3), ferroptosis(4)"
  )
})

test_that("coresh_loadings reproduces the deterministic projection", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  obj <- fake_coresh_loading_object()
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(obj),
    .package = "bulkiRNA"
  )

  local_pinned_rng()
  query <- c(1L, 3L, 5L)
  set.seed(101)
  rng_before <- .Random.seed
  first <- coresh_loadings(chunk_path, "GSE_LOAD", query, top_n = 8L)
  rng_after <- .Random.seed
  set.seed(202)
  again <- coresh_loadings(chunk_path, "GSE_LOAD", query, top_n = 8L)

  E <- obj$E1024 / 1024
  query_idx <- match(unique(query), obj$rownames)
  profile <- colSums(E[query_idx, , drop = FALSE])
  expected <- as.numeric(E %*% (profile / sqrt(sum(profile^2))))
  ordering <- order(-abs(expected), method = "radix")
  keep <- head(ordering, 8L)

  expect_identical(first, again)
  expect_identical(rng_after, rng_before)
  expect_identical(names(first), c("gse", "gpl", "entrez", "loading", "rank"))
  expect_identical(unique(first$gpl), "GPL_LOAD")
  expect_identical(first$entrez, obj$rownames[keep])
  expect_equal(first$loading, expected[keep])
  expect_identical(first$rank, seq_len(8L))
})

test_that("coresh_loadings reports duplicate query Entrez IDs", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(fake_coresh_loading_object()),
    .package = "bulkiRNA"
  )

  expect_message(
    coresh_loadings(chunk_path, "GSE_LOAD", c(1L, 2L, 3L, 1L)),
    "Dropped 1 duplicate Entrez ID from `query`"
  )
})

test_that("coresh_loadings selects an explicit platform", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- coresh_micro_fixture()
  platform_a <- fixture[[1L]]
  platform_a$gseId <- "GSE_MULTI"
  platform_a$gplId <- "GPL_A"
  platform_z <- platform_a
  platform_z$gplId <- "GPL_Z"
  platform_z$E1024 <- platform_z$E1024[
    rev(seq_len(nrow(platform_z$E1024))), , drop = FALSE
  ]
  query <- head(unique(platform_a$rownames[!is.na(platform_a$rownames)]), 8L)
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(platform_a, platform_z),
    .package = "bulkiRNA"
  )

  out <- coresh_loadings(
    chunk_path, "GSE_MULTI", query, top_n = 20L, gpl = "GPL_Z"
  )
  expected_for <- function(obj) {
    E <- obj$E1024 / 1024
    query_idx <- match(query, obj$rownames)
    profile <- colSums(E[query_idx, , drop = FALSE])
    loadings <- as.numeric(E %*% (profile / sqrt(sum(profile^2))))
    keep <- head(order(-abs(loadings), method = "radix"), 20L)
    loadings[keep]
  }

  expect_identical(unique(out$gpl), "GPL_Z")
  expect_equal(out$loading, expected_for(platform_z))
  expect_false(isTRUE(all.equal(out$loading, expected_for(platform_a))))
})

test_that("coresh_loadings warns and chooses ambiguities independently of order", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- coresh_micro_fixture()
  platform_z <- fixture[[1L]]
  platform_z$gseId <- "GSE_MULTI"
  platform_z$gplId <- "GPL_Z"
  platform_a <- platform_z
  platform_a$gplId <- "GPL_A"
  platform_a$E1024 <- platform_a$E1024[
    rev(seq_len(nrow(platform_a$E1024))), , drop = FALSE
  ]
  query <- head(unique(platform_a$rownames[!is.na(platform_a$rownames)]), 8L)
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  object_order <- list(platform_z, platform_a)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) object_order,
    .package = "bulkiRNA"
  )

  expect_warning(
    first <- coresh_loadings(chunk_path, "GSE_MULTI", query, top_n = 20L),
    "GPL_A, GPL_Z.*Using GPL_A.*radix order"
  )
  object_order <- rev(object_order)
  expect_warning(
    reversed <- coresh_loadings(
      chunk_path, "GSE_MULTI", query, top_n = 20L
    ),
    "GPL_A, GPL_Z.*Using GPL_A.*radix order"
  )

  expect_identical(first, reversed)
  expect_identical(unique(first$gpl), "GPL_A")
})

test_that("coresh_loadings keeps a single platform quiet", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- coresh_micro_fixture()
  obj <- fixture[[1L]]
  query <- head(unique(obj$rownames[!is.na(obj$rownames)]), 8L)
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(obj),
    .package = "bulkiRNA"
  )

  expect_no_warning(
    out <- coresh_loadings(
      chunk_path, as.character(obj$gseId), query, top_n = 20L
    )
  )
  expect_identical(unique(out$gpl), as.character(obj$gplId))
})

test_that("coresh_loadings names available platforms when one is absent", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- coresh_micro_fixture()
  platform_z <- fixture[[1L]]
  platform_z$gseId <- "GSE_MULTI"
  platform_z$gplId <- "GPL_Z"
  platform_a <- platform_z
  platform_a$gplId <- "GPL_A"
  query <- head(unique(platform_a$rownames[!is.na(platform_a$rownames)]), 8L)
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(platform_z, platform_a),
    .package = "bulkiRNA"
  )

  expect_error(
    coresh_loadings(
      chunk_path, "GSE_MULTI", query, top_n = 20L, gpl = "GPL_MISSING"
    ),
    "GPL_MISSING.*choose one of: GPL_A, GPL_Z"
  )
})

test_that("coresh_loadings matches the projection on real edge-case objects", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- coresh_micro_fixture()
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) fixture,
    .package = "bulkiRNA"
  )

  for (fixture_name in c("duplicate_ids", "na_ids", "pca_reduced")) {
    obj <- fixture[[fixture_name]]
    query <- head(unique(obj$rownames[!is.na(obj$rownames)]), 8L)
    out <- coresh_loadings(
      chunk_path, as.character(obj$gseId), query, top_n = 20L
    )

    E <- obj$E1024 / 1024
    query_idx <- match(query, obj$rownames)
    profile <- colSums(E[query_idx, , drop = FALSE])
    expected <- as.numeric(E %*% (profile / sqrt(sum(profile^2))))
    keep <- head(order(-abs(expected), method = "radix"), 20L)

    expect_equal(out$loading, expected[keep], info = fixture_name)
    expect_identical(out$entrez, obj$rownames[keep], info = fixture_name)
  }
})

test_that("coresh_sets removes real NA Entrez IDs before symbol mapping", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  skip_if_not_installed("qs2")
  fixture <- coresh_micro_fixture()
  obj <- fixture$na_ids
  expect_true(anyNA(obj$rownames))

  chunk_dir <- withr::local_tempdir()
  chunk_path <- file.path(chunk_dir, "001_full_objects.qs2")
  qs2::qs_save(fixture["na_ids"], chunk_path)
  query <- head(unique(obj$rownames[!is.na(obj$rownames)]), 8L)
  top_hits <- tibble::tibble(
    query_name = "q",
    gse = as.character(obj$gseId),
    gpl = as.character(obj$gplId),
    rank = 1L
  )

  mapped_entrez <- NULL
  testthat::local_mocked_bindings(
    entrez_to_gene = function(entrez, species = "human") {
      mapped_entrez <<- entrez
      if (anyNA(entrez)) stop("symbol mapping received an NA Entrez ID")
      stats::setNames(rep("SHARED_SYMBOL", length(entrez)), entrez)
    },
    .package = "bulkiRNA"
  )

  db <- coresh_sets(
    top_hits,
    list(q = query),
    chunk_dir = chunk_dir,
    species = "human",
    top_n = nrow(obj$E1024),
    min_size = 1L,
    max_size = nrow(obj$E1024)
  )

  expect_false(anyNA(mapped_entrez))
  expect_lt(length(mapped_entrez), nrow(obj$E1024))
  expect_identical(unname(db[[1L]]), "SHARED_SYMBOL")
})

test_that("coresh_loadings validates coverage and controls", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunk_path <- withr::local_tempfile(fileext = "_full_objects.qs2")
  file.create(chunk_path)
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(fake_coresh_loading_object()),
    .package = "bulkiRNA"
  )

  expect_identical(formals(coresh_loadings)$top_n, 50L)
  expect_null(formals(coresh_loadings)$gpl)
  expect_error(coresh_loadings(chunk_path, "GSE_LOAD", c(1, 2, 3)),
               "integer vector")
  expect_error(coresh_loadings(chunk_path, "GSE_LOAD", 1:2),
               "at least 3")
  expect_error(coresh_loadings(chunk_path, "GSE_OTHER", 1:3),
               "was not found")
  expect_error(coresh_loadings(chunk_path, "GSE_LOAD", 1:3, top_n = 0),
               "`top_n`")
  expect_error(coresh_loadings(chunk_path, "GSE_LOAD", 1:3, gpl = ""),
               "`gpl`")
  expect_error(coresh_loadings("missing.qs2", "GSE_LOAD", 1:3),
               "does not exist")
})

test_that("coresh_sets keeps the higher-ranked overlapping hit", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  index <- tibble::tibble(
    gse = c("GSE_LOW", "GSE_HIGH", "GSE_FAR"),
    gpl = c("GPL1", "GPL1", "GPL1"),
    chunk = c("low.qs2", "high.qs2", "far.qs2")
  )
  attr(index, "provenance") <- list(
    source = "coresh",
    snapshot = "syn-test",
    path = "/coresh/hsa",
    species = "hsa",
    n_chunks = 3L
  )
  testthat::local_mocked_bindings(
    coresh_chunks = function(...) index,
    coresh_loadings = function(chunk_path, gse_id, query, top_n = 50L,
                               gpl = NULL) {
      ids <- switch(
        gse_id,
        GSE_HIGH = 1:10,
        GSE_LOW = c(1:9, 11L),
        GSE_FAR = 20:29
      )
      tibble::tibble(
        gse = gse_id,
        gpl = gpl,
        entrez = as.integer(ids),
        loading = rev(seq_along(ids)),
        rank = seq_along(ids)
      )
    },
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )
  top_hits <- tibble::tibble(
    query_name = c("q_low", "q_far", "q_high"),
    gse = c("GSE_LOW", "GSE_FAR", "GSE_HIGH"),
    gpl = "GPL1",
    rank = c(2L, 3L, 1L)
  )
  queries <- list(
    q_low = 1:3,
    q_far = 1:3,
    q_high = 1:3,
    unused_short_query = 1:2
  )

  db <- coresh_sets(
    top_hits,
    queries,
    min_size = 1L,
    max_size = 20L,
    jaccard_threshold = 0.8
  )

  expect_s3_class(db, "gs_db")
  expect_identical(
    names(db),
    c("CORESH_q_high_GSE_HIGH", "CORESH_q_far_GSE_FAR")
  )
  expect_false("CORESH_q_low_GSE_LOW" %in% names(db))
  expect_identical(attr(db, "provenance")$snapshot, "syn-test")
  expect_identical(
    attr(db, "provenance")$queries,
    "q_low(3), q_far(3), q_high(3), unused_short_query(2)"
  )
  expect_false("queries_r" %in% names(attr(db, "provenance")))
  expect_false("selected_hits_r" %in% names(attr(db, "provenance")))
  expect_identical(attr(db, "provenance")$top_n, 50L)
  expect_identical(
    attr(db, "set_provenance")$set_name,
    names(db)
  )
  expect_identical(
    attr(db, "set_provenance")$rank_in_coresh,
    c(1L, 3L)
  )
  expect_identical(attr(db, "set_provenance")$gpl, c("GPL1", "GPL1"))

  sub <- db["CORESH_q_far_GSE_FAR"]
  expect_identical(attr(sub, "provenance"), attr(db, "provenance"))
  expect_identical(
    attr(sub, "set_provenance")$set_name,
    "CORESH_q_far_GSE_FAR"
  )
})

test_that("coresh_sets resolves a bare accession by GPL and records it", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  index_order <- tibble::tibble(
    gse = c("GSE_MULTI", "GSE_MULTI"),
    gpl = c("GPL_Z", "GPL_A"),
    chunk = c("z.qs2", "a.qs2")
  )
  attr(index_order, "provenance") <- list(
    source = "coresh", snapshot = "syn-test", path = "/coresh/hsa",
    species = "hsa", n_chunks = 2L
  )
  testthat::local_mocked_bindings(
    coresh_chunks = function(...) index_order,
    coresh_loadings = function(chunk_path, gse_id, query, top_n = 50L,
                               gpl = NULL) {
      ids <- if (gpl == "GPL_A") 1:5 else 6:10
      tibble::tibble(
        gse = gse_id, gpl = gpl, entrez = as.integer(ids),
        loading = rev(seq_along(ids)), rank = seq_along(ids)
      )
    },
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )
  top_hits <- tibble::tibble(
    query_name = "q", gse = "GSE_MULTI", rank = 1L
  )

  expect_warning(
    first <- coresh_sets(top_hits, list(q = 1:3), min_size = 1L),
    "GPL_A, GPL_Z.*Using GPL_A.*radix order"
  )
  index_provenance <- attr(index_order, "provenance")
  index_order <- index_order[2:1, , drop = FALSE]
  attr(index_order, "provenance") <- index_provenance
  expect_warning(
    reversed <- coresh_sets(top_hits, list(q = 1:3), min_size = 1L),
    "GPL_A, GPL_Z.*Using GPL_A.*radix order"
  )

  expect_identical(first, reversed)
  expect_identical(unname(first[[1L]]), paste0("G", 1:5))
  expect_identical(attr(first, "set_provenance")$gpl, "GPL_A")
  expect_identical(attr(first, "set_provenance")$chunk_path, "a.qs2")
})

test_that("coresh_sets keeps distinct platforms for one accession", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  index <- tibble::tibble(
    gse = c("GSE_MULTI", "GSE_MULTI"),
    gpl = c("GPL_A", "GPL_Z"),
    chunk = c("a.qs2", "z.qs2")
  )
  attr(index, "provenance") <- list(
    source = "coresh", snapshot = "syn-test", path = "/coresh/hsa",
    species = "hsa", n_chunks = 2L
  )
  testthat::local_mocked_bindings(
    coresh_chunks = function(...) index,
    coresh_loadings = function(chunk_path, gse_id, query, top_n = 50L,
                               gpl = NULL) {
      ids <- if (gpl == "GPL_A") 1:5 else 6:10
      tibble::tibble(
        gse = gse_id, gpl = gpl, entrez = as.integer(ids),
        loading = rev(seq_along(ids)), rank = seq_along(ids)
      )
    },
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )
  top_hits <- tibble::tibble(
    query_name = c("q", "q"),
    gse = c("GSE_MULTI", "GSE_MULTI"),
    gpl = c("GPL_Z", "GPL_A"),
    rank = c(2L, 1L)
  )

  db <- coresh_sets(
    top_hits, list(q = 1:3), min_size = 1L, jaccard_threshold = 1
  )

  expect_identical(
    names(db),
    c("CORESH_q_GSE_MULTI_GPL_A", "CORESH_q_GSE_MULTI_GPL_Z")
  )
  expect_identical(attr(db, "set_provenance")$gpl, c("GPL_A", "GPL_Z"))
})

test_that("coresh_sets distinguishes lookup, collision, size, and empty outcomes", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  index <- tibble::tibble(
    gse = c("GSE_OK", "GSE_BROKEN"),
    gpl = c("GPL1", "GPL1"),
    chunk = c("ok.qs2", "broken.qs2")
  )
  attr(index, "provenance") <- list(
    source = "coresh", snapshot = "syn-test", path = "/coresh/hsa",
    species = "hsa", n_chunks = 1L
  )
  testthat::local_mocked_bindings(
    coresh_chunks = function(...) index,
    coresh_loadings = function(chunk_path, gse_id, query, top_n = 50L,
                               gpl = NULL) {
      if (gse_id == "GSE_BROKEN") stop("broken loading extraction")
      tibble::tibble(
        gse = gse_id, gpl = gpl, entrez = 1:3,
        loading = 3:1, rank = 1:3
      )
    },
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )
  queries <- list(q = 1:3)
  partial <- tibble::tibble(
    query_name = c("q", "q"),
    gse = c("GSE_OK", "GSE_MISSING"),
    rank = 1:2
  )

  lookup_messages <- testthat::capture_messages(
    db <- coresh_sets(partial, queries, min_size = 1L)
  )
  expect_true(any(grepl("chunk lookup failed for 1 of 2", lookup_messages)))
  expect_false(any(grepl("extraction", lookup_messages)))
  expect_identical(names(db), "CORESH_q_GSE_OK")
  expect_error(
    coresh_sets(partial[2L, ], queries, min_size = 1L),
    "failed for all 1 attempted hits"
  )

  extraction_hits <- tibble::tibble(
    query_name = c("q", "q"),
    gse = c("GSE_OK", "GSE_BROKEN"),
    rank = 1:2
  )
  extraction_messages <- testthat::capture_messages(
    extracted <- coresh_sets(extraction_hits, queries, min_size = 1L)
  )
  expect_true(any(grepl(
    "loading extraction or symbol mapping failed for 1 of 2",
    extraction_messages
  )))
  expect_false(any(grepl("chunk lookup failed", extraction_messages)))
  expect_identical(names(extracted), "CORESH_q_GSE_OK")

  collision_hits <- partial[c(1L, 1L), ]
  collision_hits$rank <- 1:2
  collision_messages <- testthat::capture_messages(
    collision <- coresh_sets(collision_hits, queries, min_size = 1L)
  )
  expect_true(any(grepl(
    "set-name collision skipped 1 of 2", collision_messages
  )))
  expect_false(any(grepl("extraction", collision_messages)))
  expect_identical(names(collision), "CORESH_q_GSE_OK")

  empty_messages <- testthat::capture_messages(
    empty <- coresh_sets(partial[1L, ], queries, min_size = 4L)
  )
  expect_true(any(grepl("size filter dropped 1 of 1", empty_messages)))
  expect_true(any(grepl(
    "1 hits attempted, 0 sets produced", empty_messages
  )))
  expect_s3_class(empty, "gs_db")
  expect_length(empty, 0L)
  expect_identical(nrow(attr(empty, "set_provenance")), 0L)

  zero_messages <- testthat::capture_messages(
    zero <- coresh_sets(partial[0L, ], queries, min_size = 1L)
  )
  expect_true(any(grepl(
    "0 hits attempted, 0 sets produced", zero_messages
  )))
  expect_length(zero, 0L)
})

test_that("coresh_sets validates its inputs and parameters", {
  hits <- tibble::tibble(query_name = "q", gse = "GSE1", rank = 1L)
  queries <- list(q = 1:3)

  expect_identical(formals(coresh_sets)$top_n, 50L)
  expect_identical(formals(coresh_sets)$min_size, 15L)
  expect_identical(formals(coresh_sets)$max_size, 500L)
  expect_identical(formals(coresh_sets)$jaccard_threshold, 0.8)
  expect_error(coresh_sets(hits[-1], queries), "missing column")
  expect_error(coresh_sets(hits, list(other = 1:3)), "absent from `queries`")
  expect_error(coresh_sets(hits, list(q = 1:2)), "at least 3 unique")
  expect_error(coresh_sets(hits, queries, min_size = 10L, max_size = 5L),
               "must not exceed")
  expect_error(coresh_sets(hits, queries, top_n = 5L, min_size = 15L),
               "`min_size`.*15.*`top_n`.*5.*only retain or reduce")
  expect_error(coresh_sets(hits, queries, jaccard_threshold = 1.1),
               "`jaccard_threshold`")
  expect_error(coresh_sets(hits, queries, verbose = NA), "`verbose`")
  expect_error(
    coresh_sets(dplyr::mutate(hits, gpl = NA_character_), queries),
    "`top_hits\\$gpl`"
  )
})
