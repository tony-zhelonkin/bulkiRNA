fake_gsdb_coresh_objects <- function() {
  lapply(seq_len(4L), function(i) {
    expression <- outer(
      seq_len(8L) + i,
      c(3, -2, 5, 1),
      function(gene, component) gene * component
    )
    list(
      gseId = paste0("GSE_PROVIDER_", i),
      gplId = paste0("GPL_PROVIDER_", i),
      E1024 = expression * 1024,
      rownames = seq_len(8L),
      totalVar = sum(expression^2)
    )
  })
}

local_gsdb_coresh_chunks <- function(objects, mapper = NULL,
                                     .local_envir = parent.frame()) {
  chunk_dir <- withr::local_tempdir(.local_envir = .local_envir)
  chunk_path <- file.path(chunk_dir, "chunk_1_full_objects.qs2")
  file.create(chunk_path)
  if (is.null(mapper)) {
    mapper <- function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    }
  }
  list(
    chunk_dir = chunk_dir,
    chunk_path = chunk_path,
    objects = objects,
    mapper = mapper
  )
}

test_that("provider matches hand composition and subsets provenance", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  objects <- fake_gsdb_coresh_objects()
  fixture <- local_gsdb_coresh_chunks(objects)
  testthat::local_mocked_bindings(
    .coresh_chunk_files = function(chunk_dir) fixture$chunk_path,
    .coresh_read_chunk = function(path) fixture$objects,
    entrez_to_gene = fixture$mapper,
    .package = "bulkiRNA"
  )
  chunk_dir <- fixture$chunk_dir
  queries <- list(alpha = 1:3, beta = 4:6)

  hits <- coresh_search(
    queries,
    chunk_dir = chunk_dir,
    species = "human",
    n_cores = 1L,
    seed = 17L
  )
  manual <- coresh_sets(
    hits[hits$rank <= 2L, , drop = FALSE],
    queries,
    chunk_dir = chunk_dir,
    species = "human",
    top_n = 5L,
    min_size = 1L,
    max_size = 5L,
    jaccard_threshold = 1
  )

  set.seed(408L)
  stream_before <- .Random.seed
  provider <- gsdb_coresh(
    queries,
    chunk_dir = chunk_dir,
    species = "human",
    top_hits = 2L,
    top_n = 5L,
    min_size = 1L,
    max_size = 5L,
    jaccard_threshold = 1,
    n_cores = 1L,
    seed = 17L
  )

  expect_identical(provider, manual)
  expect_identical(.Random.seed, stream_before)
  expect_s3_class(provider, "gs_db")
  expect_identical(attr(provider, "database"), "coresh")
  expect_identical(attr(provider, "species"), "Homo sapiens")

  provenance <- attr(provider, "provenance")
  expect_identical(provenance$queries, "alpha(3), beta(3)")
  expect_false("queries_r" %in% names(provenance))
  expect_false("selected_hits_r" %in% names(provenance))
  expect_identical(provenance$top_hits, 2L)
  expect_identical(provenance$top_n, 5L)
  expect_identical(provenance$min_size, 1L)
  expect_identical(provenance$max_size, 5L)
  expect_identical(provenance$jaccard_threshold, 1)

  # Gate 4: database provenance is invariant, while set provenance follows
  # both the membership and ordering of the retained sets.
  subset <- provider[c(4L, 2L)]
  expect_identical(attr(subset, "database"), "coresh")
  expect_identical(attr(subset, "provenance"), provenance)
  expect_identical(
    attr(subset, "set_provenance"),
    attr(provider, "set_provenance")[c(4L, 2L), ]
  )
  expect_identical(attr(subset, "set_provenance")$set_name, names(subset))

  withr::local_options(width = 72L)
  printed <- capture.output(print(provider))
  provenance_lines <- grep("^(Provenance: | {12})", printed, value = TRUE)
  expect_gt(length(provenance_lines), 1L)
  expect_true(any(grepl("queries=alpha(3)", printed, fixed = TRUE)))
  expect_true(any(grepl("beta(3)", printed, fixed = TRUE)))
  expect_false(any(grepl("structure(list", printed, fixed = TRUE)))
})

test_that("the CoReSh provider works with compute, render, and persist layers", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  fixture <- local_gsdb_coresh_chunks(fake_gsdb_coresh_objects())
  testthat::local_mocked_bindings(
    .coresh_chunk_files = function(chunk_dir) fixture$chunk_path,
    .coresh_read_chunk = function(path) fixture$objects,
    entrez_to_gene = fixture$mapper,
    .package = "bulkiRNA"
  )
  chunk_dir <- fixture$chunk_dir
  db <- gsdb_coresh(
    list(alpha = 1:3, beta = 4:6),
    chunk_dir = chunk_dir,
    top_hits = 2L,
    top_n = 5L,
    min_size = 1L,
    max_size = 5L,
    jaccard_threshold = 1,
    n_cores = 1L
  )
  universe <- c(
    unique(unlist(db, use.names = FALSE)),
    paste0("BACKGROUND_", seq_len(20L))
  )
  res <- gs_test(
    utils::head(db[[1L]], 2L),
    db,
    universe = universe,
    min_size = 1L,
    max_size = 5L
  )

  expect_s3_class(res, "gs_result")
  expect_gt(nrow(res), 0L)
  expect_true(all(res$database == "coresh"))
  expect_s3_class(gs_plot_dot(res), "ggplot")
  expect_s3_class(gs_plot_bar(res), "ggplot")

  output_dir <- withr::local_tempdir()
  written <- gs_write(res, output_dir, name = "coresh_provider")
  expect_identical(as.character(written), output_dir)
  expect_true(length(attr(written, "files")) > 0L)
  expect_true(all(file.exists(attr(written, "files"))))
})

test_that("gsdb_coresh records and warns about an unidentified current", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  root <- withr::local_tempdir()
  chunk_dir <- file.path(
    root, "coresh", "current", "preprocessed_chunks", "hsa"
  )
  dir.create(chunk_dir, recursive = TRUE)
  chunk_path <- file.path(chunk_dir, "chunk_1_full_objects.qs2")
  file.create(chunk_path)
  withr::local_envvar(REFCACHE_ROOT = root)
  .clear_ref_resolutions()
  on.exit(.clear_ref_resolutions(), add = TRUE)
  testthat::local_mocked_bindings(
    .coresh_chunk_files = function(chunk_dir) chunk_path,
    .coresh_read_chunk = function(path) fake_gsdb_coresh_objects(),
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )

  expect_warning(
    db <- gsdb_coresh(
      list(alpha = 1:3),
      species = "human",
      top_hits = 1L,
      top_n = 5L,
      min_size = 1L,
      max_size = 5L,
      n_cores = 1L
    ),
    "is not a symlink"
  )
  expect_identical(attr(db, "provenance")$snapshot, "current")
})

test_that("gsdb_coresh preserves honest empty and failed outcomes", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  collapsed_mapper <- function(entrez, species = "human") {
    stats::setNames(rep("ONE_GENE", length(entrez)), entrez)
  }
  fixture <- local_gsdb_coresh_chunks(
    fake_gsdb_coresh_objects(),
    mapper = collapsed_mapper
  )
  testthat::local_mocked_bindings(
    .coresh_chunk_files = function(chunk_dir) fixture$chunk_path,
    .coresh_read_chunk = function(path) fixture$objects,
    entrez_to_gene = fixture$mapper,
    .package = "bulkiRNA"
  )
  chunk_dir <- fixture$chunk_dir

  messages <- testthat::capture_messages(
    empty <- gsdb_coresh(
      list(alpha = 1:3),
      chunk_dir = chunk_dir,
      top_hits = 1L,
      top_n = 3L,
      min_size = 2L,
      max_size = 3L,
      n_cores = 1L
    )
  )
  expect_s3_class(empty, "gs_db")
  expect_length(empty, 0L)
  expect_true(any(grepl("size filter dropped 1 of 1", messages)))
  expect_true(any(grepl("1 hits attempted, 0 sets produced", messages)))

  expect_error(
    gsdb_coresh(
      list(absent = 101:103),
      chunk_dir = chunk_dir,
      top_hits = 1L,
      top_n = 3L,
      min_size = 1L,
      max_size = 3L,
      n_cores = 1L
    ),
    "failed for all 1 attempted hits"
  )
})

test_that("gsdb_coresh validates provider controls and keeps planned defaults", {
  expect_null(formals(gsdb_coresh)$chunk_dir)
  expect_identical(formals(gsdb_coresh)$species, "human")
  expect_identical(formals(gsdb_coresh)$top_hits, 5L)
  expect_identical(formals(gsdb_coresh)$top_n, 50L)
  expect_identical(formals(gsdb_coresh)$min_size, 15L)
  expect_identical(formals(gsdb_coresh)$max_size, 500L)
  expect_identical(formals(gsdb_coresh)$jaccard_threshold, 0.8)
  expect_identical(formals(gsdb_coresh)$n_cores, 4L)
  expect_identical(formals(gsdb_coresh)$seed, 1L)
  expect_error(gsdb_coresh(list(q = 1:3), top_hits = 0L), "`top_hits`")
  expect_error(gsdb_coresh(list(q = 1:3), seed = Inf), "`seed`")
})
