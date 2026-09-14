fake_coresh_object <- function(gse = "GSE1", total_var = 10) {
  list(
    gseId = gse,
    gplId = "GPL1",
    E1024 = matrix(
      c(1024L, 0L, 9999L, 9999L, 0L, 2048L),
      nrow = 3L,
      byrow = TRUE
    ),
    rownames = c(10L, 10L, 20L),
    totalVar = total_var
  )
}

test_that("coresh_match computes pct_var from stored totalVar", {
  obj <- fake_coresh_object()
  out <- coresh_match(obj, c(10L, 20L, 99L), pvalues = FALSE)

  # match() selects rows 1 and 3; the duplicate Entrez row 2 is deliberately
  # not selected. The scaled query profile is therefore c(1, 2).
  expected <- sum(c(1, 2)^2) / 2 / 10 * 100
  expect_identical(names(out),
                   c("gse", "gpl", "pct_var", "p_value", "log2err", "size"))
  expect_equal(out$pct_var, expected)
  expect_identical(out$p_value, NA_real_)
  expect_identical(out$log2err, NA_real_)
  expect_identical(out$size, 2L)
})

test_that("coresh_match returns the zero row for an absent query", {
  out <- coresh_match(fake_coresh_object(), 999L)

  expect_equal(out$pct_var, 0)
  expect_identical(out$p_value, NA_real_)
  expect_identical(out$log2err, NA_real_)
  expect_identical(out$size, 0L)
})

fake_coresh_signal_object <- function(gse = "GSE_SIGNAL", strength = 0.5) {
  n_genes <- 80L
  n_samples <- 8L
  query_size <- 8L
  E <- outer(
    seq_len(n_genes),
    seq_len(n_samples),
    function(i, j) 0.25 * sin(i * 1.7 + j * 0.9) +
      0.25 * cos(i * j * 0.31)
  )
  latent <- c(-1.5, -1, -0.5, -0.25, 0.25, 0.5, 1, 1.5)
  E[seq_len(query_size), ] <- E[seq_len(query_size), ] + matrix(
    rep(strength * latent, each = query_size),
    nrow = query_size
  )
  E1024 <- round(E * 1024)
  ids <- seq.int(1001L, length.out = n_genes)
  ids[[n_genes]] <- ids[[1L]]
  list(
    gseId = gse,
    gplId = "GPL_SIGNAL",
    E1024 = E1024,
    rownames = ids,
    totalVar = sum(E^2)
  )
}

test_that("coresh_match computes reproducible GESECA p-values", {
  obj <- fake_coresh_signal_object()
  query <- obj$rownames[seq_len(8L)]

  without <- coresh_match(obj, query, pvalues = FALSE)
  first <- coresh_match(obj, query, pvalues = TRUE, seed = 17L)
  again <- coresh_match(obj, query, pvalues = TRUE, seed = 17L)
  different <- coresh_match(obj, query, pvalues = TRUE, seed = 29L)

  expect_identical(names(first), names(without))
  expect_equal(first$pct_var, without$pct_var)
  expect_true(is.finite(first$p_value))
  expect_gt(first$p_value, 0)
  expect_lte(first$p_value, 1)
  # `Inf` is a legitimate value here, meaning the estimate is past reliable
  # resolution, so assert only that the estimator reported a bound at all.
  expect_false(is.na(first$log2err))
  expect_identical(first$p_value, again$p_value)
  expect_false(isTRUE(all.equal(first$p_value, different$p_value)))
})

test_that("coresh_match handles missing and duplicate Entrez row names", {
  obj <- fake_coresh_signal_object()
  obj$rownames[[10L]] <- NA_integer_
  query <- obj$rownames[seq_len(8L)]

  without <- coresh_match(obj, query, pvalues = FALSE)
  with <- coresh_match(obj, query, pvalues = TRUE, seed = 17L)

  expect_equal(with$pct_var, without$pct_var)
  expect_true(is.finite(with$p_value))
  expect_gt(with$pct_var, 0)
  expect_identical(with$size, 8L)
})

test_that("coresh_match handles an almost entirely unmapped background", {
  obj <- fake_coresh_signal_object()
  query_rows <- seq_len(2L)
  background_rows <- rep(3:80, length.out = 9998L)
  obj$E1024 <- obj$E1024[c(query_rows, background_rows), , drop = FALSE]
  obj$rownames <- rep(NA_integer_, 10000L)
  obj$rownames[query_rows] <- seq.int(1001L, length.out = 2L)
  obj$totalVar <- sum((obj$E1024 / 1024)^2)
  query <- obj$rownames[query_rows]

  # Modelled on a real dataset with 9,998 unmapped rows out of 10,000.
  without <- coresh_match(obj, query, pvalues = FALSE)
  with <- coresh_match(obj, query, pvalues = TRUE, seed = 17L)

  expect_equal(with$pct_var, without$pct_var)
  expect_true(is.finite(with$p_value))
  expect_gt(with$pct_var, 0)
  expect_identical(with$size, 2L)
})

test_that("the real CoReSh micro fixture retains its edge cases", {
  fixture <- coresh_micro_fixture()

  expect_identical(
    names(fixture),
    c("plain", "duplicate_ids", "na_ids", "pca_reduced")
  )
  expect_gt(anyDuplicated(fixture$duplicate_ids$rownames), 0L)
  expect_true(anyNA(fixture$na_ids$rownames))
  expect_lt(
    ncol(fixture$pca_reduced$E1024),
    fixture$pca_reduced$nsamples
  )
})

test_that("coresh_match handles every real CoReSh micro fixture object", {
  fixture <- coresh_micro_fixture()
  result_columns <- c(
    "gse", "gpl", "pct_var", "p_value", "log2err", "size"
  )

  # The pca_reduced object's columns are components, not samples; it must use
  # the same scoring path as every other object.
  for (fixture_name in names(fixture)) {
    obj <- fixture[[fixture_name]]
    reference_ids <- unique(obj$rownames[!is.na(obj$rownames)])
    if (fixture_name == "na_ids") {
      # Missing query IDs are rejected, while missing background IDs remain
      # supported by the valid call below.
      expect_error(
        coresh_match(obj, c(NA_integer_, head(reference_ids, 7L))),
        "non-empty integer vector",
        info = fixture_name
      )
      query <- head(reference_ids, 8L)
    } else if (fixture_name == "duplicate_ids") {
      duplicated_id <- obj$rownames[
        duplicated(obj$rownames) & !is.na(obj$rownames)
      ][[1L]]
      query <- c(
        duplicated_id,
        head(setdiff(reference_ids, duplicated_id), 7L)
      )
      expect_true(
        sum(obj$rownames == duplicated_id, na.rm = TRUE) > 1L,
        info = fixture_name
      )
    } else {
      query <- head(reference_ids, 8L)
    }

    expected_size <- length(intersect(unique(query), reference_ids))
    without <- coresh_match(obj, query, pvalues = FALSE)
    with <- coresh_match(obj, query, pvalues = TRUE)

    expect_equal(nrow(without), 1L, info = fixture_name)
    expect_equal(nrow(with), 1L, info = fixture_name)
    expect_identical(names(without), result_columns, info = fixture_name)
    expect_identical(names(with), result_columns, info = fixture_name)
    expect_true(is.finite(without$pct_var), info = fixture_name)
    expect_true(without$pct_var >= 0, info = fixture_name)
    expect_identical(with$pct_var, without$pct_var, info = fixture_name)
    expect_identical(without$p_value, NA_real_, info = fixture_name)
    expect_true(is.finite(with$p_value), info = fixture_name)
    expect_true(with$p_value > 0, info = fixture_name)
    expect_true(with$p_value <= 1, info = fixture_name)
    expect_identical(
      without$size,
      as.integer(expected_size),
      info = fixture_name
    )
    expect_identical(with$size, without$size, info = fixture_name)
    if (fixture_name == "duplicate_ids") {
      # The duplicated reference row contributes one matched Entrez ID.
      expect_identical(
        without$size,
        as.integer(length(unique(query))),
        info = fixture_name
      )
    }
  }
})

test_that("coresh_match sample_size controls estimator precision", {
  # Keep the signal moderate so the smaller estimator does not saturate at
  # infinite uncertainty before the precision comparison can be made.
  obj <- fake_coresh_signal_object(strength = 0.15)
  query <- obj$rownames[seq_len(8L)]

  small <- coresh_match(
    obj, query, pvalues = TRUE, sample_size = 21L, seed = 17L
  )
  large <- coresh_match(
    obj, query, pvalues = TRUE, sample_size = 101L, seed = 17L
  )

  expect_true(is.finite(small$log2err))
  expect_lt(large$log2err, small$log2err)
})

test_that("coresh_match reports a dataset when GESECA cannot test its query", {
  obj <- fake_coresh_signal_object(gse = "GSE_ALL_ROWS")
  obj$gplId <- "GPL_ALL_ROWS"
  obj$E1024 <- obj$E1024[seq_len(8L), , drop = FALSE]
  obj$rownames <- obj$rownames[seq_len(8L)]
  obj$totalVar <- sum((obj$E1024 / 1024)^2)

  error <- tryCatch(
    coresh_match(obj, obj$rownames, pvalues = TRUE, seed = 17L),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "GSE_ALL_ROWS", fixed = TRUE)
  expect_match(conditionMessage(error), "GPL_ALL_ROWS", fixed = TRUE)
  expect_match(conditionMessage(error), "k = 8", fixed = TRUE)
  expect_match(conditionMessage(error), "nrow(E) = 8", fixed = TRUE)
  expect_false(grepl("subscript out of bounds", conditionMessage(error),
                     fixed = TRUE))
})

test_that("coresh_match de-duplicates query Entrez IDs before scoring", {
  obj <- fake_coresh_signal_object()
  deduplicated <- obj$rownames[seq_len(8L)]
  query <- c(deduplicated, deduplicated[[1L]])
  messages <- character()

  out <- withCallingHandlers(
    coresh_match(obj, query, pvalues = TRUE, seed = 17L),
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  expected <- coresh_match(obj, deduplicated, pvalues = TRUE, seed = 17L)

  expect_length(messages, 1L)
  expect_match(messages, "Dropped 1 duplicate Entrez ID", fixed = TRUE)
  expect_identical(out$size, 8L)
  expect_equal(out$pct_var, expected$pct_var)
  expect_identical(out$p_value, expected$p_value)
})

test_that("coresh_match leaves an absent query untested with pvalues", {
  out <- coresh_match(
    fake_coresh_signal_object(), 999999L, pvalues = TRUE
  )

  expect_equal(out$pct_var, 0)
  expect_identical(out$p_value, NA_real_)
  expect_identical(out$log2err, NA_real_)
  expect_identical(out$size, 0L)
})

test_that("coresh_match validates its supported arguments", {
  obj <- fake_coresh_object()
  expect_identical(formals(coresh_match)$sample_size, 21L)
  expect_identical(formals(coresh_match)$seed, 1L)
  expect_identical(formals(coresh_match)$eps, 1e-300)
  expect_error(coresh_match(obj, c(10, 20)), "integer vector")
  expect_error(coresh_match(obj, integer()), "non-empty")
  expect_error(coresh_match(obj, 10L, pvalues = NA), "`pvalues`")
  expect_error(coresh_match(obj, 10L, sample_size = 0), "`sample_size`")
  expect_error(coresh_match(obj, 10L, seed = Inf), "`seed`")
  expect_error(coresh_match(obj, 10L, eps = 0), "`eps`")

  bad <- obj
  bad$totalVar <- NULL
  expect_error(coresh_match(bad, 10L), "totalVar")
  bad <- obj
  bad$rownames <- as.character(bad$rownames)
  expect_error(coresh_match(bad, 10L), "integer Entrez")
})

test_that("coresh_chunks errors clearly for an empty directory", {
  empty <- withr::local_tempdir()
  expect_error(
    coresh_chunks(empty, cache = FALSE),
    "references/synapse-data-setup.md",
    fixed = TRUE
  )
  expect_error(coresh_chunks(empty, species = "rat"), "`species`")
  expect_error(coresh_chunks(empty, cache = NA), "`cache`")
  expect_error(coresh_chunks(1), "`chunk_dir`")
})

test_that("coresh_chunks maps species aliases and caches by resolved path", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  root <- withr::local_tempdir()
  hsa <- file.path(root, "coresh", "current", "preprocessed_chunks", "hsa")
  mmu <- file.path(root, "coresh", "current", "preprocessed_chunks", "mmu")
  dir.create(hsa, recursive = TRUE)
  dir.create(mmu, recursive = TRUE)
  file.create(file.path(hsa, "001_full_objects.qs2"))
  file.create(file.path(mmu, "001_full_objects.qs2"))
  withr::local_envvar(c(REFCACHE_ROOT = root))

  reads <- 0L
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      reads <<- reads + 1L
      gse <- if (basename(dirname(path)) == "hsa") "GSE_HSA" else "GSE_MMU"
      list(fake_coresh_object(gse))
    },
    .package = "bulkiRNA"
  )

  human <- coresh_chunks(species = "human", cache = FALSE)
  expect_identical(human$gse, "GSE_HSA")
  expect_identical(attr(human, "provenance")$species, "hsa")
  expect_identical(coresh_chunks(species = "hsa")$gse, "GSE_HSA")
  expect_identical(coresh_chunks(species = "Homo sapiens")$gse, "GSE_HSA")
  expect_identical(reads, 1L)

  mouse <- coresh_chunks(species = "mouse", cache = FALSE)
  expect_identical(mouse$gse, "GSE_MMU")
  expect_identical(attr(mouse, "provenance")$species, "mmu")
  expect_identical(coresh_chunks(species = "mmu")$gse, "GSE_MMU")
  expect_identical(reads, 2L)
})

test_that("coresh_search ranks pct_var descending without optional engines", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  file.create(file.path(chunks, "002_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      if (grepl("001_", basename(path), fixed = TRUE)) {
        list(fake_coresh_object("GSE_HIGH", total_var = 10))
      } else {
        list(fake_coresh_object("GSE_LOW", total_var = 20))
      }
    },
    .require_pkg = function(...) {
      stop("an optional engine was requested", call. = FALSE)
    },
    .package = "bulkiRNA"
  )

  out <- coresh_search(
    list(query_a = c(10L, 20L, 99L)),
    chunk_dir = chunks,
    n_cores = 1L,
    pvalues = FALSE
  )
  expect_identical(names(out), c(
    "query_name", "gse", "gpl", "pct_var", "p_value", "log2err", "size",
    "rank"
  ))
  expect_identical(out$gse, c("GSE_HIGH", "GSE_LOW"))
  expect_identical(out$rank, 1:2)
  expect_true(all(is.na(out$p_value)))
  expect_true(all(is.na(out$log2err)))
  expect_identical(attr(out, "provenance")$species, "hsa")
})

test_that("coresh_search keeps platforms as separate dataset rows", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      platform_z <- fake_coresh_object("GSE_MULTI", total_var = 10)
      platform_z$gplId <- "GPL_Z"
      platform_a <- fake_coresh_object("GSE_MULTI", total_var = 20)
      platform_a$gplId <- "GPL_A"
      list(platform_z, platform_a)
    },
    .package = "bulkiRNA"
  )

  out <- coresh_search(
    list(query_a = c(10L, 20L, 99L)), chunks,
    n_cores = 1L, pvalues = FALSE
  )

  expect_identical(nrow(out), 2L)
  expect_identical(out$gse, c("GSE_MULTI", "GSE_MULTI"))
  expect_identical(out$gpl, c("GPL_Z", "GPL_A"))
  expect_false(identical(out$pct_var[[1L]], out$pct_var[[2L]]))
})

test_that("coresh_search ranks p-values and reproduces a fixed seed", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  file.create(file.path(chunks, "002_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      if (grepl("001_", basename(path), fixed = TRUE)) {
        medium <- fake_coresh_signal_object("GSE_MEDIUM", strength = 0.35)
        high <- fake_coresh_signal_object("GSE_HIGH", strength = 0.65)
        # Stored totalVar changes pct_var but is not an input to GESECA, so
        # these factors force the variance and p-value rankings to disagree.
        high$totalVar <- high$totalVar * 1e6
        list(medium, high)
      } else {
        low <- fake_coresh_signal_object("GSE_LOW", strength = 0.2)
        low$totalVar <- low$totalVar / 1e6
        list(low)
      }
    },
    .package = "bulkiRNA"
  )
  query <- seq.int(1001L, length.out = 8L)

  first <- coresh_search(
    list(signal = query), chunks, n_cores = 1L, pvalues = TRUE, seed = 17L
  )
  again <- coresh_search(
    list(signal = query), chunks, n_cores = 1L, pvalues = TRUE, seed = 17L
  )
  variance <- coresh_search(
    list(signal = query), chunks, n_cores = 1L, pvalues = FALSE, seed = 17L
  )

  expect_identical(first$p_value, again$p_value)
  expect_identical(names(first), names(variance))
  expect_identical(first$rank, seq_len(nrow(first)))
  expect_identical(variance$rank, seq_len(nrow(variance)))
  expect_equal(first$p_value, sort(first$p_value))
  expect_equal(variance$pct_var, sort(variance$pct_var, decreasing = TRUE))
  expect_false(identical(first$gse, variance$gse))
})

test_that("coresh_search reports duplicate query IDs only once", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  file.create(file.path(chunks, "002_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      list(fake_coresh_signal_object(basename(path)))
    },
    .package = "bulkiRNA"
  )
  query <- c(seq.int(1001L, length.out = 8L), 1001L)
  messages <- character()

  out <- withCallingHandlers(
    coresh_search(
      list(signal = query), chunks, n_cores = 1L, pvalues = FALSE
    ),
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )

  expect_length(messages, 1L)
  expect_match(messages, "Dropped 1 duplicate Entrez ID", fixed = TRUE)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$size == 8L))
})

test_that("coresh_match restores the caller's RNG state", {
  obj <- fake_coresh_signal_object()
  query <- obj$rownames[seq_len(8L)]
  local_pinned_rng()

  set.seed(90210L)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  control <- runif(1L)
  assign(".Random.seed", seed_before, envir = .GlobalEnv)

  coresh_match(obj, query, pvalues = TRUE, seed = 17L)

  expect_identical(RNGkind(), kind_before)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
  expect_identical(runif(1L), control)
})

test_that("coresh_match restores an absent RNG seed", {
  obj <- fake_coresh_signal_object()
  query <- obj$rownames[seq_len(8L)]
  local_pinned_rng()
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  coresh_match(obj, query, pvalues = TRUE, seed = 17L)

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("coresh_match quietly restores a legacy sample kind", {
  obj <- fake_coresh_signal_object()
  query <- obj$rownames[seq_len(8L)]
  local_pinned_rng()
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  suppressWarnings(set.seed(90210L))
  kind_before <- RNGkind()

  expect_no_warning(coresh_match(obj, query, pvalues = TRUE, seed = 17L))
  expect_identical(RNGkind(), kind_before)
})

test_that("coresh_search p-values do not depend on parallel scheduling", {
  skip_on_os("windows")
  skip_if_not_installed("BiocParallel")
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  file.create(file.path(chunks, "002_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      strength <- if (grepl("001_", basename(path), fixed = TRUE)) 0.35 else 0.6
      list(fake_coresh_signal_object(basename(path), strength = strength))
    },
    .package = "bulkiRNA"
  )
  query <- seq.int(1001L, length.out = 8L)

  serial <- coresh_search(
    list(signal = query), chunks, n_cores = 1L, pvalues = TRUE, seed = 23L
  )
  parallel <- coresh_search(
    list(signal = query), chunks, n_cores = 2L, pvalues = TRUE, seed = 23L
  )

  serial <- serial[order(serial$gse), ]
  parallel <- parallel[order(parallel$gse), ]
  expect_gt(nrow(serial), 0L)
  expect_true(all(is.finite(serial$p_value)))
  expect_identical(serial$p_value, parallel$p_value)
})

test_that("coresh_match pins its RNG inside a BiocParallel task", {
  skip_if_not_installed("BiocParallel")
  obj <- fake_coresh_signal_object()
  query <- seq.int(1001L, length.out = 8L)

  parent <- coresh_match(obj, query, pvalues = TRUE, seed = 23L)
  nested <- BiocParallel::bplapply(
    1L,
    function(i) coresh_match(obj, query, pvalues = TRUE, seed = 23L),
    BPPARAM = BiocParallel::SerialParam()
  )[[1L]]

  expect_true(is.finite(parent$p_value))
  expect_identical(parent$p_value, nested$p_value)
})

test_that("coresh_search names invalid queries and validates controls", {
  chunks <- withr::local_tempdir()
  expect_identical(formals(coresh_search)$sample_size, 21L)
  expect_identical(formals(coresh_search)$seed, 1L)
  expect_identical(formals(coresh_search)$eps, 1e-300)
  expect_error(
    coresh_search(list(short_query = c(1L, 2L)), chunks, n_cores = 1L),
    "short_query"
  )
  expect_error(coresh_search(list(c(1L, 2L, 3L)), chunks), "named list")
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, n_cores = 0),
    "`n_cores`"
  )
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, sample_size = 0),
    "`sample_size`"
  )
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, seed = Inf),
    "`seed`"
  )
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, eps = 0),
    "`eps`"
  )
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, species = "rat"),
    "`species`"
  )
})

test_that("coresh_convergence summarizes independent top-query support", {
  ranking <- tibble::tibble(
    query_name = c("q1", "q2", "q3", "q1", "q2", "q3", "q1", "q2"),
    gse = c(
      "GSE_A", "GSE_A", "GSE_A", "GSE_B", "GSE_B", "GSE_C", "GSE_A",
      "GSE_A"
    ),
    gpl = c(
      "GPL1", "GPL1", "GPL1", "GPL1", "GPL1", "GPL1", "GPL2", "GPL2"
    ),
    pct_var = c(10, 8, 6, 9, 7, 12, 1, 2),
    rank = c(1L, 2L, 12L, 2L, 4L, 1L, 3L, 5L)
  )

  out <- coresh_convergence(ranking, top_n = 10L, min_queries = 2L)
  expect_identical(out$gse, c("GSE_A", "GSE_B", "GSE_A"))
  expect_identical(out$gpl, c("GPL1", "GPL1", "GPL2"))
  expect_identical(out$n_queries, c(2L, 2L, 2L))
  expect_identical(out$queries, rep("q1, q2", 3L))
  expect_identical(out$best_rank, c(1L, 2L, 3L))
  expect_equal(out$mean_pct_var, c(9, 8, 1.5))
})

test_that("coresh_convergence validates its inputs", {
  expect_error(coresh_convergence(list()), "data frame")
  expect_error(coresh_convergence(data.frame(gse = "GSE1")), "missing column")
  ranking <- tibble::tibble(
    query_name = "q", gse = "GSE1", gpl = "GPL1",
    pct_var = 1, rank = 1L
  )
  expect_error(coresh_convergence(ranking[-3L]), "`gpl`")
  expect_error(coresh_convergence(ranking, top_n = 0), "`top_n`")
  expect_error(coresh_convergence(ranking, min_queries = NA), "`min_queries`")
})

test_that("coresh_validate reports every preflight failure without stopping", {
  empty <- withr::local_tempdir()
  visible <- NULL
  output <- capture.output(
    visible <- withVisible(coresh_validate(empty, species = "hsa"))
  )
  out <- visible$value

  expect_s3_class(out, "tbl_df")
  expect_false(visible$visible)
  expect_identical(names(out), c("check", "ok", "detail"))
  expect_true(length(output) > 0L)
  expect_true(all(paste0("package: ", c(
    "qs2", "coresh", "BiocParallel", "org.Hs.eg.db", "org.Mm.eg.db"
  )) %in% out$check))
  coresh_row <- out[out$check == "package: coresh", ]
  expect_identical(coresh_row$ok, requireNamespace("coresh", quietly = TRUE))
  expect_match(coresh_row$detail, "no R code or callable functions", fixed = TRUE)
  expect_match(coresh_row$detail, "Not required", fixed = TRUE)
  expect_match(coresh_row$detail, "fgsea::geseca()", fixed = TRUE)
  expect_false(out$ok[out$check == "chunk files"])
  expect_match(
    out$detail[out$check == "chunk files"],
    "references/synapse-data-setup.md",
    fixed = TRUE
  )
  expect_error(coresh_validate(empty, species = "rat"), "`species`")
  expect_error(coresh_validate(1), "`chunk_dir`")
})

test_that("a real CoReSh chunk index is readable when the refcache is mounted", {
  resolved <- tryCatch(
    .ref_path("coresh", "preprocessed_chunks", "hsa"),
    error = function(e) NULL
  )
  skip_if_not(!is.null(resolved) && dir.exists(resolved),
              "CoReSh refcache is not mounted")
  skip_if_not_installed("qs2")

  paths <- sort(list.files(
    resolved,
    pattern = "_full_objects\\.qs2$",
    full.names = TRUE
  ))
  skip_if_not(length(paths) > 0L, "CoReSh refcache has no human chunks")
  one_chunk <- withr::local_tempdir()
  linked <- suppressWarnings(file.symlink(
    paths[[1L]],
    file.path(one_chunk, basename(paths[[1L]]))
  ))
  skip_if_not(isTRUE(linked), "This filesystem cannot link a test chunk")

  # Exercise one real ~120 MB chunk, not the entire ~20 GB species tree.
  out <- coresh_chunks(chunk_dir = one_chunk, species = "human", cache = FALSE)
  expect_gt(nrow(out), 0L)
  expect_identical(names(out), c("gse", "gpl", "chunk"))
  expect_true(all(file.exists(unique(out$chunk))))
})

# --- Defect found by the G5 cross-check, not by this suite ------------------
# Two of ~42,500 mmu datasets (GSE63, GSE1457) carry totalVar = NA. Validation
# stop()ed on them inside a bplapply(), so one bad row killed the whole sweep:
# "1 remote error, 61 unevaluated". These tests pin the skip contract.

test_that("coresh_search skips an invalid dataset, and names it", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))

  broken <- fake_coresh_object("GSE_BROKEN")
  broken$totalVar <- NA_real_

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      list(fake_coresh_object("GSE_OK", total_var = 10), broken)
    },
    .require_pkg = function(...) stop("an optional engine was requested",
                                     call. = FALSE),
    .package = "bulkiRNA"
  )

  # Collect messages by hand: expect_message() returns the condition, not the
  # value, and the returned object is what carries the skip report.
  msgs <- character()
  out <- withCallingHandlers(
    coresh_search(
      list(q = c(10L, 20L, 99L)), chunk_dir = chunks, species = "human",
      n_cores = 1L, pvalues = FALSE
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("skipped 1 dataset", msgs, fixed = TRUE)))
  expect_true(any(grepl("GSE_BROKEN", msgs, fixed = TRUE)))

  # The valid dataset is still scored: the run completes rather than aborting.
  expect_identical(out$gse, "GSE_OK")

  skipped <- attr(out, "skipped")
  expect_s3_class(skipped, "data.frame")
  expect_identical(nrow(skipped), 1L)
  expect_identical(skipped$gse, "GSE_BROKEN")
  expect_identical(skipped$element, 2L)
  expect_identical(skipped$chunk, "001_full_objects.qs2")
  expect_match(skipped$reason, "totalVar")
})

test_that("coresh_search reports a zero-row skip frame for a clean sweep", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(fake_coresh_object("GSE_OK")),
    .require_pkg = function(...) stop("an optional engine was requested",
                                      call. = FALSE),
    .package = "bulkiRNA"
  )

  out <- coresh_search(
    list(q = c(10L, 20L, 99L)), chunk_dir = chunks, species = "human",
    n_cores = 1L, pvalues = FALSE
  )
  expect_identical(nrow(attr(out, "skipped")), 0L)
  expect_identical(
    names(attr(out, "skipped")),
    c("chunk", "element", "gse", "reason")
  )
})

test_that("a chunk where nothing validates is still an error", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))

  broken <- fake_coresh_object("GSE_BROKEN")
  broken$totalVar <- NA_real_

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(broken, broken),
    .require_pkg = function(...) stop("an optional engine was requested",
                                      call. = FALSE),
    .package = "bulkiRNA"
  )

  # An all-invalid chunk is indistinguishable from a broken snapshot, so it
  # must not degrade into a valid empty result.
  expect_error(
    coresh_search(
      list(q = c(10L, 20L, 99L)), chunk_dir = chunks, species = "human",
      n_cores = 1L, pvalues = FALSE
    ),
    "No dataset in the chunk at `path` passed validation"
  )
})
