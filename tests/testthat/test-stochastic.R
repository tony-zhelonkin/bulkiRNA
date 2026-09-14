test_that("every exported seed formal is declared stochastic", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  exports <- getNamespaceExports("bulkiRNA")
  has_seed <- vapply(exports, function(name) {
    "seed" %in% names(formals(getExportedValue("bulkiRNA", name)))
  }, logical(1L))

  expect_s3_class(registry, "tbl_df")
  expect_identical(
    names(registry),
    c("name", "seed_arg", "seed_default", "source_of_randomness", "note")
  )
  expect_setequal(exports[has_seed], registry$name[!is.na(registry$seed_arg)])
})

test_that("declared seed arguments and defaults match their formals", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  declared <- registry[!is.na(registry$seed_arg), ]

  for (i in seq_len(nrow(declared))) {
    fun <- getExportedValue("bulkiRNA", declared$name[[i]])
    formal_names <- names(formals(fun))
    seed_arg <- declared$seed_arg[[i]]
    seed_default <- paste(deparse(formals(fun)[[seed_arg]]), collapse = "")

    expect_true(seed_arg %in% formal_names, info = declared$name[[i]])
    expect_identical(
      declared$seed_default[[i]], seed_default,
      info = declared$name[[i]]
    )
  }
})

test_that("gs_test exposes its seed only through dots and uses it", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  through_dots <- registry[is.na(registry$seed_arg), ]

  expect_identical(through_dots$name, "gs_test")
  expect_true("..." %in% names(formals(gs_test)))
  expect_false("seed" %in% names(formals(gs_test)))
  adapter_default <- paste(
    deparse(formals(bulkiRNA:::.gs_fgsea)$seed), collapse = ""
  )
  expect_identical(through_dots$seed_default, adapter_default)

  first <- gs_test(
    fake_ranks(), fake_gs_db(), min_size = 5L, max_size = 50L,
    n_perm_simple = 1000L, seed = 101L
  )
  second <- gs_test(
    fake_ranks(), fake_gs_db(), min_size = 5L, max_size = 50L,
    n_perm_simple = 1000L, seed = 202L
  )

  expect_false(identical(first, second))
})

test_that("the API stochastic flag is derived from the registry", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  api <- bulkirna_api(quiet = TRUE)

  expect_setequal(api$name[api$stochastic], registry$name)
  expect_false(anyNA(api$stochastic))
})

test_that("only rng.R mutates RNG state", {
  r_dir <- NULL
  for (path in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(path) && length(list.files(path, "[.]R$"))) {
      r_dir <- path
      break
    }
  }
  if (is.null(r_dir)) {
    skip(paste(
      "Source-tree-only RNG ownership audit: package R/ sources are",
      "unavailable when tests run against the installed package."
    ))
  }

  call_name <- function(call) {
    head <- call[[1L]]
    if (is.symbol(head)) return(as.character(head))
    if (is.call(head) && length(head) == 3L &&
        as.character(head[[1L]]) %in% c("::", ":::")) {
      return(as.character(head[[3L]]))
    }
    NA_character_
  }
  literal_random_seed <- function(call, name) {
    args <- as.list(call)[-1L]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
    if (name == "assign") {
      target <- if ("x" %in% arg_names) {
        args[[which(arg_names == "x")[[1L]]]]
      } else if (length(args)) {
        args[[1L]]
      } else {
        NULL
      }
      return(is.character(target) && identical(target, ".Random.seed"))
    }
    candidates <- args[arg_names %in% c("", "list")]
    any(vapply(candidates, function(x) {
      is.character(x) && any(x == ".Random.seed")
    }, logical(1L)))
  }
  rng_mutation <- function(call) {
    name <- call_name(call)
    if (is.na(name)) return(NA_character_)
    if (name == "set.seed") return("set.seed()")
    if (name == "RNGkind" && length(call) > 1L) return("RNGkind(...)")
    if (name %in% c("assign", "rm") && literal_random_seed(call, name)) {
      return(paste0(name, "(\".Random.seed\", ...)"))
    }
    NA_character_
  }
  collect_mutations <- function(node) {
    found <- character()
    walk <- function(x) {
      if (is.call(x)) {
        mutation <- rng_mutation(x)
        if (!is.na(mutation)) found <<- c(found, mutation)
      }
      if (is.recursive(x)) lapply(as.list(x), walk)
      invisible(NULL)
    }
    walk(node)
    found
  }

  files <- list.files(r_dir, "[.]R$", full.names = TRUE)
  mutations <- unlist(lapply(files, function(path) {
    calls <- collect_mutations(parse(path, keep.source = TRUE))
    stats::setNames(calls, rep(basename(path), length(calls)))
  }))
  owner_mutations <- unname(mutations[names(mutations) == "rng.R"])
  # paste0() recycles a zero-length vector to "" when another argument is
  # longer, so building this unconditionally yields ": " on the passing case.
  outside <- which(names(mutations) != "rng.R")
  offenders <- if (length(outside)) {
    paste0(names(mutations)[outside], ": ", unname(mutations)[outside])
  } else {
    character(0L)
  }

  expect_true(
    "set.seed()" %in% owner_mutations,
    info = "R/rng.R must contain the package's set.seed() mutation"
  )
  expect_identical(
    offenders,
    character(0L),
    info = paste(
      "RNGkind() with no arguments is a permitted read;",
      "an RNG-state mutation belongs in .with_pinned_seed() in R/rng.R"
    )
  )
})

test_that("session provenance records RNG structure and stochastic names", {
  .clear_ref_resolutions()
  path <- tempfile(fileext = ".txt")
  write_session_provenance(path)
  lines <- readLines(path)
  registry <- bulkirna_stochastic(quiet = TRUE)

  rng_line <- grep("^RNGkind\\(\\):", lines, value = TRUE)
  expect_length(rng_line, 1L)
  expect_true(all(c("kind =", "normal.kind =", "sample.kind =") %in%
                    regmatches(
                      rng_line,
                      gregexpr("(?:normal\\.|sample\\.)?kind =", rng_line,
                               perl = TRUE)
                    )[[1L]]))
  for (name in registry$name) {
    expect_true(any(grepl(paste0("^  ", name, " = "), lines)))
  }
})

test_that("every declared stochastic function is classified and reproducible", {
  declared <- bulkirna_stochastic(quiet = TRUE)$name

  # Functions this test can call with nothing but the package's Imports.
  exercised <- c("gs_test", "coresh_match", "gs_coregulation")
  # Functions it cannot, each with the reason and where they are covered.
  covered_elsewhere <- c(
    gatom_module  = "needs gatom and mwcsr; covered in test-gatom.R",
    coresh_search = "needs a chunk tree; covered by the mocked tests above",
    gsdb_coresh = paste(
      "needs a chunk tree; composition and RNG state are covered in",
      "test-gsdb-coresh.R"
    )
  )

  # The point of this assertion: declaring another stochastic function fails
  # here until somebody decides which list it belongs in. Coverage cannot
  # silently lag the registry.
  expect_setequal(c(exercised, names(covered_elsewhere)), declared)
  expect_length(intersect(exercised, names(covered_elsewhere)), 0L)
  expect_true(all(nzchar(covered_elsewhere)))

  # State is captured inside, immediately around the calls, so that building a
  # fixture beforehand -- which itself seeds -- is not mistaken for the function
  # under test disturbing the caller's stream.
  same_seed_agrees <- function(call_with_seed, label) {
    kind_before <- RNGkind()
    stream_before <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }

    a <- call_with_seed(11L)
    b <- call_with_seed(11L)
    c2 <- call_with_seed(12L)

    expect_identical(a, b, info = label)
    expect_false(isTRUE(all.equal(a, c2)), info = label)

    expect_identical(RNGkind(), kind_before, info = label)
    stream_after <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    expect_identical(stream_after, stream_before, info = label)
  }

  for (fn in exercised) {
    if (fn == "gs_test") {
      # Deliberately mid-range p-values over many sets: the multilevel
      # estimator escalates only in the tail, so a tiny fixture with saturated
      # p-values would compare equal across seeds and prove nothing.
      set.seed(4)
      ranks <- sort(stats::setNames(stats::rnorm(2000), paste0("g", 1:2000)),
                    decreasing = TRUE)
      db <- fake_gs_db(sets = stats::setNames(
        lapply(1:16, function(i) sample(names(ranks), 40L)),
        paste0("SET", 1:16)
      ))
      same_seed_agrees(function(s) {
        gs_test(ranks, db, method = "fgsea", seed = s)$p_value
      }, fn)
    } else if (fn == "coresh_match") {
      set.seed(4)
      n <- 300L
      m <- matrix(stats::rnorm(n * 12L), nrow = n)
      latent <- stats::rnorm(12L)
      m[1:10, ] <- m[1:10, ] + 0.35 * matrix(latent, nrow = 10L, ncol = 12L,
                                             byrow = TRUE)
      m <- m - rowMeans(m)
      obj <- list(
        gseId = "GSE_LOOP", gplId = "GPL_LOOP",
        E1024 = round(m * 1024), rownames = seq_len(n),
        samples = paste0("s", 1:12), nsamples = 12L,
        totalVar = sum((round(m * 1024) / 1024)^2)
      )
      same_seed_agrees(function(s) {
        coresh_match(obj, 1:10, pvalues = TRUE, seed = s)$p_value
      }, fn)
    } else {
      fixture <- fake_coregulation_input()
      same_seed_agrees(function(s) {
        gs_coregulation(
          fixture$expr, fixture$db,
          min_size = 10L, max_size = 50L,
          sample_size = 21L, seed = s
        )$p_value
      }, fn)
    }

  }
})
