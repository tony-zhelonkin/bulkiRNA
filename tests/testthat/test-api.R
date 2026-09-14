# These helpers normalise zero-length input so a comparison of two empty
# structures passes rather than failing on a names-or-type mismatch. That is not
# defensive padding. At v1.0.0 the 21 deprecated exports go and
# `bulkirna_api(lifecycle = "deprecated")` returns zero rows, so the shim
# allowlists below become permanently empty -- at which point the empty case
# stops being the untested path and becomes the only path, and this
# normalisation is the whole assertion. Do not delete these as dead code or
# inline them into a comparison that only works on non-empty input.
api_expect_set_allowlist <- function(observed, expected, info = NULL) {
  normalise <- function(x) sort(unique(as.character(x)))
  expect_identical(normalise(observed), normalise(expected), info = info)
}

api_expect_named_allowlist <- function(observed, expected, info = NULL) {
  normalise <- function(x) {
    if (!length(x)) {
      return(stats::setNames(character(0L), character(0L)))
    }
    x
  }
  expect_identical(normalise(observed), normalise(expected), info = info)
}

api_expect_machine_readable_targets <- function(api, deprecated, message_only) {
  machine_readable <- setdiff(
    as.character(deprecated),
    as.character(message_only)
  )
  checked <- character(0L)
  for (name in machine_readable) {
    expect_identical(
      api$superseded_by[api$name == name],
      bulkiRNA:::.bulkirna_deprecation_target(name)
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, machine_readable)
}

test_that("the API registry covers the complete namespace exactly once", {
  api <- bulkirna_api(quiet = TRUE)
  exports <- getNamespaceExports("bulkiRNA")

  expect_s3_class(api, "tbl_df")
  expect_identical(
    names(api),
    c("name", "layer", "lifecycle", "frozen", "stochastic", "superseded_by",
      "removed_in")
  )
  expect_equal(nrow(api), length(exports))
  # 64 at 1.1.0: + bulki_palettes, gs_palette, gs_plot_size, gs_scale_fonts
  # (stable) and coresh_labels (experimental).
  expect_equal(nrow(api), 64L)
  expect_equal(anyDuplicated(api$name), 0L)
  expect_setequal(api$name, exports)
  expect_identical(api$name, sort(api$name))
})

test_that("each export has one stability lifecycle and one layer", {
  api <- bulkirna_api(quiet = TRUE)

  expect_false(anyNA(api$lifecycle))
  expect_false(anyNA(api$layer))
  expect_true(all(nzchar(api$layer)))
  expect_true(all(api$lifecycle %in% c("stable", "experimental", "deprecated")))
  expect_equal(sum(api$lifecycle == "stable"), 51L)
  expect_equal(sum(api$lifecycle == "experimental"), 13L)
  expect_equal(sum(api$lifecycle == "deprecated"), 0L)
})

test_that("experimental status covers CoReSh, coregulation and gene-id helpers", {
  api <- bulkirna_api(quiet = TRUE)
  experimental <- c(
    # coresh_labels is new at 1.1.0 and stays experimental with the rest of the
    # CoReSh layer: a verb written this release has not earned a semver freeze.
    "coresh_chunks", "coresh_convergence", "coresh_labels", "coresh_loadings",
    "coresh_match",
    "coresh_search", "coresh_sets", "coresh_validate", "gs_coregulation",
    "gsdb_coresh", "entrez_to_gene", "filter_confounder_genes",
    "gene_to_entrez"
  )

  expect_setequal(api$name[api$lifecycle == "experimental"], experimental)
  expect_true(all(api$layer[grepl("^coresh_", api$name)] == "coresh"))
  expect_identical(api$layer[api$name == "gs_coregulation"], "gs")
  expect_identical(api$layer[api$name == "gsdb_coresh"], "gsdb")
})

test_that("the remaining historical signature freeze is recorded", {
  api <- bulkirna_api(quiet = TRUE)
  frozen <- c("build_dge", "ensure_dir", "format_pathway_name")

  expect_equal(sum(api$frozen), length(frozen))
  expect_setequal(api$name[api$frozen], frozen)
  expect_true(all(api$lifecycle[api$frozen] == "stable"))
  api_expect_set_allowlist(
    api$name[api$frozen & api$lifecycle == "stable"],
    c("build_dge", "ensure_dir", "format_pathway_name")
  )
  api_expect_set_allowlist(
    character(0L),
    NULL,
    info = "The frozen-and-stable allowlist must be shape-stable when empty."
  )
})

test_that("deprecated metadata columns remain empty at v1.0.0", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api[api$lifecycle == "deprecated", ]
  maintained <- api[api$lifecycle != "deprecated", ]

  expect_equal(nrow(deprecated), 0L)
  expect_true(all(is.na(api$superseded_by)))
  expect_true(all(is.na(api$removed_in)))
  expect_true(all(!is.na(deprecated$superseded_by)))
  expect_true(all(nzchar(deprecated$superseded_by)))
  expect_true(all(deprecated$removed_in == "1.0.0"))
  expect_true(all(is.na(maintained$superseded_by)))
  expect_true(all(is.na(maintained$removed_in)))

  no_deprecated <- deprecated[0L, , drop = FALSE]
  expect_true(all(!is.na(no_deprecated$superseded_by)))
  expect_true(all(nzchar(no_deprecated$superseded_by)))
  expect_true(all(no_deprecated$removed_in == "1.0.0"))
})

test_that("machine-readable successor checks accept the empty tier", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]
  api_expect_machine_readable_targets(api, deprecated, character(0L))
  api_expect_machine_readable_targets(
    api[0L, , drop = FALSE],
    character(0L),
    character(0L)
  )
})

test_that("message-only successor checks accept the empty tier", {
  api_expect_named_allowlist(
    character(0L),
    NULL
  )
  api_expect_named_allowlist(
    character(0L),
    NULL,
    info = "The message-only shim allowlist must be shape-stable when empty."
  )
})

test_that("bulkirna_api is itself a stable top-level export", {
  row <- bulkirna_api(quiet = TRUE)
  row <- row[row$name == "bulkirna_api", ]

  expect_equal(nrow(row), 1L)
  expect_identical(row$lifecycle, "stable")
  expect_identical(row$layer, "top-level")
  expect_false(row$frozen)
  expect_true(is.na(row$superseded_by))
  expect_true(is.na(row$removed_in))
})

test_that("G4 widens a stable vocabulary with an experimental verb", {
  api <- bulkirna_api(quiet = TRUE)
  stat_types <- api[api$name == "gs_stat_types", ]
  coregulation <- api[api$name == "gs_coregulation", ]

  expect_identical(stat_types$lifecycle, "stable")
  expect_false(stat_types$frozen)
  expect_identical(coregulation$lifecycle, "experimental")
  expect_identical(coregulation$layer, "gs")
  expect_false(coregulation$frozen)
})

api_expect_deprecated_first <- function(deprecated) {
  checked <- character(0L)
  for (name in deprecated) {
    # Namespace lookup, not getExportedValue(): since 1.0.0 these are internal
    # fixtures, which getExportedValue() cannot see.
    fun <- get(name, envir = asNamespace("bulkiRNA"), inherits = FALSE)
    first <- body(fun)[[2L]]
    expect_true(
      is.call(first) && identical(first[[1L]], as.name(".Deprecated")),
      info = name
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, deprecated)
}

api_expect_deprecated_warnings <- function(deprecated) {
  checked <- character(0L)
  for (name in deprecated) {
    # Namespace lookup, not getExportedValue(): since 1.0.0 these are internal
    # fixtures, which getExportedValue() cannot see.
    fun <- get(name, envir = asNamespace("bulkiRNA"), inherits = FALSE)
    # The test above establishes that this is the shim's first statement.
    # Install a body-truncated copy under the real shim name. This preserves
    # .Deprecated()'s call-derived `old` value without allowing delegated work
    # in this or a future shim to run.
    deprecation_call <- body(fun)[[2L]]
    warning_fun <- fun
    body(warning_fun) <- as.call(list(as.name("{"), deprecation_call))
    shim_env <- new.env(parent = environment(fun))
    assign(name, warning_fun, envir = shim_env)
    # Both shim forms name themselves, in two spellings: `.Deprecated("target")`
    # produces `'name' is deprecated` from the call, and the six message-only
    # shims write `` `name()` is deprecated `` in their own prose. Matching the
    # bare name covers both, and still fails if the truncated body loses the
    # call-derived `old` value and the warning reads `'.Deprecated' is
    # deprecated` instead.
    expect_warning(
      eval(call(name), envir = shim_env),
      name,
      fixed = TRUE,
      class = "deprecatedWarning"
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, deprecated)
}

# The 21 names demoted to internal fixtures at 1.0.0. A closed historical set:
# it can never grow, so a literal is the anchor. The golden baseline at 752481f
# runs through these, so their deprecation contract still has to hold.
API_LEGACY_FIXTURES <- c(
  "convert_human_to_mouse", "create_MD_plot", "create_standard_volcano",
  "custom_minimal_theme_with_grid", "download_gatom_references",
  "empty_gsea_tibble", "filter_by_size", "gsea_barplot", "gsea_dotplot",
  "gsea_dotplot_facet", "gsea_running_sum_plot", "list_reference_dbs",
  "list_to_term2gene", "load_reference_db", "normalize_gsea_results",
  "parse_gmx", "parse_mitoxplorer", "plot_all_gsea_results", "run_gsea",
  "run_gsea_analysis", "save_gsea_log"
)

test_that("the legacy fixtures are internal and exactly the expected set", {
  expect_length(API_LEGACY_FIXTURES, 21L)
  exports <- getNamespaceExports("bulkiRNA")
  ns <- asNamespace("bulkiRNA")

  # Present in the namespace, absent from the public surface. Fails if a
  # demotion is reverted or a fixture is deleted outright.
  for (name in API_LEGACY_FIXTURES) {
    expect_true(exists(name, envir = ns, inherits = FALSE), info = name)
    expect_false(name %in% exports, info = name)
  }
  expect_false(any(API_LEGACY_FIXTURES %in% bulkirna_api(quiet = TRUE)$name))
})

test_that("every legacy fixture calls .Deprecated first", {
  api_expect_deprecated_first(API_LEGACY_FIXTURES)
  api_expect_deprecated_first(character(0L))
})

test_that("every legacy fixture warns before doing any work", {
  api_expect_deprecated_warnings(API_LEGACY_FIXTURES)
  api_expect_deprecated_warnings(character(0L))
})

test_that("a name's layer is invariant to lifecycle selection", {
  api <- bulkirna_api(quiet = TRUE)
  for (lifecycle in c("stable", "experimental", "deprecated")) {
    selected <- bulkirna_api(lifecycle, quiet = TRUE)
    expect_identical(
      selected$layer,
      api$layer[match(selected$name, api$name)]
    )
  }
  selected <- api[0L, , drop = FALSE]
  expect_identical(
    selected$layer,
    api$layer[match(selected$name, api$name)]
  )
})

test_that("lifecycle selection follows the dependency-report interface", {
  api <- bulkirna_api(quiet = TRUE)
  stable <- bulkirna_api("stable", quiet = TRUE)
  selected <- bulkirna_api(c("experimental", "deprecated"), quiet = TRUE)

  expect_true(all(stable$lifecycle == "stable"))
  api_expect_set_allowlist(
    unique(selected$lifecycle),
    intersect(c("experimental", "deprecated"), unique(api$lifecycle))
  )
  api_expect_set_allowlist(
    character(0L),
    NULL,
    info = "Lifecycle selection must be shape-stable when a tier is empty."
  )
  expect_error(
    bulkirna_api("not-a-lifecycle", quiet = TRUE),
    "should be one of"
  )
  expect_error(
    bulkirna_api(quiet = NA),
    "`quiet` must be a single non-missing logical value"
  )
  expect_invisible(bulkirna_api("experimental", quiet = FALSE))
})
