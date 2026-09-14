name_audit_r_dir <- function() {
  for (path in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(path) && length(list.files(path, "[.]R$"))) return(path)
  }
  NULL
}

name_audit_exports <- function() {
  r_dir <- name_audit_r_dir()
  if (is.null(r_dir)) return(NULL)
  namespace <- file.path(dirname(r_dir), "NAMESPACE")
  if (!file.exists(namespace)) return(NULL)
  lines <- readLines(namespace, warn = FALSE)
  sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", lines, value = TRUE))
}

name_audit_method_formals <- function(exports) {
  namespace <- asNamespace("bulkiRNA")
  method_table <- get(
    ".__S3MethodsTable__.",
    envir = namespace,
    inherits = FALSE
  )
  method_names <- ls(envir = method_table, all.names = TRUE)

  unique(unlist(lapply(method_names, function(method_name) {
    method <- get(method_name, envir = method_table, inherits = FALSE)
    exported_generics <- exports[vapply(exports, function(generic) {
      prefix <- paste0(generic, ".")
      if (!startsWith(method_name, prefix)) return(FALSE)
      class <- substring(method_name, nchar(prefix) + 1L)
      registered <- getS3method(
        generic,
        class,
        optional = TRUE,
        envir = namespace
      )
      identical(registered, method)
    }, logical(1L))]
    if (!length(exported_generics)) return(character(0L))

    # A complete exported-name prefix, verified against the registered method,
    # avoids truncating the generic at the method name's first dot. Prefer the
    # longest verified prefix if one generic name prefixes another.
    generic <- exported_generics[[which.max(nchar(exported_generics))]]
    class <- substring(method_name, nchar(generic) + 2L)
    names(formals(getS3method(generic, class, envir = namespace)))
  }), use.names = FALSE))
}

name_audit_expect_allowlist <- function(observed, allowlist, info = NULL) {
  expected <- names(allowlist)
  if (is.null(expected)) expected <- character(0L)
  expect_identical(sort(as.character(observed)), sort(expected), info = info)
}

test_that("every live export has a layer prefix or a reasoned exception", {
  exports <- name_audit_exports()
  if (is.null(exports)) {
    skip(paste(
      "Source-tree-only name audit: R/ and NAMESPACE are unavailable",
      "when tests run against the installed package."
    ))
  }
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])

  top_level <- c(
    # Original I/O verbs describe the action more clearly than an io_ prefix.
    read_counts_matrix = "specific file-reading verb",
    read_metadata = "specific file-reading verb",
    ensure_dir = "frozen utility used directly by two unmigrated consumers",
    write_session_provenance = "specific provenance-writing verb",
    # These construct the DE object or its annotation before the de_ renderers.
    build_dge = "frozen constructor predating the de_ renderer layer",
    annotate_genes = "annotation verb used before differential expression",
    # Identifier conversion and filtering are cross-layer gene operations.
    gene_to_entrez = "identifier conversion is not owned by one layer",
    entrez_to_gene = "identifier conversion is not owned by one layer",
    filter_confounder_genes = "cross-layer gene filtering verb",
    # Presentation names follow R/ggplot vocabulary and remain intentionally short.
    format_pathway_name = "general label formatter",
    theme_bulki = "house theme parallel to ggplot2 theme names",
    bulki_palettes = "palette table parallel to the house theme, not one layer's",
    # The bulkirna_ prefix marks package-wide metadata and dependency checks.
    bulkirna_api = "package metadata",
    bulkirna_check_deps = "package-wide dependency inspection",
    bulkirna_stochastic = "package metadata"
  )
  prefixed <- grepl("^(gsdb_|gs_|de_|gatom_|coresh_)", live)
  observed <- sort(live[!prefixed])

  name_audit_expect_allowlist(observed, top_level)
  name_audit_expect_allowlist(
    character(0L),
    NULL,
    info = "The top-level exception check must be shape-stable when empty."
  )
})

test_that("live export formals use the complete audited vocabulary", {
  exports <- name_audit_exports()
  if (is.null(exports)) {
    skip(paste(
      "Source-tree-only argument audit: R/ and NAMESPACE are unavailable",
      "when tests run against the installed package."
    ))
  }
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])

  concept <- function(formals, reason) {
    list(formals = formals, reason = reason)
  }
  single_use <- function(formals, reason) {
    stats::setNames(
      lapply(formals, function(formal) concept(formal, reason)),
      formals
    )
  }

  argument_concepts <- c(
    list(
      result_limit = concept(
        "top_n",
        paste(
          "This is the canonical spelling for the number of results to",
          "retain; see 04_NAME_AUDIT.md."
        )
      ),
      colour_palette = concept(
        "palette",
        paste(
          "This is the canonical spelling for a colour palette; see",
          "04_NAME_AUDIT.md."
        )
      ),
      b_statistic_cutoff = concept(
        "B_cutoff",
        paste(
          "limma names this statistic `B`, and the argument reads as that",
          "statistic's cutoff. Not frozen; kept to match upstream."
        )
      ),
      mean_expression = concept(
        "baseMean",
        "GATOM's DE input contract uses the upstream DESeq2 column spelling."
      ),
      log2_fold_change = concept(
        "log2FC",
        "GATOM's DE input contract uses the upstream GATOM column spelling."
      ),
      label_overlap_limit = concept(
        "max.overlaps",
        "The plotting API forwards ggrepel's upstream formal unchanged."
      ),
      sample_group = concept(
        "group",
        paste(
          "This renderer input names the sample annotation used to facet a",
          "score heatmap."
        )
      ),
      sample_selection = concept(
        "samples",
        paste(
          "This renderer input selects and orders the samples displayed in a",
          "score heatmap."
        )
      ),
      optional_schema_columns = concept(
        "optional",
        paste(
          "This schema accessor input chooses whether the optional columns are",
          "included, matching whether gs_to_master() received an entity_type."
        )
      )
    ),
    single_use(
      c(
        "species", "db", "contrast", "seed", "quiet", "verbose", "path",
        "min_size", "max_size", "n_cores", "dir"
      ),
      paste(
        "This is the sole canonical spelling of a required cross-layer",
        "concept."
      )
    ),
    single_use(
      c(
        "...", "biomart_host", "biomart_version", "by", "by_contrast",
        "by_direction", "cache", "center", "chunk_dir", "chunk_path", "coef",
        "collapse", "collection", "compare", "count_mat", "data", "database",
        "database_label", "db_species", "de", "de_results", "df", "dge",
        "direction", "download", "drop", "drop_empty", "ens_ids",
        "ensembl_version",
        "entity_type", "entrez", "eps", "error", "expr", "fc_cutoff",
        "fdr_cutoff", "features", "fit", "gene2reaction_extra", "genes",
        "genes_df", "genome_build", "gpl", "gse_id", "gsea_param", "id",
        "input_gene_name", "jaccard_threshold", "k_gene", "k_met", "kcdf",
        "lifecycle", "m", "max_genes", "met_de", "method", "metric",
        "alias_fallback",
        "ids", "titles", "type",
        "min_genes", "min_queries", "multi_vals", "n", "name", "network",
        "networks", "norm_method", "obj", "overwrite", "p_cutoff", "p_value",
        "padj", "pathway_id", "pathway_names", "pathways", "pattern", "per",
        "prefix", "pval", "pvalues", "queries", "query", "ranking", "ranks",
        "rebuild", "refs", "required_cols", "res", "round_nonint",
        "sample_col_candidates", "sample_data", "sample_size", "samples_df",
        "schema_version", "sets", "solver", "stat", "subcollection", "symbols",
        "table", "unique_genes", "universe", "use_biomart", "x"
      ),
      paste(
        "This API-specific input has one spelling and does not compete with",
        "another audited concept."
      )
    ),
    single_use(
      c(
        "aes_x", "annotate_counts", "base_family", "base_font_size",
        "base_size", "base_theme",
        "caption", "colour_by", "database_labels", "decision_by", "dpi",
        "es_ylim", "linewidth", "rescale_fonts", "rug_ylabels", "xticks",
        "style", "title_width",
        "facet", "fixed_p_boundary", "formats", "grid", "height", "highlight",
        "highlight_gene", "keep_first_caption", "label", "label_method",
        "label_size", "labels", "legend_pos", "legend_position", "limits",
        "max_name_length", "metric_label", "orientation", "overview",
        "padj_max", "panel_heights", "plot", "plots", "point_size", "prune",
        "scale", "shape_by", "show_grid", "show_quadrant_counts", "size_range",
        "sort_by", "stat_as_nes", "strip_prefix", "subtitle", "symbol_by",
        "text", "title",
        "top_hits", "use_formatting", "width", "wrap_width", "x_breaks",
        "xlim_abs", "y_padding", "ylim_abs"
      ),
      paste(
        "This renderer or presentation input has one spelling and controls",
        "one distinct aesthetic or display choice."
      )
    )
  )
  # Empty on purpose, and it is the goal state rather than an unfinished list:
  # every concept now has exactly one spelling. A new entry here is a decision
  # to accept a second spelling for one concept, and needs a reason saying why.
  multiple_spelling_exceptions <- list()
  non_snake_case_exceptions <- c(
    B_cutoff = "established limma B-statistic spelling",
    baseMean = "upstream DESeq2 column spelling",
    log2FC = "upstream GATOM column spelling",
    max.overlaps = "upstream ggrepel formal"
  )
  required <- c(
    "species", "db", "contrast", "seed", "quiet", "verbose", "path",
    "min_size", "max_size", "n_cores", "dir"
  )
  export_formals <- unlist(lapply(live, function(name) {
    names(formals(getExportedValue("bulkiRNA", name)))
  }), use.names = FALSE)
  method_formals <- name_audit_method_formals(live)
  observed <- sort(unique(c(export_formals, method_formals)))
  assignments <- unlist(lapply(argument_concepts, `[[`, "formals"),
                        use.names = FALSE)
  audited <- sort(unique(assignments))
  unexpected <- setdiff(observed, audited)
  unused <- setdiff(audited, observed)
  unused_concepts <- names(argument_concepts)[vapply(
    argument_concepts,
    function(x) any(x$formals %in% unused),
    logical(1L)
  )]
  # DO NOT DELETE `named_list()` AS DEAD CODE. It is the only reason the
  # assertion below passes today.
  #
  # Both sides must be a *named* list even when empty. Subsetting an empty
  # named list by `sort(names(.))` yields one; `lapply()` over an empty
  # unnamed list does not, so the two compared unequal on the passing case --
  # which only appeared once the last multi-spelling exception went away.
  #
  # Since the argument unification there is one concept per formal, all
  # singletons, so both sides are permanently empty and this normalisation *is*
  # the assertion rather than a defensive edge case. It stops being load-bearing
  # only if somebody reintroduces a second spelling for one concept -- which is
  # the thing the assertion exists to prevent. So an empty-versus-empty
  # comparison here is the success state, not a leftover.
  named_list <- function(x) {
    if (!length(x)) return(stats::setNames(list(), character(0L)))
    x[sort(names(x))]
  }
  actual_multiple <- named_list(lapply(
    argument_concepts[lengths(lapply(argument_concepts, `[[`, "formals")) > 1L],
    function(x) sort(x$formals)
  ))
  expected_multiple <- named_list(lapply(
    multiple_spelling_exceptions,
    function(x) sort(x$formals)
  ))
  multiple_names <- union(names(actual_multiple), names(expected_multiple))
  multiple_drift <- multiple_names[!vapply(multiple_names, function(name) {
    identical(actual_multiple[[name]], expected_multiple[[name]])
  }, logical(1L))]
  multiple_details <- vapply(multiple_drift, function(name) {
    actual <- actual_multiple[[name]]
    expected <- expected_multiple[[name]]
    if (is.null(actual)) actual <- character(0L)
    if (is.null(expected)) expected <- character(0L)
    paste0(
      "`", name, "` has [", paste(actual, collapse = ", "),
      "]; its exception records [", paste(expected, collapse = ", "), "]."
    )
  }, character(1L))
  non_snake_case <- audited[
    audited != "..." & !grepl("^[a-z][a-z0-9]*(_[a-z0-9]+)*$", audited)
  ]
  reasoned_groups <- c(argument_concepts, multiple_spelling_exceptions)

  expect_true(
    all(vapply(reasoned_groups, function(x) {
      is.character(x$reason) && length(x$reason) == 1L &&
        nzchar(x$reason) && !grepl("[\r\n]", x$reason)
    }, logical(1L))) &&
      all(nzchar(non_snake_case_exceptions)) &&
      !any(grepl("[\r\n]", non_snake_case_exceptions)),
    info = paste(
      "Every argument concept and spelling exception needs a one-line",
      "reason."
    )
  )
  expect_identical(
    sort(unique(assignments[duplicated(assignments)])),
    character(0L),
    info = "Each audited formal must belong to exactly one argument concept."
  )

  expect_identical(
    unexpected,
    character(0L),
    info = paste(
      "Unclassified formal(s):", paste(unexpected, collapse = ", "),
      "Classify each under an existing concept or add a new reasoned concept."
    )
  )
  expect_identical(
    unused,
    character(0L),
    info = paste(
      "Remove unused spelling(s) from concept(s):",
      paste(unused_concepts, collapse = ", ")
    )
  )
  expect_identical(
    actual_multiple,
    expected_multiple,
    info = paste(
      paste(multiple_details, collapse = " "),
      "A concept may have multiple spellings only when that exact set is a",
      "reasoned exception. Reject the new spelling or deliberately revise",
      "the exception."
    )
  )
  name_audit_expect_allowlist(
    non_snake_case,
    non_snake_case_exceptions,
    info = paste(
      "Every non-snake-case formal needs a recorded upstream/frozen",
      "reason."
    )
  )
  name_audit_expect_allowlist(
    character(0L),
    NULL,
    info = paste(
      "The non-snake-case exception check must be shape-stable when no",
      "exceptions remain."
    )
  )
  expect_true(
    all(required %in% observed),
    info = "The canonical cross-layer formals must remain present."
  )
})
