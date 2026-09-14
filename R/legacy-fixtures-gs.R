#' Legacy gene-set baseline fixtures
#'
#' These non-exported implementations are retained because the golden baseline
#' captured at `752481f` verifies current behavior against them.
#'
#' @name deprecated-gs
#' @keywords internal
NULL

#' Deprecated: list bundled reference gene-set databases
#'
#' Superseded by [gsdb_list()]. This shim reorders and relabels
#' [gsdb_list()]'s output to match the old column layout.
#'
#' @param toolkit_dir Character. Accepted and ignored: the new API resolves
#'   bundled data solely through `system.file(package = "bulkiRNA")`, so there
#'   is no toolkit checkout path to point at.
#' @return A data frame with columns `database`, `name`, `bundled`,
#'   `description`, `species`.
#' @keywords internal
list_reference_dbs <- function(toolkit_dir = NULL) {
  .Deprecated("gsdb_list")
  out <- gsdb_list()
  out$species[!out$bundled] <- "(not bundled)"
  out[c("database", "name", "bundled", "description", "species")]
}

#' Deprecated: load a bundled reference gene-set database
#'
#' Superseded by [gsdb_load()]. Returns the legacy `list(T2G =, T2N =,
#' source =, created =)` shape a `clusterProfiler::GSEA()` call expects,
#' rebuilt from the [gs_db()] that [gsdb_load()] now returns.
#'
#' @param database Character. Database key; see [gsdb_list()].
#' @param species Character. Species directory name (default:
#'   `"Mus_musculus"`). Forwarded as-is; [gsdb_load()] accepts both the
#'   underscore and spaced forms.
#' @param toolkit_dir Character. Accepted and ignored, as in
#'   `list_reference_dbs()`.
#' @param rebuild Logical. If `TRUE`, [gsdb_load()] always errors: an
#'   installed package never ships the raw source files a rebuild needs, even
#'   when this shim happens to run from a source checkout. The old function
#'   could actually rebuild in that situation; this is a documented behaviour
#'   change, not a bug in the shim.
#' @return List with `T2G`, `T2N`, `source`, `created`, matching the object
#'   the pre-refactor `load_reference_db()` returned.
#' @keywords internal
load_reference_db <- function(
    database,
    species = "Mus_musculus",
    toolkit_dir = NULL,
    rebuild = FALSE
) {
  .Deprecated("gsdb_load")
  db <- gsdb_load(database, species = species, rebuild = rebuild)
  legacy <- .gsdb_as_t2g(db)

  # Best-effort recovery of the original `source`/`created` fields: read them
  # back off the processed RDS this database resolves to, when there is one.
  reg <- .gsdb_registry()
  entry <- reg[[database]]
  src <- NULL
  crt <- NULL
  if (!is.null(entry)) {
    sp_dir <- .gsdb_species_dir(species)
    path <- .gsdb_extdata(entry$dir, "processed", sp_dir, entry$rds_file)
    if (nzchar(path) && file.exists(path)) {
      raw <- tryCatch(readRDS(path), error = function(e) NULL)
      src <- raw$source
      crt <- raw$created
    }
  }
  legacy$source <- src %||% attr(db, "database_label")
  legacy$created <- crt %||% Sys.time()
  legacy
}

#' Deprecated: filter a T2G/T2N pair by gene-set size
#'
#' Superseded by the internal `filter_by_size()` in `R/gs-db.R`, which
#' operates on a [gs_db()] rather than a `list(T2G =, T2N =)` pair. The 5/500
#' defaults live here; the internal version defaults to `NULL`/`NULL` ("no
#' bound").
#'
#' **Naming collision, not resolved here:** `R/gs-db.R` already defines an
#' internal (non-exported, but still top-level) `filter_by_size(db, min_size
#' = NULL, max_size = NULL, verbose = FALSE)`. A package namespace has one
#' binding per name, so this shim's own `filter_by_size(result, min_size =
#' 5, max_size = 500)` would silently shadow -- or be shadowed by -- that
#' internal one purely based on file collation order, corrupting whichever
#' one loses. This implementation therefore does **not** call the internal
#' `filter_by_size()` at all; it reimplements the old
#' `build_reference_databases.R` size filter directly on the `T2G`/`T2N`
#' data frames. The collision itself (two distinct top-level `filter_by_size
#' <- function(...)` bindings in one package) is a package-level defect that
#' this shim cannot fix without editing `R/gs-db.R`; see the handback report.
#'
#' @param result List with `T2G` (`gs_name`, `gene_symbol`) and `T2N`
#'   (`gs_name`, `description`).
#' @param min_size Integer. Minimum genes per set (default: 5).
#' @param max_size Integer. Maximum genes per set (default: 500).
#' @return List with filtered `T2G` and `T2N`, same shape as `result`.
#' @keywords internal
filter_by_size <- function(result, min_size = 5, max_size = 500) {
  # Names a public path, not an internal. There is no exported successor
  # function: size filtering became an *argument* on every provider, so the
  # replacement for a separate filtering step is not calling one. (The former
  # target, an internal also called `filter_by_size()`, no longer exists under
  # that name -- it was renamed `.gs_filter_size()` precisely because sharing a
  # name with this export silently shadowed it.)
  .Deprecated(msg = paste(
    "`filter_by_size()` is deprecated. Size filtering is now an argument on",
    "the providers: gsdb_msigdb(min_size=, max_size=), gsdb_load(...),",
    "gsdb_from_file(...)."))
  gs_sizes <- table(result$T2G$gs_name)
  keep <- names(gs_sizes[gs_sizes >= min_size & gs_sizes <= max_size])
  result$T2G <- result$T2G[result$T2G$gs_name %in% keep, , drop = FALSE]
  result$T2N <- result$T2N[result$T2N$gs_name %in% keep, , drop = FALSE]
  result
}

#' Local, collision-free set-size filter for the shims in this file
#'
#' Used by `parse_gmx()` and `parse_mitoxplorer()` instead of calling
#' `filter_by_size()` by name, for the reason documented on that function:
#' the identifier is ambiguous package-wide.
#'
#' @param db A [gs_db()].
#' @param min_size,max_size Integer(1) or `NULL`.
#' @return A `gs_db` with out-of-range sets removed.
#' @keywords internal
.dep_gs_filter_size <- function(db, min_size, max_size) {
  if (is.null(min_size) && is.null(max_size)) return(db)
  sizes <- vapply(db, length, integer(1L))
  lo <- if (is.null(min_size)) 1L else as.integer(min_size)
  hi <- if (is.null(max_size)) Inf else as.numeric(max_size)
  db[sizes >= lo & sizes <= hi]
}

#' Deprecated: convert a named gene-set list to TERM2GENE format
#'
#' Superseded by [gsdb_register()], which accepts a named list of gene sets and
#' returns a registered `gs_db`. This shim additionally converts that database
#' to the old TERM2GENE shape. Set names with no genes are dropped and genes
#' are de-duplicated within a set; the old function preserved empty/duplicate
#' rows verbatim.
#'
#' @param geneset_list Named list where names are pathway IDs and values are
#'   character vectors of genes.
#' @param gs_col Character. Name for the gene-set column (default:
#'   `"gs_name"`).
#' @param gene_col Character. Name for the gene column (default:
#'   `"gene_symbol"`).
#' @return Data frame in TERM2GENE format with columns `gs_col`, `gene_col`.
#' @keywords internal
list_to_term2gene <- function(
    geneset_list,
    gs_col = "gs_name",
    gene_col = "gene_symbol"
) {
  .Deprecated("gsdb_register")
  if (!is.list(geneset_list) || is.null(names(geneset_list))) {
    stop("`geneset_list` must be a named list; name each element with its ",
         "pathway ID.", call. = FALSE)
  }
  db <- gs_db(geneset_list, database = "legacy_list", species = "Mus musculus")
  t2g <- .gsdb_as_t2g(db)$T2G
  colnames(t2g) <- c(gs_col, gene_col)
  t2g
}

#' Deprecated: parse a GMX file to gene sets
#'
#' Superseded by [gsdb_from_file()]. This shim calls the same internal GMX
#' parser [gsdb_from_file()] uses, then rebuilds the old
#' `list(T2G =, T2N =, source =, created =)` shape.
#'
#' @param file Character. Path to GMX file.
#' @param prefix Character. Prefix to add to gene-set names (optional).
#' @param min_size Integer. Minimum genes per set to retain (default: 5).
#' @param max_size Integer. Maximum genes per set to retain (default: 500).
#' @return List with `T2G`, `T2N`, `source`, `created`.
#' @keywords internal
parse_gmx <- function(file, prefix = NULL, min_size = 5, max_size = 500) {
  .Deprecated("gsdb_from_file")
  if (!file.exists(file)) {
    stop("`file` not found: ", file,
         ". Supply the path to an existing GMX file.", call. = FALSE)
  }
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parsed <- .gsdb_parse_gmx(lines)
  ids <- names(parsed$sets)
  if (!is.null(prefix)) {
    ids <- paste0(prefix, "_", ids)
    names(parsed$sets) <- ids
    names(parsed$labels) <- ids
  }
  db <- gs_db(parsed$sets, database = "legacy_gmx", species = "Mus musculus",
              pathway_names = parsed$labels)
  db <- .dep_gs_filter_size(db, min_size, max_size)
  if (!length(db)) {
    warning("No gene sets passed size filters")
    return(list(
      T2G = data.frame(gs_name = character(), gene_symbol = character(),
                        stringsAsFactors = FALSE),
      T2N = data.frame(gs_name = character(), description = character(),
                        stringsAsFactors = FALSE),
      source = basename(file),
      created = Sys.time()
    ))
  }
  legacy <- .gsdb_as_t2g(db)
  legacy$source <- basename(file)
  legacy$created <- Sys.time()
  legacy
}

#' Deprecated: parse a mitoXplorer gene-function table
#'
#' Superseded by the internal `.gsdb_parse_mitoxplorer()` in
#' `R/gsdb-rebuild.R`, which this shim calls directly, then converts back to
#' the old `list(T2G =, T2N =, stats =, source =, created =)` shape.
#'
#' @param file Character. Path to mitoXplorer `mouse_gene_function.txt`.
#' @param prefix Character. Prefix for gene-set names (default:
#'   `"MITOXPLORER"`).
#' @param gene_col Character. Column name for gene symbols (default:
#'   `"MGI_symbol"`).
#' @param process_col Character. Column name for pathway/process (default:
#'   `"mito_process"`).
#' @param min_size Integer. Minimum genes per set to retain (default: 5).
#' @param max_size Integer. Maximum genes per set to retain (default: 500).
#' @return List with `T2G`, `T2N`, `stats`, `source`, `created`.
#' @keywords internal
parse_mitoxplorer <- function(
    file,
    prefix = "MITOXPLORER",
    gene_col = "MGI_symbol",
    process_col = "mito_process",
    min_size = 5,
    max_size = 500
) {
  # The parser itself is internal, so naming it told users to call something
  # they cannot reach. Point at the public route to the same sets instead.
  .Deprecated(msg = paste(
    "`parse_mitoxplorer()` is deprecated. Load the bundled sets with",
    "gsdb_load(\"mitoxplorer\"), or parse an arbitrary file with",
    "gsdb_from_file(). Rebuilding from a raw mitoXplorer export is a",
    "maintainer task and has no exported equivalent."))
  if (!file.exists(file)) {
    stop("`file` not found: ", file,
         ". Supply the path to an existing mitoXplorer file.", call. = FALSE)
  }
  db <- .gsdb_parse_mitoxplorer(file, species = "Mus musculus",
                                prefix = prefix, gene_col = gene_col,
                                process_col = process_col)
  db <- .dep_gs_filter_size(db, min_size, max_size)
  if (!length(db)) {
    warning("No gene sets passed size filters")
    return(list(
      T2G = data.frame(gs_name = character(), gene_symbol = character(),
                        stringsAsFactors = FALSE),
      T2N = data.frame(gs_name = character(), description = character(),
                        stringsAsFactors = FALSE),
      source = basename(file),
      created = Sys.time()
    ))
  }
  legacy <- .gsdb_as_t2g(db)
  legacy$stats <- list(
    n_genesets = length(db),
    n_genes = length(unique(unlist(db, use.names = FALSE))),
    source_file = basename(file)
  )
  legacy$source <- "mitoXplorer3.0"
  legacy$created <- Sys.time()
  legacy
}

#' Deprecated: convert human gene symbols in a T2G table to mouse orthologs
#'
#' Superseded by the internal `.gsdb_human_to_mouse()` in
#' `R/gsdb-rebuild.R`, which operates on a [gs_db()] rather than a raw T2G
#' data frame. This shim round-trips through that shape.
#'
#' @param T2G Data frame with gene-set/gene mappings (`gs_name`,
#'   `gene_symbol` columns).
#' @param drop_unmapped Logical. Accepted, but ignored when `FALSE`: the new
#'   ortholog mapper always drops genes it cannot map, so there is no
#'   supported way to keep them. A warning is raised in that case.
#' @param verbose Logical. Print mapping statistics (default: TRUE).
#' @return `T2G` with mouse gene symbols.
#' @keywords internal
convert_human_to_mouse <- function(T2G, drop_unmapped = TRUE, verbose = TRUE) {
  # The mapper is internal; the user-facing equivalent is asking a provider for
  # mouse sets from the human MSigDB, which maps orthologs itself.
  .Deprecated(msg = paste(
    "`convert_human_to_mouse()` is deprecated. Ask for the species you want",
    "directly: gsdb_msigdb(species = \"Mus musculus\", db_species = \"HS\")",
    "maps human sets to mouse orthologs as it loads them."))
  if (isFALSE(drop_unmapped)) {
    warning("`drop_unmapped = FALSE` has no effect: the new ortholog mapper ",
            "always drops genes it cannot map.", call. = FALSE)
  }
  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"
  sets <- split(as.character(T2G[[gene_col]]), as.character(T2G$gs_name))
  db <- gs_db(sets, database = "legacy_human", species = "Homo sapiens")
  mouse_db <- .gsdb_human_to_mouse(db, verbose = verbose)
  out <- .gsdb_as_t2g(mouse_db)$T2G
  colnames(out)[colnames(out) == "gene_symbol"] <- gene_col
  out
}

#' Deprecated: an empty normalized GSEA tibble
#'
#' A `gs_result` comes back from [gs_test()]; filter one to zero rows when an
#' empty result is needed. There is no exported constructor, by design. The
#' old function returned a zero-row tibble in `normalize_gsea_results()`'s
#' bespoke schema (`pathway_id`,
#' `pathway_name`, `database`, `contrast`, `NES`, `pvalue`, `padj`,
#' `set_size`, `leading_edge_size`, `gene_ratio`, `core_enrichment`,
#' `genes_full_set`, `direction`, `neg_log_padj`); this shim instead returns a
#' zero-row `gs_result`, the schema `gs_test()` now produces directly. The
#' column set is therefore genuinely different, not aliased -- that is the
#' point of "superseded by", not "identical to".
#'
#' @return A zero-row `gs_result`.
#' @keywords internal
empty_gsea_tibble <- function() {
  .Deprecated(msg = paste(
    "`empty_gsea_tibble()` is deprecated. gs_test() returns a gs_result;",
    "filter its result to zero rows when an empty result is needed. There is",
    "no exported constructor, by design."))
  core <- .gs_empty_core()
  core$database <- character(0)
  core$contrast <- character(0)
  core$method <- character(0)
  core$stat_type <- character(0)
  core$direction <- character(0)
  gs_result(core)
}
#' Drop columns the legacy shims never returned
#'
#' The shims exist so that unmigrated callers keep the old shape, and such a
#' caller may index the table positionally. Additions to `gs_result` therefore
#' stop at this boundary: anyone who wants `log2err` wants [gs_test()].
#'
#' @param x A `gs_result`.
#' @return `x` without the post-legacy columns.
#' @keywords internal
.drop_legacy_extras <- function(x) {
  x[, setdiff(names(x), "log2err"), drop = FALSE]
}


#' Deprecated: run preranked GSEA against an MSigDB collection
#'
#' Superseded by [gs_ranks()] + [gs_test()]. The old `run_gsea()` returned a
#' `clusterProfiler` S4 `gseaResult`, which carried both the ranked vector
#' (`@geneList`) and the tested gene sets (`@geneSets`) alongside its results
#' table. A [gs_result()] carries neither, so this shim attaches them as the
#' `ranks` and `gene_sets` attributes that [gs_plot_running()] falls back to
#' when its own `ranks`/`db` arguments are absent -- do not drop this when
#' touching this function again.
#'
#' @param DE_results A data frame of DE results. Must have gene identifiers
#'   as rownames and a column named by `rank_metric`.
#' @param rank_metric Character. Column of `DE_results` to rank by (default:
#'   `"t"`).
#' @param species Character. Species for MSigDB gene sets (default:
#'   `"Mus musculus"`).
#' @param db_species Character. MSigDB source database (`"MM"`/`"HS"`), or
#'   `NULL` (default) for the old "legacy" mode, which this shim maps to
#'   `"HS"` -- the same human-database-with-ortholog-mapping behaviour the
#'   old species-only branch used.
#' @param collection Character. MSigDB collection code (default: `"H"`).
#' @param subcollection Character. MSigDB subcollection code (default: `""`,
#'   meaning "no subcollection"; mapped to `NULL` for [gsdb_msigdb()]).
#' @param pvalue_cutoff Numeric. Nominal p-value cutoff applied to the result
#'   after testing (default: 1, meaning no filtering).
#' @param padj_method Character. Accepted and ignored: [gs_test()]'s fgsea
#'   adapter always adjusts with `"BH"` (what `"fdr"` aliases to in
#'   [stats::p.adjust()]), so a mismatched value only warns.
#' @param nperm Integer. Forwarded to [gs_test()] as `n_perm_simple` (default:
#'   100000).
#' @param seed Integer. Random seed set before testing (default: 123).
#' @return A [gs_result()] with `ranks` and `gene_sets` attributes set.
#' @keywords internal
run_gsea <- function(
    DE_results,
    rank_metric   = "t",
    species       = "Mus musculus",
    db_species    = NULL,
    collection    = "H",
    subcollection = "",
    pvalue_cutoff = 1,
    padj_method   = "fdr",
    nperm         = 100000,
    seed          = 123
) {
  .Deprecated("gs_ranks() and gs_test()")
  if (!is.data.frame(DE_results)) {
    stop("`DE_results` must be a data frame.", call. = FALSE)
  }
  if (is.null(rownames(DE_results))) {
    stop("`DE_results` must have rownames (gene identifiers).",
         call. = FALSE)
  }
  if (!rank_metric %in% colnames(DE_results)) {
    stop("`rank_metric` names a column absent from `DE_results`; got ",
         encodeString(rank_metric, quote = "\""), ". Choose one of: ",
         paste0("`", colnames(DE_results), "`", collapse = ", "), ".",
         call. = FALSE)
  }
  if (!identical(padj_method, "fdr")) {
    warning("`padj_method` is ignored: gs_test() always adjusts with \"BH\".",
            call. = FALSE)
  }

  ranks <- gs_ranks(DE_results, metric = rank_metric)
  sub <- if (is.null(subcollection) || !nzchar(subcollection)) {
    NULL
  } else {
    subcollection
  }
  db <- gsdb_msigdb(species = species, collection = collection,
                    subcollection = sub, db_species = db_species %||% "HS")
  res <- gs_test(ranks, db, n_perm_simple = nperm, seed = seed)
  # A shim's contract is the old shape. `gs_test()` now carries fgsea's
  # `log2err`, which is useful and which no legacy caller asked for -- and an
  # unmigrated caller may index this table positionally. Drop it here; anyone
  # who wants the uncertainty bound wants `gs_test()`.
  res <- .drop_legacy_extras(res)
  if (!is.null(pvalue_cutoff) && pvalue_cutoff < 1) {
    res <- gs_filter(res, p_value = pvalue_cutoff)
  }
  attr(res, "ranks") <- ranks
  attr(res, "gene_sets") <- db
  res
}

#' Deprecated: normalize a GSEA result to a standard tibble
#'
#' Superseded by [gs_result()] -- which [gs_test()] (and this file's
#' `run_gsea()` shim) now returns directly, so there is usually nothing left
#' to normalize. This shim accepts a [gs_result()] or a plain data frame and
#' stamps in `database`/`contrast`, matching the old function's job of
#' attaching those labels to a per-database result.
#'
#' @param gsea_obj A [gs_result()] or data frame; `NULL` returns
#'   `empty_gsea_tibble()`.
#' @param database Character. Name of the source database.
#' @param contrast Character. Name of the contrast.
#' @param padj_cutoff Numeric. Filter to `padj < padj_cutoff` (default: 1, no
#'   filter).
#' @param format_names Logical. Apply [format_pathway_name()] to
#'   `pathway_name` (default: TRUE).
#' @param max_name_length Integer. Truncate `pathway_name` to this length
#'   (default: 80).
#' @param atlas_universe Character vector or `NULL`. Accepted and ignored,
#'   with a warning when supplied: a [gs_result()] carries no gene-set
#'   membership, so there is nothing to intersect against a universe. Use
#'   [gs_leading_edge()] with the original `gs_db` instead.
#' @return A [gs_result()].
#' @keywords internal
normalize_gsea_results <- function(
    gsea_obj,
    database,
    contrast,
    padj_cutoff = 1,
    format_names = TRUE,
    max_name_length = 80,
    atlas_universe = NULL
) {
  .Deprecated("gs_result, produced directly by gs_test()")
  if (is.null(gsea_obj)) {
    warning(sprintf(
      "NULL gsea_obj provided for database '%s', returning empty tibble",
      database
    ))
    return(empty_gsea_tibble())
  }
  if (inherits(gsea_obj, "gs_result")) {
    out <- gsea_obj
  } else if (is.data.frame(gsea_obj)) {
    out <- gs_result(gsea_obj, database = database, contrast = contrast)
  } else {
    warning(sprintf("Unexpected object type for database '%s': %s",
                     database, class(gsea_obj)[1]))
    return(empty_gsea_tibble())
  }
  if (!nrow(out)) {
    message(sprintf("No results in gsea_obj for database '%s'", database))
    return(empty_gsea_tibble())
  }

  out$database <- database
  out$contrast <- contrast

  if (!is.null(atlas_universe)) {
    warning("`atlas_universe` is ignored: a gs_result carries no gene-set ",
            "membership. Use gs_leading_edge() with the original gs_db.",
            call. = FALSE)
  }
  if (!is.null(padj_cutoff) && padj_cutoff < 1) {
    out <- gs_filter(out, padj = padj_cutoff)
    if (!nrow(out)) {
      message(sprintf("No significant results (padj < %.2f) for database '%s'",
                       padj_cutoff, database))
      return(empty_gsea_tibble())
    }
  }

  if (isTRUE(format_names)) {
    out$pathway_name <- format_pathway_name(out$pathway_id)
  }
  too_long <- nchar(out$pathway_name) > max_name_length
  if (any(too_long)) {
    out$pathway_name[too_long] <- paste0(
      substr(out$pathway_name[too_long], 1, max_name_length - 3), "..."
    )
  }
  .drop_legacy_extras(out)
}

#' Deprecated: run a multi-database GSEA analysis and save its plots
#'
#' Superseded by looping [gs_ranks()]/[gsdb_msigdb()]/[gs_test()] over the
#' databases of interest and rendering with [gs_plot_dot()], [gs_plot_bar()],
#' [gs_plot_running()] and [gs_save()]. This shim reproduces the old formals
#' and default database list, and forwards each database's run to the
#' `run_gsea()` shim above (so the `ranks`/`gene_sets` attributes it sets are
#' available for the running-sum plot) rather than sourcing helper scripts --
#' there is no `source()` step in the new API.
#'
#' Simplifications relative to the old pipeline, made deliberately thin per
#' the module brief: the up/down/facet dotplot triptych collapses to one
#' [gs_plot_dot()] call with `direction = "both"`; `sample_annotation` and
#' `sample_order` (the per-database sample x pathway heatmap) are accepted but
#' ignored, with a warning, because [gs_plot_heatmap()] draws directly from a
#' [gs_result()] or [gs_matrix()] and has no equivalent single-NES-vector
#' knockoff heatmap.
#'
#' @param de_table Differential expression results table.
#' @param analysis_name Name of the analysis, used in plot filenames.
#' @param rank_metric Column to rank genes by (default: `"t"`).
#' @param species Species for MSigDB gene sets (default: `"Mus musculus"`).
#' @param n_pathways Number of top pathways to display (default: 30).
#' @param padj_cutoff Adjusted p-value cutoff for highlighting (default:
#'   0.05).
#' @param save_plots Logical. Save generated plots (default: TRUE).
#' @param output_dir Directory to save plots.
#' @param databases List of database configs (default: NULL, uses the
#'   predefined MSigDB set).
#' @param nperm Number of permutations for GSEA (default: 100000).
#' @param pvalue_cutoff P-value cutoff for storing GSEA results (default: 1,
#'   stores ALL pathways).
#' @param sample_annotation Accepted and ignored (with a warning when
#'   supplied); see Details.
#' @param sample_order Accepted and ignored (with a warning when supplied);
#'   see Details.
#' @param helper_root Accepted and ignored: the new API needs no
#'   `source()`-based helper resolution.
#' @return List of [gs_result()] objects for each database, invisibly.
#' @keywords internal
run_gsea_analysis <- function(
    de_table,
    analysis_name,
    rank_metric = "t",
    species = "Mus musculus",
    n_pathways = 30,
    padj_cutoff = 0.05,
    save_plots = TRUE,
    output_dir = "./GSEA_Plots",
    databases = NULL,
    nperm = 100000,
    pvalue_cutoff = 1,
    sample_annotation = NULL,
    sample_order = NULL,
    helper_root = NULL
) {
  .Deprecated(paste("gs_ranks(), gsdb_msigdb(), gs_test(), gs_plot_dot(),",
                     "gs_plot_bar(), gs_plot_running() and gs_save()"))
  if (!is.null(sample_annotation) || !is.null(sample_order)) {
    warning("`sample_annotation`/`sample_order` (the per-database sample x ",
            "pathway heatmap) are not reproduced by this shim; use ",
            "gs_plot_heatmap() directly on the gs_result.", call. = FALSE)
  }

  dbsp <- if (grepl("sapiens", species, ignore.case = TRUE)) "HS" else "MM"
  if (is.null(databases)) {
    databases <- list(
      hallmark = list(name = "Hallmark", db_species = dbsp, collection = "H", subcollection = ""),
      canon    = list(name = "Canonical Pathways", db_species = dbsp, collection = "C2", subcollection = "CP"),
      gobp     = list(name = "GO BP", db_species = dbsp, collection = "C5", subcollection = "GO:BP"),
      gomf     = list(name = "GO MF", db_species = dbsp, collection = "C5", subcollection = "GO:MF"),
      gocc     = list(name = "GO CC", db_species = dbsp, collection = "C5", subcollection = "GO:CC"),
      kegg     = list(name = "KEGG", db_species = dbsp, collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
      reactome = list(name = "Reactome", db_species = dbsp, collection = "C2", subcollection = "CP:REACTOME"),
      wiki     = list(name = "WikiPath", db_species = dbsp, collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
      cgp      = list(name = "Chem-Genetic Perturbations", db_species = dbsp, collection = "C2", subcollection = "CGP"),
      tf       = list(name = "GTRD", db_species = dbsp, collection = "C3", subcollection = "TFT:GTRD")
    )
  }

  gsea_results <- list()

  for (db_name in names(databases)) {
    cfg <- databases[[db_name]]
    message("  ", cfg$name)

    res <- tryCatch(
      run_gsea(
        DE_results = de_table,
        rank_metric = rank_metric,
        species = species,
        db_species = cfg$db_species,
        collection = cfg$collection,
        subcollection = cfg$subcollection,
        nperm = nperm,
        pvalue_cutoff = pvalue_cutoff
      ),
      error = function(e) {
        warning(e)
        NULL
      }
    )

    gsea_results[[db_name]] <- res
    if (is.null(res) || !nrow(res)) next
    if (!isTRUE(save_plots)) next

    db_dir <- file.path(output_dir, db_name)
    ensure_dir(db_dir)

    gs_save(
      gs_plot_dot(res, top_n = n_pathways, direction = "both",
                  highlight = padj_cutoff,
                  title = sprintf("%s %s", analysis_name, db_name)),
      file.path(db_dir, sprintf("%s_%s_dot", analysis_name, db_name))
    )
    gs_save(
      gs_plot_bar(res, top_n = n_pathways, highlight = padj_cutoff,
                  title = sprintf("%s %s NES", analysis_name, db_name)),
      file.path(db_dir, sprintf("%s_%s_nes_bar", analysis_name, db_name))
    )

    top5 <- gs_top(res, n = 5L, by = "padj", per = character(0))$pathway_id
    gs_save(
      gs_plot_running(res, pathways = top5, top_n = length(top5)),
      file.path(db_dir, sprintf("%s_%s_running_sum", analysis_name, db_name))
    )
  }

  invisible(gsea_results)
}
