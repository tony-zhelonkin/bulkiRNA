# Contract tests for the frozen GSEA/db deprecation shims in R/deprecated-gs.R.
# Each test checks (a) the old formals still work, (b) a deprecation warning
# fires, and (c) the forwarded new-API result is sane.

test_that("list_reference_dbs() forwards to gsdb_list() with the old layout", {
  out <- suppressWarnings(list_reference_dbs())
  expect_identical(names(out), c("database", "name", "bundled",
                                 "description", "species"))
  expect_true(out$bundled[out$database == "mitopathways"])
  expect_warning(list_reference_dbs(), class = "deprecatedWarning")
  # toolkit_dir is accepted and ignored
  expect_warning(list_reference_dbs(toolkit_dir = "/nonexistent"),
                 class = "deprecatedWarning")
})

test_that("load_reference_db() returns the legacy T2G/T2N/source/created shape", {
  out <- suppressWarnings(load_reference_db("mitopathways"))
  expect_type(out, "list")
  expect_true(all(c("T2G", "T2N", "source", "created") %in% names(out)))
  expect_true(all(c("gs_name", "gene_symbol") %in% names(out$T2G)))
  expect_true(all(c("gs_name", "description") %in% names(out$T2N)))
  expect_warning(load_reference_db("mitopathways"), class = "deprecatedWarning")
})

test_that("load_reference_db() rebuild=TRUE now always errors", {
  expect_error(suppressWarnings(load_reference_db("mitopathways", rebuild = TRUE)),
               "source checkout")
})

test_that("filter_by_size() filters a T2G/T2N pair with the old defaults", {
  # RESOLVED COLLISION. C1 found that R/gs-db.R defined its own top-level
  # `filter_by_size(db, min_size = NULL, max_size = NULL, verbose = FALSE)`.
  # A namespace has one binding per name and collation puts "gs-db.R" after
  # "deprecated-gs.R", so that definition -- not this shim -- was what
  # `bulkiRNA::filter_by_size()` resolved to: no `.Deprecated()` warning, and
  # the frozen 5/500 defaults silently replaced by NULL/NULL. C1 correctly
  # refused to edit a file it did not own and skipped this test instead. The
  # integrator renamed the internal to `.gs_filter_size()`, so the shim now
  # owns its frozen name and the assertions below are reachable.
  t2g <- data.frame(
    gs_name = c(rep("SMALL", 2), rep("BIG", 8)),
    gene_symbol = paste0("G", 1:10),
    stringsAsFactors = FALSE
  )
  t2n <- data.frame(gs_name = c("SMALL", "BIG"),
                    description = c("s", "b"), stringsAsFactors = FALSE)
  out <- suppressWarnings(filter_by_size(list(T2G = t2g, T2N = t2n),
                                         min_size = 5, max_size = 500))
  expect_identical(unique(out$T2G$gs_name), "BIG")
  expect_identical(out$T2N$gs_name, "BIG")
  expect_warning(filter_by_size(list(T2G = t2g, T2N = t2n)),
                 class = "deprecatedWarning")
})

test_that("list_to_term2gene() renames columns and forwards through gs_db", {
  sets <- list(SET_A = c("A1", "A2"), SET_B = c("B1"))
  out <- suppressWarnings(list_to_term2gene(sets))
  expect_identical(names(out), c("gs_name", "gene_symbol"))
  expect_setequal(out$gs_name, c("SET_A", "SET_A", "SET_B"))
  out2 <- suppressWarnings(list_to_term2gene(sets, gs_col = "term",
                                             gene_col = "gene"))
  expect_identical(names(out2), c("term", "gene"))
  expect_error(suppressWarnings(list_to_term2gene(unname(sets))), "named list")
})

test_that("parse_gmx() parses a GMX file into the legacy shape", {
  f <- tempfile(fileext = ".gmx")
  writeLines(c(
    "desc1\tdesc2",
    "SET_ONE\tSET_TWO",
    paste(paste0("g", 1:6), paste0("h", 1:6), sep = "\t")
  ), f)
  out <- suppressWarnings(parse_gmx(f, min_size = 1, max_size = 10))
  expect_setequal(unique(out$T2G$gs_name), c("SET_ONE", "SET_TWO"))
  expect_identical(out$source, basename(f))
  expect_s3_class(out$created, "POSIXct")
  expect_warning(parse_gmx(f, min_size = 1, max_size = 10),
                 class = "deprecatedWarning")
  expect_error(suppressWarnings(parse_gmx(tempfile())), "not found")
})

test_that("parse_mitoxplorer() parses a gene-function table with stats", {
  f <- tempfile(fileext = ".txt")
  d <- data.frame(
    MGI_symbol = paste0("G", 1:12),
    mito_process = rep(c("Metabolism", "Transport"), each = 6),
    stringsAsFactors = FALSE
  )
  write.table(d, f, sep = "\t", row.names = FALSE, quote = FALSE)
  out <- suppressWarnings(parse_mitoxplorer(f, min_size = 1, max_size = 100))
  expect_setequal(unique(out$T2G$gs_name),
                  c("MITOXPLORER_METABOLISM", "MITOXPLORER_TRANSPORT"))
  expect_identical(out$stats$n_genesets, 2L)
  expect_identical(out$source, "mitoXplorer3.0")
})

test_that("convert_human_to_mouse() forwards through the ortholog mapper", {
  skip_if_not_installed("homologene")
  t2g <- data.frame(
    gs_name = c("SET_A", "SET_A"),
    gene_symbol = c("TP53", "EGFR"),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(convert_human_to_mouse(t2g, verbose = FALSE))
  expect_true(all(c("gs_name", "gene_symbol") %in% names(out)))
  expect_warning(convert_human_to_mouse(t2g, drop_unmapped = FALSE,
                                        verbose = FALSE),
                 "always drops")
})

test_that("empty_gsea_tibble() is a zero-row gs_result, not the old schema", {
  out <- suppressWarnings(empty_gsea_tibble())
  expect_s3_class(out, "gs_result")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("pathway_id", "database", "contrast", "stat",
                    "stat_type", "direction") %in% names(out)))
})

test_that("run_gsea() sets both `ranks` and `gene_sets` attributes", {
  # >= 10 genes per set: gs_test()'s fgsea adapter defaults to min_size = 10,
  # matching the old clusterProfiler::GSEA() default this shim preserves.
  de <- data.frame(t = seq(20, -20, length.out = 40),
                    row.names = paste0("G", 1:40))
  db <- structure(
    list(SET_UP = paste0("G", 1:15), SET_DOWN = paste0("G", 26:40)),
    pathway_names = c(SET_UP = "Up set", SET_DOWN = "Down set"),
    database = "testdb", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
  called_with <- NULL
  testthat::with_mocked_bindings(
    gsdb_msigdb = function(...) { called_with <<- list(...); db },
    {
      res <- suppressWarnings(run_gsea(de, seed = 1, nperm = 1000,
                                       collection = "H"))
      expect_s3_class(res, "gs_result")
      # the contract test-gs-plot-running.R documents: gs_plot_running()
      # falls back to these attributes when its own args are absent.
      expect_true(is.numeric(attr(res, "ranks")))
      expect_identical(names(attr(res, "ranks")),
                       names(sort(stats::setNames(de$t, rownames(de)),
                                  decreasing = TRUE)))
      expect_identical(attr(res, "gene_sets"), db)
      expect_s3_class(gs_plot_running(res, top_n = 1), "ggplot")
    }
  )
  expect_warning(
    testthat::with_mocked_bindings(
      gsdb_msigdb = function(...) db,
      run_gsea(de, seed = 1, nperm = 1000)
    ),
    class = "deprecatedWarning"
  )
})

test_that("run_gsea() maps db_species = NULL to \"HS\" and subcollection = \"\" to NULL", {
  de <- data.frame(t = c(1, -1), row.names = c("A", "B"))
  db <- structure(
    list(SET_A = c("A", "B")),
    pathway_names = c(SET_A = "Set A"),
    database = "testdb", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
  seen <- NULL
  testthat::with_mocked_bindings(
    gsdb_msigdb = function(...) { seen <<- list(...); db },
    suppressWarnings(run_gsea(de, seed = 1, nperm = 100))
  )
  expect_identical(seen$db_species, "HS")
  expect_null(seen$subcollection)
})

test_that("normalize_gsea_results() stamps database/contrast onto a gs_result", {
  db <- structure(
    list(A = c("g1", "g2")), pathway_names = c(A = "Set A"),
    database = "d", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
  res <- bulkiRNA:::gs_result(
    data.frame(pathway_id = "A", pathway_name = "A", n_genes = 2L,
              n_genes_tested = 2L, stat = 1.5, p_value = 0.01, padj = 0.02,
              stringsAsFactors = FALSE),
    database = "old", contrast = "old", method = "fgsea", stat_type = "NES"
  )
  out <- suppressWarnings(normalize_gsea_results(res, database = "Hallmark",
                                                 contrast = "KO_vs_WT"))
  expect_identical(unique(out$database), "Hallmark")
  expect_identical(unique(out$contrast), "KO_vs_WT")
  expect_warning(normalize_gsea_results(res, database = "Hallmark",
                                        contrast = "KO_vs_WT",
                                        atlas_universe = c("g1")),
                 "atlas_universe.*ignored")

  expect_s3_class(suppressWarnings(normalize_gsea_results(NULL, database = "x",
                                                          contrast = "y")),
                 "gs_result")
})

test_that("run_gsea_analysis() reproduces the default database list and formals", {
  fmls <- names(formals(run_gsea_analysis))
  expect_identical(
    fmls,
    c("de_table", "analysis_name", "rank_metric", "species", "n_pathways",
      "padj_cutoff", "save_plots", "output_dir", "databases", "nperm",
      "pvalue_cutoff", "sample_annotation", "sample_order", "helper_root")
  )
  de <- data.frame(t = c(3, -3), row.names = c("A", "B"))
  db <- structure(
    list(SET_A = c("A", "B")), pathway_names = c(SET_A = "Set A"),
    database = "d", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
  called <- 0L
  out <- suppressWarnings(testthat::with_mocked_bindings(
    gsdb_msigdb = function(...) { called <<- called + 1L; db },
    run_gsea_analysis(de, "demo", databases = list(
      one = list(name = "One", db_species = "HS", collection = "H", subcollection = "")
    ), save_plots = FALSE)
  ))
  expect_identical(called, 1L)
  expect_named(out, "one")
  expect_s3_class(out$one, "gs_result")
})

test_that("run_gsea_analysis() warns when sample_annotation/sample_order are given", {
  de <- data.frame(t = c(1, -1), row.names = c("A", "B"))
  expect_warning(
    run_gsea_analysis(
      de, "demo", databases = list(), sample_annotation = data.frame(x = 1)
    ),
    "not reproduced"
  )
})
