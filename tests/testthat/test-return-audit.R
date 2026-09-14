return_audit_r_dir <- function() {
  for (path in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(path) && length(list.files(path, "[.]R$"))) return(path)
  }
  NULL
}

return_audit_exports <- function() {
  r_dir <- return_audit_r_dir()
  if (is.null(r_dir)) return(NULL)
  namespace <- file.path(dirname(r_dir), "NAMESPACE")
  if (!file.exists(namespace)) return(NULL)
  lines <- readLines(namespace, warn = FALSE)
  sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", lines, value = TRUE))
}

return_audit_rows <- function(names, category, expected, reason = "") {
  data.frame(
    name = names,
    category = category,
    expected = expected,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

return_audit_contract <- function() {
  rbind(
    return_audit_rows(
      c(
        "annotate_genes", "bulkirna_api", "bulkirna_check_deps",
        "bulkirna_stochastic", "coresh_chunks", "coresh_convergence",
        "coresh_loadings", "coresh_match", "coresh_search",
        "coresh_validate", "gs_to_master", "gs_validate_master"
      ),
      "tabular", "tbl_df"
    ),
    return_audit_rows(
      c("gs_coregulation", "gs_filter", "gs_read", "gs_test", "gs_top"),
      "tabular", "gs_result"
    ),
    return_audit_rows(
      c(
        "de_bfc_plot", "de_md_plot", "de_pca", "de_volcano",
        "de_volcano_grid", "gs_plot_bar", "gs_plot_dot",
        "gs_plot_heatmap", "gs_scale_fonts"
      ),
      "renderer", "ggplot"
    ),
    # Pinned as patchwork, not ggplot. A patchwork inherits "ggplot", so the
    # looser assertion would still pass if this regressed to a single plot --
    # and that regression is exactly what broke every consumer figure at 1.0.0.
    return_audit_rows("gs_plot_running", "renderer", "patchwork"),
    return_audit_rows(
      c(
        "ensure_dir", "gatom_download_refs", "gatom_save_html", "gs_save",
        "gs_write", "write_session_provenance"
      ),
      "writer", "character"
    ),
    return_audit_rows(
      c(
        "coresh_sets", "gsdb_coresh", "gsdb_from_file", "gsdb_load",
        "gsdb_msigdb", "gsdb_register"
      ),
      "domain object", "gs_db"
    ),
    return_audit_rows("gs_score", "domain object", "gs_matrix"),
    return_audit_rows("build_dge", "domain object", "DGEList"),
    return_audit_rows("gatom_de", "domain object", "gatom_de"),
    return_audit_rows("gatom_module", "domain object", "igraph"),
    return_audit_rows("gatom_refs", "domain object", "gatom_refs"),
    return_audit_rows("de_pca_3d", "interactive renderer", "plotly"),
    return_audit_rows("theme_bulki", "renderer component", "theme"),
    return_audit_rows("gsdb_list", "metadata", "data.frame"),
    return_audit_rows(c("bulki_palettes", "gs_plot_size"), "metadata", "list"),
    return_audit_rows("gsdb_info", "metadata", "list"),
    return_audit_rows("read_counts_matrix", "reader", "matrix"),
    return_audit_rows("read_metadata", "reader", "data.frame"),
    return_audit_rows(
      c(
        "coresh_labels", "entrez_to_gene", "filter_confounder_genes",
        "format_pathway_name",
        "gatom_genes", "gs_master_columns", "gs_palette", "gs_stat_types"
      ),
      "vector", "character"
    ),
    return_audit_rows("gene_to_entrez", "vector", "integer"),
    return_audit_rows("gs_ranks", "vector", "numeric"),
    return_audit_rows("gs_leading_edge", "extractor", "list"),
    return_audit_rows("gs_split", "extractor", "list")
  )
}

return_audit_probe_exceptions <- c(
  annotate_genes = paste(
    "Requires the optional AnnotationDbi and organism database; its focused",
    "test asserts the tibble result when those packages are installed."
  ),
  build_dge = paste(
    "Requires optional edgeR; test-dge.R asserts the DGEList contract in the",
    "dependency-enabled test path."
  ),
  de_pca = "Requires optional edgeR; test-de-pca.R asserts the ggplot result.",
  de_pca_3d = paste(
    "Requires optional edgeR and plotly; test-de-pca.R asserts the plotly",
    "result in the dependency-enabled path."
  ),
  de_volcano_grid = paste(
    "Requires optional patchwork; test-de-volcano.R asserts the patchwork",
    "object, which is also a ggplot."
  ),
  entrez_to_gene = paste(
    "Requires optional AnnotationDbi and an organism database; the focused",
    "gene-ID tests assert its character-vector contract."
  ),
  gatom_genes = paste(
    "Requires optional igraph; test-gatom.R asserts its character-vector",
    "contract on a module fixture."
  ),
  gatom_module = paste(
    "Requires the optional GATOM, mwcsr, and igraph stack; test-gatom.R mocks",
    "the stochastic third-party pipeline and asserts an igraph result."
  ),
  gatom_save_html = paste(
    "Requires optional GATOM plus a working pandoc executable; its focused",
    "test asserts an invisible path after a mocked HTML write."
  ),
  gene_to_entrez = paste(
    "Requires optional AnnotationDbi and an organism database; the focused",
    "gene-ID tests assert its integer-vector contract."
  ),
  gs_coregulation = paste(
    "Runs stochastic GESECA; test-gs-coregulation.R uses the discriminating",
    "600-gene fixture and asserts a gs_result."
  ),
  gs_score = paste(
    "Requires optional GSVA; test-gs-score.R asserts a gs_matrix on its",
    "dependency-enabled fixture."
  ),
  coresh_search = paste(
    "Requires optional BiocParallel; test-coresh.R uses the shipped chunk",
    "boundary and asserts the search tibble in that dependency-enabled path."
  ),
  gsdb_coresh = paste(
    "Composes the optional BiocParallel/refcache path; test-gsdb-coresh.R",
    "mocks the shipped chunk boundary and asserts a gs_db."
  ),
  gsdb_msigdb = paste(
    "The provider is network-backed; test-gsdb-msigdb.R mocks msigdbr and",
    "asserts a gs_db without making this audit network-dependent."
  )
)

return_audit_expect_class <- function(value, expected, info) {
  expect_true(inherits(value, expected), info = info)
}

return_audit_expect_allowlist <- function(observed, allowlist, info = NULL) {
  expected <- names(allowlist)
  if (is.null(expected)) expected <- character(0L)
  expect_identical(sort(as.character(observed)), sort(expected), info = info)
}

test_that("every live export has one asserted return class or one reason", {
  exports <- return_audit_exports()
  if (is.null(exports)) {
    skip(paste(
      "Source-tree-only return audit: R/ and NAMESPACE are unavailable",
      "when tests run against the installed package."
    ))
  }
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])
  contract <- return_audit_contract()

  expect_identical(
    sort(contract$name), sort(live),
    info = paste(
      "Every live export must be categorized exactly once; classify a new",
      "export before its return contract can enter the package."
    )
  )
  expect_false(any(duplicated(contract$name)))
  expect_true(all(nzchar(contract$category)))
  expect_true(all(nzchar(contract$expected)))

  exception_names <- names(return_audit_probe_exceptions)
  if (is.null(exception_names)) exception_names <- character(0L)
  expect_true(
    all(nzchar(return_audit_probe_exceptions)) &&
      !any(grepl("[\r\n]", return_audit_probe_exceptions)),
    info = "Every unprobed return needs its own one-line reason."
  )

  res <- fake_gs_result(6L)
  res$leading_edge <- rep(list(c("G1", "G2")), nrow(res))
  plot_res <- fake_plot_result(6L)
  mat <- fake_gs_matrix()
  db <- fake_gs_db()
  ranks <- fake_ranks()
  de <- fake_de_table()
  root <- withr::local_tempdir()

  gmt <- file.path(root, "sets.gmt")
  writeLines("SET_A\tSet A\tG1\tG2\tG3\tG4", gmt)
  counts <- file.path(root, "counts.tsv")
  utils::write.table(
    data.frame(Geneid = c("G1", "G2"), S1 = 1:2, S2 = 3:4),
    counts, sep = "\t", quote = FALSE, row.names = FALSE
  )
  metadata <- file.path(root, "metadata.csv")
  utils::write.table(
    data.frame(Sample_ID = c("S1", "S2"), group = c("WT", "KO")),
    metadata, sep = ",", quote = FALSE, row.names = FALSE
  )

  refs_dir <- file.path(root, "gatom-refs")
  dir.create(refs_dir)
  saveRDS(list(marker = "network"), file.path(refs_dir, "network.kegg.rds"))
  saveRDS(list(marker = "metdb"), file.path(refs_dir, "met.kegg.db.rds"))
  saveRDS(list(marker = "anno"),
          file.path(refs_dir, "org.Hs.eg.gatom.anno.rds"))

  download_dir <- file.path(root, "gatom-download")
  dir.create(download_dir)
  download_files <- c(
    "network.kegg.rds", "met.kegg.db.rds",
    "gene2reaction.kegg.hsa.eg.tsv", "org.Hs.eg.gatom.anno.rds"
  )
  vapply(file.path(download_dir, download_files), function(path) {
    writeBin(as.raw(1L), path)
    TRUE
  }, logical(1L))

  chunk_dir <- file.path(root, "chunks")
  dir.create(chunk_dir)
  chunk_path <- file.path(chunk_dir, "chunk_1_full_objects.qs2")
  file.create(chunk_path)
  expression <- outer(1:12, c(3, -2, 5, 1), `*`)
  chunk_object <- list(
    gseId = "GSE_A", gplId = "GPL_A", E1024 = expression * 1024,
    rownames = 1:12, totalVar = sum(expression^2)
  )
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) list(chunk_object),
    entrez_to_gene = function(entrez, species = "human") {
      stats::setNames(paste0("G", entrez), entrez)
    },
    .package = "bulkiRNA"
  )
  query <- as.integer(1:4)
  queries <- list(query_a = query)
  chunks <- coresh_chunks(chunk_dir, cache = FALSE)
  hit <- tibble::tibble(
    query_name = "query_a", gse = "GSE_A", gpl = "GPL_A", rank = 1L
  )

  tables_dir <- file.path(root, "tables")
  gs_write(res, tables_dir)
  master <- gs_to_master(res)
  gatom_input <- data.frame(
    symbol = c("A", "B", "C"), p = c(0.001, 0.02, 0.4),
    fc = c(2, -1, 0.5), mean = c(100, 80, 60)
  )

  probes <- list(
    bulkirna_api = function() bulkirna_api(quiet = TRUE),
    bulkirna_check_deps = function() bulkirna_check_deps(quiet = TRUE),
    bulkirna_stochastic = function() bulkirna_stochastic(quiet = TRUE),
    coresh_chunks = function() chunks,
    coresh_convergence = function() coresh_convergence(tibble::tibble(
      query_name = c("a", "b"), gse = "GSE_A", gpl = "GPL_A",
      pct_var = c(2, 1), p_value = c(NA_real_, NA_real_),
      log2err = c(NA_real_, NA_real_), size = c(3L, 3L), rank = c(1L, 1L)
    )),
    coresh_loadings = function() coresh_loadings(
      chunk_path, "GSE_A", query, top_n = 5L
    ),
    coresh_match = function() coresh_match(chunk_object, query),
    coresh_sets = function() coresh_sets(
      hit, queries, chunk_dir = chunk_dir, top_n = 5L,
      min_size = 1L, max_size = 10L, jaccard_threshold = 1
    ),
    coresh_validate = function() suppressMessages(coresh_validate(chunk_dir)),
    de_bfc_plot = function() de_bfc_plot(de),
    # de_md_plot()'s first argument is a limma fit, not a DE table, and `coef`
    # has no default. Passing `de_results` skips the topTable() call, and with a
    # character `coef` and an AveExpr column `fit` is never touched.
    de_md_plot = function() de_md_plot(NULL, "treated", de_results = de),
    de_volcano = function() de_volcano(de),
    ensure_dir = function() ensure_dir(file.path(root, "created")),
    filter_confounder_genes = function() filter_confounder_genes(
      c("ACTB", "RPL10", "GENE1")
    ),
    bulki_palettes = function() bulki_palettes(),
    coresh_labels = function() coresh_labels(
      tibble::tibble(gse = "GSE1", gpl = "GPL1", pct_var = 1, query_size = 3L)
    ),
    format_pathway_name = function() format_pathway_name("HALLMARK_MTORC1"),
    gs_palette = function() gs_palette(c("SET_A", "SET_B")),
    gs_plot_size = function() gs_plot_size("bar"),
    gs_scale_fonts = function() gs_scale_fonts(
      gs_plot_bar(plot_res), width = 7, height = 5
    ),
    gatom_de = function() gatom_de(
      gatom_input, id = symbol, pval = p, log2FC = fc, baseMean = mean
    ),
    gatom_download_refs = function() suppressMessages(gatom_download_refs(
      download_dir, species = "human", networks = "kegg"
    )),
    gatom_refs = function() gatom_refs("human", dir = refs_dir),
    gs_filter = function() gs_filter(res, padj = 0.2),
    gs_leading_edge = function() gs_leading_edge(res),
    gs_master_columns = function() gs_master_columns(),
    gs_plot_bar = function() gs_plot_bar(plot_res),
    gs_plot_dot = function() gs_plot_dot(plot_res),
    gs_plot_heatmap = function() gs_plot_heatmap(mat),
    gs_plot_running = function() gs_plot_running(
      gs_test(ranks, db, min_size = 5L, max_size = 50L),
      ranks = ranks, db = db, top_n = 1L
    ),
    gs_ranks = function() gs_ranks(ranks),
    gs_read = function() gs_read(tables_dir),
    gs_save = function() gs_save(
      gs_plot_bar(plot_res), file.path(root, "figure"),
      formats = character(0L), table = FALSE
    ),
    gs_split = function() gs_split(res),
    gs_stat_types = function() gs_stat_types(),
    gs_test = function() gs_test(ranks, db, min_size = 5L, max_size = 50L),
    gs_to_master = function() master,
    gs_top = function() gs_top(res),
    gs_validate_master = function() gs_validate_master(master, error = FALSE),
    gs_write = function() gs_write(res, file.path(root, "tables-2")),
    gsdb_from_file = function() gsdb_from_file(gmt, min_size = 1L),
    gsdb_info = function() gsdb_info(db),
    gsdb_list = function() gsdb_list(),
    gsdb_load = function() gsdb_load(
      "mitopathways", min_size = 1L, max_size = Inf
    ),
    gsdb_register = function() gsdb_register(
      list(SET_A = c("A", "B", "C")), database = "audit"
    ),
    read_counts_matrix = function() read_counts_matrix(counts),
    read_metadata = function() read_metadata(metadata),
    theme_bulki = function() theme_bulki(),
    write_session_provenance = function() suppressMessages(
      write_session_provenance(file.path(root, "provenance.txt"))
    )
  )

  probed <- names(probes)
  expect_identical(
    sort(c(probed, exception_names)), sort(live),
    info = paste(
      "Every live export needs an executable return probe or its own reasoned",
      "exception."
    )
  )
  return_audit_expect_allowlist(
    intersect(probed, exception_names),
    NULL,
    info = "A live export cannot be both probed and excepted."
  )
  return_audit_expect_allowlist(
    character(0L),
    NULL,
    info = "The return exception comparison must stay shape-stable when empty."
  )

  writer_names <- contract$name[contract$category == "writer"]
  for (name in probed) {
    observed <- withVisible(probes[[name]]())
    expected <- contract$expected[match(name, contract$name)]
    return_audit_expect_class(observed$value, expected, info = name)
    if (name %in% writer_names) {
      expect_false(observed$visible, info = paste(name, "must return invisibly"))
    }
  }
})
