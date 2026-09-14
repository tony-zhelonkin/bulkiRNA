# gatom_* module. The gatom_de() validation tests encode the traps and run
# everywhere; the pipeline tests need the gatom stack and skip without it.

test_that("pinned seeding preserves draws under R's default generator", {
  local_pinned_rng()

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(42)
  expected <- rnorm(3)

  # Protects the two GATOM resets from moving results for default-RNG callers.
  actual <- bulkiRNA:::.with_pinned_seed(42, rnorm(3))
  expect_identical(actual, expected)
})

test_that("a NULL pinned seed evaluates once without touching RNG state", {
  local_pinned_rng()

  set.seed(818L)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  evaluations <- 0L

  value <- bulkiRNA:::.with_pinned_seed(NULL, {
    evaluations <- evaluations + 1L
    "evaluated"
  })

  expect_identical(value, "evaluated")
  expect_identical(evaluations, 1L)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

test_that("pinned seeding restores RNG state when evaluation fails", {
  local_pinned_rng()

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(919L)
  kind_before <- RNGkind()
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  expect_error(
    bulkiRNA:::.with_pinned_seed(42L, {
      runif(1L)
      stop("seeded failure", call. = FALSE)
    }),
    "seeded failure"
  )

  expect_identical(RNGkind(), kind_before)
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

fake_tt <- function(n = 6) {
  data.frame(
    symbol  = c("IDO1", "KMO", "KYNU", "HAAO", "QPRT", "TDO2")[seq_len(n)],
    P.Value = seq(1e-6, 0.05, length.out = n),
    adj.P.Val = seq(1e-4, 0.2, length.out = n),
    logFC   = seq(-2, 2, length.out = n),
    AveExpr = seq(4, 10, length.out = n),
    stringsAsFactors = FALSE
  )
}

# ---- gatom_de(): trap 1, raw p-values -------------------------------------

test_that("gatom_de() builds the four gatom columns, sorted and deduplicated", {
  tt <- fake_tt()
  tt <- rbind(tt, tt[1, ])
  tt$P.Value[nrow(tt)] <- 0.5     # duplicate IDO1 with a worse p-value
  de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                 baseMean = 2^AveExpr)
  expect_s3_class(de, "gatom_de")
  expect_identical(names(de), c("ID", "pval", "log2FC", "baseMean"))
  expect_false(any(duplicated(de$ID)))
  expect_identical(de$pval, sort(de$pval))
  expect_equal(de$pval[de$ID == "IDO1"], 1e-6)  # lowest p-value wins
})

test_that("gatom_de() rejects an adjusted p-value column by name", {
  tt <- fake_tt()
  expect_error(
    gatom_de(tt, id = symbol, pval = adj.P.Val, log2FC = logFC,
             baseMean = 2^AveExpr),
    "adjusted p-value"
  )
  tt$FDR <- tt$adj.P.Val
  expect_error(
    gatom_de(tt, id = symbol, pval = FDR, log2FC = logFC, baseMean = 2^AveExpr),
    "adjusted p-value"
  )
  tt$padj <- tt$adj.P.Val
  expect_error(
    gatom_de(tt, id = symbol, pval = padj, log2FC = logFC,
             baseMean = 2^AveExpr),
    "adjusted p-value"
  )
  tt$qval <- tt$adj.P.Val
  expect_error(
    gatom_de(tt, id = symbol, pval = qval, log2FC = logFC,
             baseMean = 2^AveExpr),
    "adjusted p-value"
  )
  tt$q.value <- tt$adj.P.Val
  expect_error(
    gatom_de(tt, id = symbol, pval = q.value, log2FC = logFC,
             baseMean = 2^AveExpr),
    "adjusted p-value"
  )
})

test_that("gatom_de() rejects p-values outside [0, 1] and all-NA p-values", {
  tt <- fake_tt()
  tt$score <- -log10(tt$P.Value)
  expect_error(
    gatom_de(tt, id = symbol, pval = score, log2FC = logFC,
             baseMean = 2^AveExpr),
    "\\[0, 1\\]"
  )
  tt$neg <- -tt$P.Value
  expect_error(
    gatom_de(tt, id = symbol, pval = neg, log2FC = logFC,
             baseMean = 2^AveExpr),
    "\\[0, 1\\]"
  )
  tt$allna <- NA_real_
  expect_error(
    gatom_de(tt, id = symbol, pval = allna, log2FC = logFC,
             baseMean = 2^AveExpr),
    "entirely NA"
  )
})

test_that("gatom_de() requires a numeric p-value column", {
  tt <- fake_tt()
  tt$chr <- as.character(tt$P.Value)
  expect_error(
    gatom_de(tt, id = symbol, pval = chr, log2FC = logFC,
             baseMean = 2^AveExpr),
    "`pval` must be numeric"
  )
})

# ---- gatom_de(): trap 2, linear baseMean ----------------------------------

test_that("gatom_de() warns when baseMean looks log-scaled", {
  tt <- fake_tt()
  expect_warning(
    gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
             baseMean = AveExpr),          # log scale, max 10
    "LINEAR"
  )
  tt$negexpr <- tt$AveExpr - 8
  expect_warning(
    gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
             baseMean = negexpr),
    "negative"
  )
  expect_silent(
    gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
             baseMean = 2^AveExpr)
  )
})

test_that("gatom_de() errors on an all-NA baseMean", {
  tt <- fake_tt()
  expect_error(
    gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
             baseMean = NA_real_),
    "entirely NA"
  )
})

# ---- gatom_de(): general contract -----------------------------------------

test_that("gatom_de() accepts a constant baseMean placeholder", {
  tt <- fake_tt()
  de <- suppressWarnings(
    gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC, baseMean = 1)
  )
  expect_true(all(de$baseMean == 1))
  expect_equal(nrow(de), nrow(tt))
})

test_that("gatom_de() drops rows with missing id or pval", {
  tt <- fake_tt()
  tt$symbol[1] <- NA
  tt$P.Value[2] <- NA
  expect_message(
    de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                   baseMean = 2^AveExpr),
    "dropped"
  )
  expect_equal(nrow(de), nrow(tt) - 2L)
})

test_that("gatom_de() validates its inputs", {
  expect_error(gatom_de(list(a = 1), id = a, pval = a, log2FC = a,
                        baseMean = a),
               "must be a data frame")
  tt <- fake_tt()[0, ]
  expect_error(gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                        baseMean = 1),
               "zero rows")
  tt <- fake_tt()
  expect_error(gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                        baseMean = rep(1, 3)),
               "length 1 or nrow")
  expect_error(gatom_de(tt, id = symbol, pval = nope, log2FC = logFC,
                        baseMean = 1),
               "could not be evaluated")
})

# ---- gatom_refs(): trap 4 and species handling ----------------------------

test_that("gatom_refs() validates species", {
  expect_error(gatom_refs("Rattus norvegicus"), "`species` must be one of")
  expect_error(
    gatom_refs("Homo sapiens", network = "rhea"),
    "`network` must be one of \"kegg\", \"combined\"",
    fixed = TRUE
  )
  expect_error(gatom_refs("Homo sapiens", dir = c("a", "b")), "`dir` must be")
  expect_error(gatom_refs("Homo sapiens", download = NA), "`download` must be")
})

test_that("gatom_refs() names the missing file and how to fetch it", {
  empty <- tempfile("gatomrefs"); dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  testthat::local_mocked_bindings(
    .gatom_search_dirs = function(dir, download, default_dir) dir,
    .package = "bulkiRNA"
  )
  err <- tryCatch(gatom_refs("Mus musculus", dir = empty),
                  error = function(e) conditionMessage(e))
  expect_match(err, "org.Mm.eg.gatom.anno.rds")
  expect_match(err, "gatom_download_refs\\(species = \"Mus_musculus\"")
  expect_match(err, empty, fixed = TRUE)
})

test_that("gatom_refs() finds files in an explicit dir", {
  d <- tempfile("gatomrefs"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  saveRDS(list(marker = "network"), file.path(d, "network.kegg.rds"))
  saveRDS(list(marker = "metdb"), file.path(d, "met.kegg.db.rds"))
  saveRDS(list(marker = "anno"), file.path(d, "org.Hs.eg.gatom.anno.rds"))
  refs <- gatom_refs("Homo sapiens", dir = d)
  expect_s3_class(refs, "gatom_refs")
  expect_identical(refs$species, "Homo sapiens")
  expect_identical(refs$met_db$marker, "metdb")     # trap 4: met.db carried
  expect_identical(refs$org_anno$marker, "anno")
  expect_output(print(refs), "gatom_refs")
})

test_that("gatom_refs() resolves the combined network filenames", {
  d <- tempfile("gatomrefs"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  saveRDS(list(marker = "combined-network"),
          file.path(d, "network.combined.rds"))
  saveRDS(list(marker = "combined-metdb"),
          file.path(d, "met.combined.db.rds"))
  saveRDS(list(marker = "anno"), file.path(d, "org.Hs.eg.gatom.anno.rds"))

  refs <- gatom_refs("Homo sapiens", dir = d, network = "combined")

  expect_identical(refs$network$marker, "combined-network")
  expect_identical(refs$met_db$marker, "combined-metdb")
  expect_identical(
    unname(basename(refs$files)),
    c("network.combined.rds", "met.combined.db.rds",
      "org.Hs.eg.gatom.anno.rds")
  )
})

# ---- guards and argument validation on the pipeline entry points ----------

test_that("gatom_module() validates its arguments before touching gatom", {
  fake_refs <- structure(list(network = 1, met_db = 2, org_anno = 3,
                              species = "Homo sapiens", files = character()),
                         class = "gatom_refs")
  de <- data.frame(ID = "IDO1", pval = 0.01, log2FC = 1, baseMean = 100)
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  expect_error(gatom_module(de, refs = list()), "`refs` must be")
  expect_error(gatom_module("nope", fake_refs), "`de` must be a data frame")
  expect_error(gatom_module(de[, 1:2], fake_refs), "missing column")
  expect_error(gatom_module(de, fake_refs, k_gene = -1), "`k_gene` must be")
  expect_error(gatom_module(de, fake_refs, k_met = 50), "`k_met` was supplied")
  err <- tryCatch(
    gatom_module(de, fake_refs,
                 gene2reaction_extra = data.frame(gene = "1")),
    error = function(e) conditionMessage(e)
  )
  expect_match(
    err,
    "`gene2reaction_extra` must be NULL or a data frame with `gene` and `reaction` columns",
    fixed = TRUE
  )
  expect_match(err, "data.table::fread()", fixed = TRUE)
  expect_error(gatom_module(de, fake_refs, seed = NA), "`seed` must be")
  expect_error(gatom_module(de, fake_refs, solver = "cplex"),
               "`solver` must be one of")
})

test_that("gatom_module() passes extra mappings and preserves NULL omission", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  skip_if_not_installed("igraph")
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))

  seen <- list()
  graph <- igraph::make_graph(~ A - B)
  testthat::local_mocked_bindings(
    makeMetabolicGraph = function(..., gene2reaction.extra = NULL) {
      seen[[length(seen) + 1L]] <<- list(value = gene2reaction.extra)
      graph
    },
    scoreGraph = function(g, ...) g,
    .package = "gatom"
  )
  testthat::local_mocked_bindings(
    rnc_solver = function(...) NULL,
    solve_mwcsp = function(solver, instance) {
      list(graph = instance, weight = 1)
    },
    .package = "mwcsr"
  )

  fake_refs <- structure(
    list(network = 1, met_db = 2, org_anno = 3,
         species = "Homo sapiens", files = character()),
    class = "gatom_refs"
  )
  de <- data.frame(ID = "IDO1", pval = 0.01, log2FC = 1, baseMean = 100)
  extra <- data.frame(gene = "3620", reaction = "R01920")

  # The separate scoring and solver resets must not replace the caller's draw.
  set.seed(606L)
  control <- runif(1L)
  set.seed(606L)
  with_extra <- gatom_module(de, fake_refs, gene2reaction_extra = extra)
  omitted <- gatom_module(de, fake_refs)
  explicit_null <- gatom_module(de, fake_refs, gene2reaction_extra = NULL)
  after <- runif(1L)

  expect_s3_class(with_extra, "igraph")
  expect_identical(seen[[1L]]$value, extra)
  expect_null(seen[[2L]]$value)
  expect_null(seen[[3L]]$value)
  expect_identical(omitted, explicit_null)
  expect_identical(after, control)
})

test_that("gatom_module() appends its optional mapping argument", {
  fmls <- formals(gatom_module)
  expect_identical(
    names(fmls)[seq_len(8L)],
    c("de", "refs", "k_gene", "k_met", "met_de", "seed", "solver", "verbose")
  )
  expect_identical(names(fmls)[[9L]], "gene2reaction_extra")
  expect_null(fmls$gene2reaction_extra)
})

test_that("GATOM solution weights have a stable numeric representation", {
  solution_weight <- bulkiRNA:::.gatom_solution_weight
  expect_identical(solution_weight(list(weight = 12.5)), 12.5)
  expect_identical(solution_weight(list()), NA_real_)
})

test_that("gatom_genes() rejects non-igraph input", {
  skip_if_not_installed("igraph")
  expect_error(gatom_genes(data.frame(a = 1)), "must be an igraph module")
})

test_that("gatom_save_html() rejects non-igraph input", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("igraph")
  expect_error(gatom_save_html(data.frame(a = 1), "x.html"),
               "must be an igraph module")
})

test_that("entry points guard Suggests with an actionable message", {
  skip_if(requireNamespace("gatom", quietly = TRUE),
          "gatom is installed; the guard cannot fire")
  de <- data.frame(ID = "IDO1", pval = 0.01, log2FC = 1, baseMean = 100)
  refs <- structure(list(species = "Homo sapiens"), class = "gatom_refs")
  expect_error(gatom_module(de, refs), "BiocManager::install\\(\"gatom\"\\)")
  expect_error(gatom_save_html(structure(list(), class = "igraph"), "x.html"),
               "BiocManager::install\\(\"gatom\"\\)")
})

test_that("gatom_genes() reads Symbol from EDGES, not vertices (trap 3)", {
  skip_if_not_installed("igraph")
  g <- igraph::make_graph(~ A - B, B - C, C - D)
  igraph::E(g)$Symbol <- c("IDO1", "KMO", "IDO1")
  igraph::V(g)$Symbol <- rep("WRONG", igraph::vcount(g))
  expect_identical(gatom_genes(g), c("IDO1", "KMO"))

  h <- igraph::make_graph(~ A - B)
  expect_error(gatom_genes(h), "no `Symbol` column")
})

# ---- the real pipeline ----------------------------------------------------

gatom_stack_ready <- function() {
  requireNamespace("gatom", quietly = TRUE) &&
    requireNamespace("mwcsr", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE) &&
    dir.exists("/opt/gatom-refs")
}

test_that("the full GATOM pipeline runs and is seed-stable", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  skip_if_not_installed("igraph")
  skip_if_not(dir.exists("/opt/gatom-refs"), "staged GATOM references absent")

  refs <- gatom_refs("Homo sapiens", dir = "/opt/gatom-refs")
  set.seed(1)
  genes <- c("IDO1", "KMO", "KYNU", "HAAO", "QPRT", "TDO2", "AFMID", "NADSYN1",
             "NMNAT1", "ACMSD", "IDO2", "TPH1", "DDC", "MAOA", "ALDH2",
             "KYAT1", "KYAT3", "NAMPT", "NAPRT", "NMNAT3")
  tt <- data.frame(
    symbol  = genes,
    P.Value = 10^-seq(6, 2, length.out = length(genes)),
    logFC   = seq(2, -2, length.out = length(genes)),
    AveExpr = seq(9, 5, length.out = length(genes)),
    stringsAsFactors = FALSE
  )
  de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                 baseMean = 2^AveExpr)
  expect_equal(nrow(de), length(genes))

  m <- gatom_module(de, refs, k_gene = 50, seed = 42)
  expect_s3_class(m, "igraph")
  expect_gt(igraph::vcount(m), 0)
  expect_identical(attr(m, "k_gene"), 50)
  expect_identical(attr(m, "seed"), 42)
  expect_identical(attr(m, "solver"), "rnc")
  expect_type(attr(m, "solution_weight"), "double")
  expect_length(attr(m, "solution_weight"), 1L)
  expect_identical(attr(m, "n_nodes"), igraph::vcount(m))
  expect_identical(attr(m, "n_edges"), igraph::ecount(m))

  # trap 3 -- genes come off the edges, and vertices carry no Symbol
  mg <- gatom_genes(m)
  expect_type(mg, "character")
  expect_gt(length(mg), 0)
  expect_true(all(mg %in% igraph::as_data_frame(m, "edges")$Symbol))

  # seed stability: same seed, same module size
  m2 <- gatom_module(de, refs, k_gene = 50, seed = 42)
  expect_identical(igraph::vcount(m2), igraph::vcount(m))
  expect_identical(igraph::ecount(m2), igraph::ecount(m))
  expect_identical(gatom_genes(m2), mg)
})

test_that("a realistic p-value distribution gives a non-degenerate BUM fit", {
  # The test above deliberately feeds an all-significant toy list: it proves
  # the plumbing and seed stability, but BUM cannot fit a mixture with no
  # null component (threshold 1.0, "parameters on the limit"). This test
  # feeds what the real pipeline feeds: the whole enzyme universe, mostly
  # null, with a planted kynurenine-pathway signal. Runs in ~17 s.
  skip_if_not(gatom_stack_ready(), "gatom stack or staged references absent")

  refs <- gatom_refs("Homo sapiens", dir = "/opt/gatom-refs")
  anno <- refs$org_anno
  enz <- unique(anno$genes$symbol[anno$genes$gene %in%
                                    unique(anno$gene2enzyme$gene)])
  enz <- enz[!is.na(enz) & nzchar(enz)]
  skip_if(length(enz) < 5000, "annotation has too few enzyme symbols")

  set.seed(20260810)
  planted <- intersect(
    c("IDO1", "KMO", "KYNU", "HAAO", "QPRT", "TDO2", "AFMID", "ACMSD",
      "NADSYN1", "NMNAT1", "KYAT1", "KYAT3", "NAMPT", "NAPRT"), enz)
  nulls <- setdiff(enz, planted)
  tt <- data.frame(
    symbol  = c(planted, nulls),
    P.Value = c(10^-runif(length(planted), 6, 9),
                rbeta(length(nulls), 0.7, 1)),   # null bulk, graded
    logFC   = c(rnorm(length(planted), 2, 0.3),
                rnorm(length(nulls), 0, 0.5)),
    AveExpr = runif(length(planted) + length(nulls), 4, 12),
    stringsAsFactors = FALSE
  )
  de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                 baseMean = 2^AveExpr)

  warnings_seen <- character()
  msgs <- testthat::capture_messages(
    m <- withCallingHandlers(
      gatom_module(de, refs, k_gene = 50, seed = 42),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )

  # BUM actually fitted: no boundary warning, no zeroed edge scores
  expect_false(any(grepl("limit of the defined parameter space",
                         warnings_seen)))
  expect_false(any(grepl("inappropriate p-value distribution",
                         warnings_seen)))
  # ... and the reported gene p-value threshold is well below 1
  thr_line <- grep("Gene p-value threshold", msgs, value = TRUE)
  if (length(thr_line)) {
    thr <- as.numeric(sub(".*threshold:\\s*", "", thr_line[[1]]))
    expect_lt(thr, 0.5)
  }

  # a real selection, not everything and not nothing
  expect_gt(igraph::vcount(m), 10)
  expect_lt(igraph::vcount(m), attr(m, "graph_nodes"))
  genes <- gatom_genes(m)
  expect_gt(length(genes), 5)
  expect_true(any(planted %in% genes))   # the planted signal is recovered

  # still seed-stable at realistic scale
  m2 <- gatom_module(de, refs, k_gene = 50, seed = 42)
  expect_identical(igraph::vcount(m2), igraph::vcount(m))
  expect_identical(gatom_genes(m2), genes)
})

test_that("gatom_save_html() writes a self-contained file and makes its dir", {
  skip_if_not(gatom_stack_ready(), "gatom stack or staged references absent")
  refs <- gatom_refs("Homo sapiens", dir = "/opt/gatom-refs")
  tt <- data.frame(
    symbol  = c("IDO1", "KMO", "KYNU", "HAAO", "QPRT", "TDO2"),
    P.Value = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3),
    logFC   = c(2, 1.5, 1, -1, -1.5, -2),
    AveExpr = c(9, 8.5, 8, 7.5, 7, 6.5),
    stringsAsFactors = FALSE
  )
  de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                 baseMean = 2^AveExpr)
  m <- gatom_module(de, refs, k_gene = 50, seed = 42)

  tmp <- tempfile("gatomhtml"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  out <- file.path(tmp, "nested", "module.html")
  expect_invisible(gatom_save_html(m, out, name = "Kynurenine"))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 1000)
})

# --- regression from the IO review -------------------------------------------

test_that("a fresh download is searched before the staged container copy", {
  # `dir = NULL, download = TRUE` fetched references into 00_data/references
  # /gatom and then loaded the older /opt/gatom-refs ones, so the download had
  # no visible effect.
  f <- bulkiRNA:::.gatom_search_dirs
  default_dir <- "00_data/references/gatom"

  expect_identical(
    f(NULL, TRUE, default_dir),
    c(default_dir, "/opt/gatom-refs")
  )
  expect_identical(
    f("/my/refs", TRUE, default_dir),
    c("/my/refs", "/opt/gatom-refs", default_dir)
  )
  # Without a download the staged copy still wins over the default dir.
  expect_identical(
    f(NULL, FALSE, default_dir),
    c("/opt/gatom-refs", default_dir)
  )
})
