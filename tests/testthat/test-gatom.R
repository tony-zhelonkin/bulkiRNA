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

# ---- gatom_de(): baseMean is a rank ----------------------------------------

test_that("gatom_de() takes baseMean on any monotone scale without warning", {
  # gatom only ranks genes by baseMean (gene.keep.top), so a log-scale column
  # is as good as a linear one.
  tt <- fake_tt()
  expect_silent(gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                         baseMean = AveExpr))
  tt$negexpr <- tt$AveExpr - 8
  expect_silent(gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC,
                         baseMean = negexpr))
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
  de <- gatom_de(tt, id = symbol, pval = P.Value, log2FC = logFC, baseMean = 1)
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

# ---- reference files ---------------------------------------------------------

test_that("one file table names every network's files", {
  mm <- bulkiRNA:::.species("Mus musculus")
  files <- bulkiRNA:::.gatom_ref_files
  expect_identical(
    unname(files("kegg", mm)),
    c("network.kegg.rds", "met.kegg.db.rds", "org.Mm.eg.gatom.anno.rds")
  )
  expect_identical(
    unname(files("combined", mm)),
    c("network.combined.rds", "met.combined.db.rds",
      "gene2reaction.combined.mmu.eg.tsv", "org.Mm.eg.gatom.anno.rds")
  )
  expect_identical(
    unname(files("rhea", bulkiRNA:::.species("Homo sapiens"))),
    c("network.rhea.rds", "met.rhea.db.rds", "gene2reaction.rhea.hsa.eg.tsv",
      "org.Hs.eg.gatom.anno.rds")
  )
  # The vignette's prose says met.rhea.lipids.db.rds; the server has this.
  expect_identical(
    unname(files("lipids", mm)),
    c("network.rhea.lipids.rds", "met.lipids.db.rds",
      "gene2reaction.rhea.mmu.eg.tsv", "org.Mm.eg.gatom.anno.rds")
  )
})

test_that("every file in the table exists on the GATOM server", {
  skip_on_cran()
  skip_if_offline("artyomovlab.wustl.edu")
  index <- tryCatch(
    readLines("https://artyomovlab.wustl.edu/publications/supp_materials/GATOM/",
              warn = FALSE),
    error = function(e) NULL
  )
  skip_if(is.null(index), "GATOM server index unreachable")
  served <- regmatches(index, gregexpr('href="[^"/?]+"', index))
  served <- gsub('href="|"', "", unlist(served))
  wanted <- unlist(lapply(c("Mus musculus", "Homo sapiens"), function(sp) {
    lapply(bulkiRNA:::.gatom_networks, bulkiRNA:::.gatom_ref_files,
           species = bulkiRNA:::.species(sp))
  }))
  expect_true(all(wanted %in% served),
              info = paste(setdiff(wanted, served), collapse = ", "))
})

test_that("gatom_refs() validates its arguments", {
  expect_error(gatom_refs("Rattus norvegicus"), "`species` must be one of")
  expect_error(gatom_refs("Homo sapiens", network = "reactome"),
               "`network` must be one of \"kegg\", \"rhea\"", fixed = TRUE)
  expect_error(gatom_refs("Homo sapiens", dir = c("a", "b")), "`dir` must be")
  expect_error(gatom_refs("Homo sapiens", dir = NULL), "`dir` must be")
  expect_error(gatom_refs("Homo sapiens", download = NA), "`download` must be")
})

test_that("gatom_refs() reads one directory and names what is missing there", {
  empty <- withr::local_tempdir()
  err <- tryCatch(gatom_refs("Mus musculus", dir = empty),
                  error = function(e) conditionMessage(e))
  expect_match(err, "org.Mm.eg.gatom.anno.rds", fixed = TRUE)
  expect_match(err, empty, fixed = TRUE)
  expect_match(err, "gatom_download_refs(dir = ", fixed = TRUE)
  # No second search location: the old staged path is never mentioned.
  expect_false(grepl("/opt/gatom-refs", err, fixed = TRUE))
})

write_marker_refs <- function(d, network = "kegg", species = "Homo sapiens") {
  files <- bulkiRNA:::.gatom_ref_files(network, bulkiRNA:::.species(species))
  for (nm in names(files)) {
    path <- file.path(d, files[[nm]])
    if (nm == "gene2reaction_extra") {
      writeLines(c("reaction\tgene", "R00001\t0123"), path)
    } else {
      saveRDS(list(marker = nm), path)
    }
  }
  files
}

test_that("gatom_refs() carries met.db and the gene2reaction table", {
  d <- withr::local_tempdir()
  write_marker_refs(d, "kegg")
  refs <- gatom_refs("Homo sapiens", dir = d)
  expect_s3_class(refs, "gatom_refs")
  expect_identical(refs$met_db$marker, "met_db")   # trap 4: met.db carried
  expect_null(refs$gene2reaction_extra)            # KEGG has none
  expect_identical(refs$checksums, "no SHA256SUMS")
  expect_output(print(refs), "gatom_refs")

  write_marker_refs(d, "combined")
  refs <- gatom_refs("Homo sapiens", dir = d, network = "combined")
  expect_identical(refs$network$marker, "network")
  expect_identical(refs$network_name, "combined")
  # Read as character, so a leading zero survives.
  expect_identical(refs$gene2reaction_extra$gene, "0123")
})

test_that("gatom_refs() verifies files against SHA256SUMS", {
  skip_if_not(exists("sha256sum", asNamespace("tools")), "needs R >= 4.5")
  d <- withr::local_tempdir()
  files <- write_marker_refs(d, "kegg")
  bulkiRNA:::.gatom_write_sums(d)
  refs <- gatom_refs("Homo sapiens", dir = d)
  expect_identical(refs$checksums, "verified")
  expect_identical(unname(refs$sha256),
                   unname(tools::sha256sum(file.path(d, files))))

  saveRDS(list(marker = "tampered"), file.path(d, files[["met_db"]]))
  expect_error(gatom_refs("Homo sapiens", dir = d),
               "do not match its SHA256SUMS")

  writeLines(readLines(file.path(d, "SHA256SUMS"))[-1L],
             file.path(d, "SHA256SUMS"))
  expect_error(gatom_refs("Homo sapiens", dir = d), "does not list")
})

test_that("gatom_refs(download = TRUE) is deprecated and loads from dir", {
  d <- withr::local_tempdir()
  write_marker_refs(d, "kegg")
  called <- NULL
  testthat::local_mocked_bindings(
    gatom_download_refs = function(dir, species, networks, ...) {
      called <<- list(dir = dir, networks = networks)
      invisible(character())
    },
    .package = "bulkiRNA"
  )
  expect_warning(refs <- gatom_refs("Homo sapiens", dir = d, download = TRUE),
                 "`download` is deprecated")
  expect_identical(called, list(dir = d, networks = "kegg"))
  expect_identical(dirname(unname(refs$files)), rep(d, 3L))
})

# ---- argument validation, before gatom is touched --------------------------

fake_refs <- function() {
  structure(list(network = 1, met_db = 2, org_anno = 3,
                 gene2reaction_extra = NULL, species = "Homo sapiens",
                 network_name = "kegg", files = character()),
            class = "gatom_refs")
}

test_that("gatom_module() validates its arguments before touching gatom", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  skip_if_not_installed("igraph")
  de <- data.frame(ID = "IDO1", pval = 0.01, log2FC = 1, baseMean = 100)
  refs <- fake_refs()
  expect_error(gatom_module(de, refs = list()), "`refs` must be")
  expect_error(gatom_module("nope", refs), "`de` must be a data frame")
  expect_error(gatom_module(de[, 1:2], refs), "missing column")
  expect_error(gatom_module(NULL, refs), "both NULL")
  expect_error(gatom_module(de, refs, k_gene = -1), "`k_gene` must be")
  expect_error(gatom_module(de, refs, k_met = 50), "`k_met` was supplied")
  expect_error(gatom_module(de, refs, seed = NA), "`seed` must be")
  expect_error(gatom_module(de, refs, verbose = NA), "`verbose` must be")
  expect_error(gatom_module(de, refs, topology = "reactions"),
               "should be one of")
  expect_error(gatom_module(de, refs, solver = "cplex"),
               "`solver` must be \"rnc\", \"virgo\"", fixed = TRUE)
  expect_error(
    gatom_module(de, refs, gene2reaction_extra = data.frame(gene = "1")),
    "`gene2reaction_extra` must be NULL or a data frame"
  )
})

test_that("gatom_module() keeps its stable formals and appends new ones", {
  fmls <- formals(gatom_module)
  expect_identical(
    names(fmls),
    c("de", "refs", "k_gene", "k_met", "met_de", "seed", "solver", "verbose",
      "gene2reaction_extra", "topology", "keep_reactions_without_enzymes")
  )
  expect_identical(fmls$k_gene, 50)
  expect_identical(fmls$seed, 42)
  expect_identical(fmls$solver, "rnc")
})

test_that("the layers refuse an input from the wrong layer", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  skip_if_not_installed("igraph")
  g <- igraph::make_graph(~ A - B)
  expect_error(gatom_score(g), "`g` must be the result of gatom_graph()",
               fixed = TRUE)
  expect_error(gatom_solve(g), "`gs` must be the result of gatom_score()",
               fixed = TRUE)
  attr(g, "gatom_layer") <- "graph"
  expect_error(gatom_score(g, k_gene = NULL, k_met = NULL), "both NULL")
  expect_error(gatom_solve(g), "`gs` must be the result of gatom_score()",
               fixed = TRUE)
})

test_that("gatom_solver() builds the vignette's solvers and refuses the rest", {
  skip_if_not_installed("mwcsr")
  rnc <- gatom_solver("rnc")
  expect_s3_class(rnc, "mwcs_solver")
  expect_identical(attr(rnc, "gatom_solver"), list(type = "rnc"))
  expect_error(gatom_solver("rnc", threads = 2), "`threads` is a virgo setting")
  expect_error(gatom_solver("annealing"), "should be one of")

  # Never the approximate mode: no CPLEX means an error naming CPLEX_HOME.
  expect_error(gatom_solver("virgo", cplex_dir = ""), "`CPLEX_HOME`")
  missing_dir <- file.path(tempdir(), "no-cplex-here")
  expect_error(gatom_solver("virgo", cplex_dir = missing_dir),
               "`cplex_dir` does not exist")
  empty <- withr::local_tempdir()
  expect_error(gatom_solver("virgo", cplex_dir = empty),
               "not a usable CPLEX installation")
  withr::local_envvar(CPLEX_HOME = "")
  expect_error(gatom_solver("virgo"), "`CPLEX_HOME`")
})

test_that("deprecated solver strings warn and still solve", {
  skip_if_not_installed("mwcsr")
  resolve <- bulkiRNA:::.gatom_resolve_solver
  expect_warning(s <- resolve("annealing"), "deprecated")
  expect_s3_class(s, "mwcs_solver")
  expect_identical(attr(s, "gatom_solver")$type, "annealing")
  expect_warning(resolve("rmwcs"), "deprecated")
  # Any mwcs_solver passes through untouched.
  direct <- mwcsr::rnc_solver(max_iterations = 10L)
  expect_identical(resolve(direct), direct)
})

test_that("GATOM solution weights have a stable numeric representation", {
  solution_weight <- bulkiRNA:::.gatom_solution_weight
  expect_identical(solution_weight(list(weight = 12.5)), 12.5)
  expect_identical(solution_weight(list()), NA_real_)
})

test_that("the exact BUM fit is recovered from scores, and checked", {
  exact <- bulkiRNA:::.gatom_exact_fit
  a <- 0.13; t <- 4e-8
  p <- c(1e-12, 1e-9, 1e-5, 0.2, NA)
  s <- (a - 1) * (log(p) - log(t))
  # gatom prints six decimals: this threshold prints as 0.
  printed <- c(round(a, 6), round(t, 6))
  expect_identical(printed[[2L]], 0)
  expect_equal(exact(p, s, printed, "gene"), c(a, t))
  expect_identical(exact(c(0.1, 0.1), c(1, 1), printed, "gene"), printed)
  expect_error(exact(p, s + c(0, 0, 1, 0, 0), printed, "gene"),
               "no longer follow")
  expect_error(exact(p, s, c(0.5, 0), "gene"), "no longer follow")
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
  expect_error(gatom_module(de, fake_refs()),
               "BiocManager::install\\(\"gatom\"\\)")
  expect_error(gatom_graph(de, fake_refs()),
               "BiocManager::install\\(\"gatom\"\\)")
  expect_error(gatom_save_html(structure(list(), class = "igraph"), "x.html"),
               "BiocManager::install\\(\"gatom\"\\)")
})

test_that("gatom_genes() reads the edge label, not vertices", {
  skip_if_not_installed("igraph")
  g <- igraph::make_graph(~ A - B, B - C, C - D)
  igraph::E(g)$label <- c("IDO1", "KMO", "IDO1")
  igraph::V(g)$label <- rep("WRONG", igraph::vcount(g))
  expect_identical(gatom_genes(g), c("IDO1", "KMO"))

  h <- igraph::make_graph(~ A - B)
  expect_error(gatom_genes(h), "no `label` column")
})

# ---- the real pipeline, on gatom's own example data -----------------------
# gatom ships a 617-gene RefSeq DE table, a KEGG subnetwork, its met.db, a
# mouse annotation and a metabolite DE table. Written to disk under the real
# file names, they exercise every layer offline and in the built package.

gatom_example_refs <- function(env = parent.frame()) {
  skip_if_not_installed("gatom")
  skip_if_not_installed("mwcsr")
  skip_if_not_installed("igraph")
  ex <- new.env()
  utils::data(list = c("networkEx", "met.kegg.dbEx", "org.Mm.eg.gatom.annoEx",
                       "gene.de.rawEx", "met.de.rawEx"),
              package = "gatom", envir = ex)
  d <- withr::local_tempdir(.local_envir = env)
  saveRDS(ex$networkEx, file.path(d, "network.kegg.rds"))
  saveRDS(ex$met.kegg.dbEx, file.path(d, "met.kegg.db.rds"))
  saveRDS(ex$org.Mm.eg.gatom.annoEx, file.path(d, "org.Mm.eg.gatom.anno.rds"))
  refs <- gatom_refs("Mus musculus", dir = d)
  de <- gatom_de(ex$gene.de.rawEx, id = ID, pval = pval, log2FC = log2FC,
                 baseMean = baseMean)
  list(refs = refs, de = de, met_de = ex$met.de.rawEx)
}

quiet_gatom <- function(expr) {
  # BioNet warns when a BUM parameter sits on its bound; that is gatom's own
  # behaviour on the example data and not what these tests examine.
  suppressMessages(withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("limit of the defined parameter space", conditionMessage(w),
                fixed = TRUE)) invokeRestart("muffleWarning")
    }
  ))
}

test_that("gatom_graph() reports the gene.keep.top cut-off", {
  ex <- gatom_example_refs()
  expect_message(g <- gatom_graph(ex$de, ex$refs),
                 "617 genes in `de`, 617 kept by gatom's gene.keep.top (12000)",
                 fixed = TRUE)
  expect_identical(attr(g, "gatom_layer"), "graph")
  expect_identical(attr(g, "id_type"), "RefSeq")
  expect_identical(attr(g, "network"), "kegg")
  expect_identical(attr(g, "gene_keep_top"), 12000)
  expect_setequal(attr(g, "graph_genes"),
                  unique(as.character(igraph::E(g)$gene)))
  expect_identical(attr(g, "graph_edges"), igraph::ecount(g))

  # The cut-off is gatom's: rank by baseMean, keep the top 12,000.
  big <- ex$de[rep(seq_len(nrow(ex$de)), 25L), ]
  big$ID <- paste0(big$ID, "_", seq_len(nrow(big)))
  big$baseMean <- seq_len(nrow(big))
  big$pval <- seq(1e-8, 0.9, length.out = nrow(big))
  g_big <- quiet_gatom(gatom_graph(rbind(ex$de, big), ex$refs))
  expect_identical(attr(g_big, "genes_in_de"), nrow(ex$de) + nrow(big))
  expect_identical(attr(g_big, "genes_kept"), 12000L)
})

test_that("gatom_graph() stops on a graph with no edges, naming the id type", {
  ex <- gatom_example_refs()
  testthat::local_mocked_bindings(
    makeMetabolicGraph = function(...) stop("No edges in the graph!"),
    .package = "gatom"
  )
  expect_error(quiet_gatom(gatom_graph(ex$de, ex$refs)),
               "no edges: gatom read `de$ID` as RefSeq ids", fixed = TRUE)
})

test_that("gatom_graph() is deterministic and leaves the caller's stream", {
  # gatom samples 1,000 rows to detect the id type of a longer table.
  ex <- gatom_example_refs()
  long <- ex$de[rep(seq_len(nrow(ex$de)), 3L), ]
  long$ID <- c(ex$de$ID, paste0(ex$de$ID, "_a"), paste0(ex$de$ID, "_b"))
  set.seed(31L)
  control <- stats::runif(1L)
  set.seed(31L)
  g1 <- quiet_gatom(gatom_graph(long, ex$refs))
  expect_identical(stats::runif(1L), control)
  g2 <- quiet_gatom(gatom_graph(long, ex$refs))
  expect_identical(attr(g1, "id_type"), "RefSeq")
  expect_identical(igraph::ecount(g1), igraph::ecount(g2))
})

test_that("the layers compose to gatom_module() and record the run", {
  ex <- gatom_example_refs()
  g <- quiet_gatom(gatom_graph(ex$de, ex$refs))
  gs <- quiet_gatom(gatom_score(g, k_gene = 25))
  m <- gatom_solve(gs, "rnc")

  expect_identical(attr(gs, "gatom_layer"), "scored")
  expect_identical(attr(m, "gatom_layer"), "module")
  expect_lt(attr(gs, "gene_threshold"), 1)
  expect_gt(attr(gs, "gene_bum_alpha"), 0)
  expect_lte(attr(gs, "gene_fdr"), 0.1)
  expect_true(is.na(attr(gs, "met_fdr")))
  for (nm in c("network", "genes_kept", "graph_genes", "k_gene",
               "gene_threshold")) {
    expect_identical(attr(m, nm), attr(gs, nm), info = nm)
  }
  expect_identical(attr(m, "solver"), "rnc")
  expect_identical(attr(m, "seed"), 42)
  expect_false(attr(m, "solved_to_optimality"))   # rnc is a heuristic
  expect_identical(attr(m, "n_edges"), igraph::ecount(m))
  expect_gt(igraph::ecount(m), 0)
  expect_lt(igraph::ecount(m), igraph::ecount(g))

  one_call <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25))
  expect_identical(gatom_genes(one_call), gatom_genes(m))
  expect_identical(attr(one_call, "solution_weight"),
                   attr(m, "solution_weight"))
})

test_that("larger k_gene gives a larger threshold and no smaller module", {
  ex <- gatom_example_refs()
  g <- quiet_gatom(gatom_graph(ex$de, ex$refs))
  gs <- lapply(c(10, 25), function(k) quiet_gatom(gatom_score(g, k_gene = k)))
  expect_lt(attr(gs[[1L]], "gene_threshold"), attr(gs[[2L]], "gene_threshold"))
  m <- lapply(gs, gatom_solve)
  expect_lte(igraph::ecount(m[[1L]]), igraph::ecount(m[[2L]]))
})

test_that("metabolite data is scored when given, alone or with genes", {
  ex <- gatom_example_refs()
  both <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25,
                                   met_de = ex$met_de, k_met = 25))
  expect_false(is.na(attr(both, "met_threshold")))
  expect_false(is.na(attr(both, "gene_threshold")))

  met_only <- quiet_gatom(gatom_module(NULL, ex$refs, k_gene = NULL,
                                       met_de = ex$met_de, k_met = 25))
  expect_true(is.na(attr(met_only, "gene_threshold")))
  expect_identical(attr(met_only, "genes_in_de"), 0L)
  expect_gt(igraph::vcount(met_only), 0)

  g <- quiet_gatom(gatom_graph(ex$de, ex$refs))
  expect_error(gatom_score(g, k_gene = 25, k_met = 25),
               "graph carries no metabolite p-values")
})

test_that("the metabolite topology runs and still reads genes off edges", {
  ex <- gatom_example_refs()
  m <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25,
                                topology = "metabolites"))
  expect_identical(attr(m, "topology"), "metabolites")
  expect_gt(length(gatom_genes(m)), 0)
})

test_that("a p-value distribution BUM cannot fit is an error, not a warning", {
  ex <- gatom_example_refs()
  de <- ex$de
  set.seed(1)
  de$pval <- stats::runif(nrow(de))
  expect_error(quiet_gatom(gatom_module(de, ex$refs, k_gene = 25)),
               "do not fit a beta-uniform mixture")
})

test_that("scoring and solving are seed-stable and leave the caller's stream", {
  ex <- gatom_example_refs()
  g <- quiet_gatom(gatom_graph(ex$de, ex$refs))
  set.seed(606L)
  control <- stats::runif(1L)
  set.seed(606L)
  a <- gatom_solve(quiet_gatom(gatom_score(g, k_gene = 25, seed = 7)),
                   seed = 7)
  b <- gatom_solve(quiet_gatom(gatom_score(g, k_gene = 25, seed = 7)),
                   seed = 7)
  expect_identical(stats::runif(1L), control)
  expect_identical(gatom_genes(a), gatom_genes(b))
  expect_identical(attr(a, "solution_weight"), attr(b, "solution_weight"))
})

test_that("gatom_module(gene2reaction_extra =) is deprecated and replaces refs'", {
  ex <- gatom_example_refs()
  seen <- NULL
  real <- gatom::makeMetabolicGraph
  testthat::local_mocked_bindings(
    makeMetabolicGraph = function(..., gene2reaction.extra = NULL) {
      seen <<- gene2reaction.extra
      real(..., gene2reaction.extra = gene2reaction.extra)
    },
    .package = "gatom"
  )
  extra <- data.frame(gene = "66925", reaction = "R00001")
  expect_warning(
    quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25,
                             gene2reaction_extra = extra)),
    "`gene2reaction_extra` is deprecated"
  )
  expect_identical(seen, extra)
  quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25))
  expect_null(seen)
})

test_that("gatom_pathways() is the vignette's fora call plus collapse", {
  ex <- gatom_example_refs()
  m <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25))
  universe <- attr(m, "graph_genes")
  module_genes <- unique(igraph::E(m)$gene)
  # The example annotation carries no pathways, so give it a planted one
  # (the module itself) among random ones drawn from the universe.
  set.seed(2)
  pathways <- c(list(planted = module_genes),
                stats::setNames(lapply(1:20, function(i) sample(universe, 8)),
                                paste0("random", 1:20)))
  refs <- ex$refs
  refs$org_anno$pathways <- pathways

  p <- gatom_pathways(m, refs)
  direct <- fgsea::fora(pathways = pathways, genes = module_genes,
                        universe = universe, minSize = 5)
  expect_identical(p$pathway, direct$pathway)
  expect_identical(p$pval, direct$pval)
  expect_identical(p$main[p$pathway == "planted"], TRUE)
  expect_false(any(p$main[p$padj >= 0.05]))
  expect_identical(attr(p, "n_universe"), length(universe))

  expect_true(all(is.na(gatom_pathways(m, refs, collapse = FALSE)$main)))
  expect_error(gatom_pathways(m, ex$refs), "carries no pathways")
  expect_error(gatom_pathways(m, refs, universe = character()),
               "`universe` must be")
})

test_that("gatom_save_html() writes a self-contained file and makes its dir", {
  ex <- gatom_example_refs()
  m <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25))
  pandoc_before <- Sys.getenv("RSTUDIO_PANDOC", unset = NA)
  out <- file.path(withr::local_tempdir(), "nested", "module.html")
  expect_invisible(gatom_save_html(m, out, name = "Example"))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 1000)
  # The pandoc location is scoped to the call.
  expect_identical(Sys.getenv("RSTUDIO_PANDOC", unset = NA), pandoc_before)
})

test_that("gatom_save_pdf() writes the vignette's seeded layout", {
  ex <- gatom_example_refs()
  m <- quiet_gatom(gatom_module(ex$de, ex$refs, k_gene = 25))
  out <- file.path(withr::local_tempdir(), "nested", "module.pdf")
  set.seed(7)
  before <- .Random.seed
  # ggplot2 drops the unlabelled points of gatom's own layout with a warning.
  expect_invisible(suppressWarnings(gatom_save_pdf(m, out, name = "Example")))
  expect_identical(.Random.seed, before)
  expect_identical(readBin(out, "raw", 4L), charToRaw("%PDF"))

  expect_error(gatom_save_pdf(m, out, n_iter = 0), "`n_iter` must be")
  expect_error(gatom_save_pdf(m, out, force = -1), "`force` must be")
  expect_error(gatom_save_pdf(data.frame(a = 1), out),
               "must be an igraph module")
})

test_that("a failed gatom_save_pdf() closes gatom's device and its file", {
  skip_if_not_installed("gatom")
  skip_if_not_installed("igraph")
  # gatom opens the PDF device first and fails while drawing.
  local_mocked_bindings(
    saveModuleToPdf = function(module, file, name, n_iter, force) {
      grDevices::pdf(file)
      stop("missing value where TRUE/FALSE needed")
    },
    .package = "gatom"
  )
  m <- igraph::make_graph(~ a - b)
  out <- file.path(withr::local_tempdir(), "module.pdf")
  devices <- grDevices::dev.list()
  expect_error(gatom_save_pdf(m, out), "could not draw `m`: missing value")
  expect_identical(grDevices::dev.list(), devices)
  expect_false(file.exists(out))
})
