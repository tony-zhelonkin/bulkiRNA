# Golden-output capture for the bulkiRNA refactor.
#
# Runs every externally-called toolkit function against tests/fixtures/ and
# serializes a normalized form of each result. Step C of the refactor re-runs
# verify_golden.R against these to prove behaviour did not drift.
#
# Compute results are stored as-is; plots are stored as ggplot_build()$data —
# the computed layer data, not pixels, so a silent numerical change is caught
# without being fragile to renderer/patch-version churn.
#
# Run:  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
#         -v <cache>:/cache -v <repo>:/pkg -w /pkg \
#         scdock-r-dev:v0.5.11 Rscript tests/golden/capture_golden.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(limma)
})

# --capture (default) writes the goldens; --verify recomputes and diffs.
MODE <- if (any(grepl("--verify", commandArgs(TRUE)))) "verify" else "capture"

# --cases=a,b,c restricts a capture to named cases. Step C re-captures exactly
# the handful of goldens whose change is sanctioned and written down, and a
# blanket re-capture would silently bless every other drift at the same time.
# With a filter set, the manifest is left alone: it describes the full run.
ONLY <- {
  a <- grep("^--cases=", commandArgs(TRUE), value = TRUE)
  if (length(a)) strsplit(sub("^--cases=", "", a[1]), ",")[[1]] else NULL
}

FIX  <- "tests/fixtures"
OUT  <- "tests/golden/data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- load the toolkit -------------------------------------------------------
# Step C: this used to source() every file under scripts/. It now loads the
# installed package, so the goldens captured from the old script library are
# compared against the *new* implementations reached through the deprecation
# shims in R/deprecated-gs.R and R/deprecated-plot.R. Every case below still
# calls a frozen name, unchanged -- that is the point. Until this switch, a
# green run only proved the old code still worked.
# export_all = TRUE since 1.0.0: the legacy names left the public API and are
# now non-exported fixtures, kept precisely so these cases stay runnable.
suppressWarnings(suppressMessages(
  pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
))
# The shims all emit .Deprecated() by design; the harness calls them on purpose.
options(warn = -1)

# Harness self-test: perturb one function so --verify MUST report a failure.
# Proves the gate is not vacuous. Never set in normal use.
if (nzchar(Sys.getenv("GOLDEN_SELFTEST_PERTURB"))) {
  .orig_fpn <- format_pathway_name
  format_pathway_name <- function(...) toupper(.orig_fpn(...))
  message("SELFTEST: format_pathway_name perturbed")
}

# ---- fixture ----------------------------------------------------------------
rd      <- function(x) readRDS(file.path(FIX, x))
tt      <- rd("toptable.rds")
counts  <- rd("counts.rds")
meta    <- rd("metadata.rds")
gsets   <- rd("genesets.rds")
de_real <- rd("de_real_symbols.rds")

# ---- normalizers ------------------------------------------------------------
# Reduce any result to a stable, comparable structure.
normalize <- function(x) {
  if (inherits(x, "patchwork")) {
    # a patchwork is a ggplot with attached panels; build each
    return(lapply(c(list(x), x$patches$plots), function(p)
      try(ggplot2::ggplot_build(p)$data, silent = TRUE)))
  }
  if (inherits(x, "ggplot"))  return(ggplot2::ggplot_build(x)$data)
  if (inherits(x, "gtable"))  return(lapply(x$grobs, class))
  if (inherits(x, "theme"))   return(names(x))
  if (isS4(x)) {
    sl <- methods::slotNames(x)
    if ("result" %in% sl) return(methods::slot(x, "result"))
    return(lapply(stats::setNames(sl, sl), function(s) methods::slot(x, s)))
  }
  x
}

# Column order inside a built layer carries no meaning -- ggplot2 emits aesthetic
# columns in the order the aes() mapping happened to name them, so an equivalent
# rewrite of a geom call reorders them without changing a single value. all.equal()
# compares data frames positionally, which turns that into a spurious "target is
# character, current is numeric". Canonicalising the order is not a loosened
# tolerance: every value is still compared exactly, and comparing like-named
# columns is strictly stricter than comparing by position.
canon <- function(x) {
  if (is.data.frame(x)) return(x[, order(names(x)), drop = FALSE])
  if (is.list(x)) return(lapply(x, canon))
  x
}

CASES  <- list()
case   <- function(name, expr) CASES[[name]] <<- substitute(expr)

# ---- helpers / small utilities ---------------------------------------------
case("empty_gsea_tibble",   empty_gsea_tibble())
case("list_to_term2gene",   list_to_term2gene(gsets$alpha))
case("filter_by_size",      filter_by_size(gsets$alpha, min_size = 20, max_size = 60))
case("format_pathway_name", format_pathway_name(
       c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "GOBP_T_CELL_ACTIVATION",
         "REACTOME_CELL_CYCLE_MITOTIC", "ALPHA_SIGNAL_UP")))
case("ensure_dir", { d <- file.path(tempdir(), "golden_ensure_dir_probe")
                     unlink(d, recursive = TRUE); ensure_dir(d); dir.exists(d) })

# ---- data / IO --------------------------------------------------------------
genes_df <- rd("gene_annotation.rds")
rownames(genes_df) <- genes_df$ensembl_id
case("build_dge", { g <- build_dge(counts, meta, genes_df[rownames(counts), ])
                    list(dim = dim(g), samples = colnames(g),
                         norm_factors = round(g$samples$norm.factors, 6)) })

# ---- reference databases (real bundled data) --------------------------------
case("list_reference_dbs", list_reference_dbs(toolkit_dir = "."))
case("load_reference_db_mitopathways", {
       db <- load_reference_db("mitopathways", "Mus_musculus", toolkit_dir = ".")
       lapply(db, function(x) utils::head(x, 50)) })
case("load_reference_db_transportdb", {
       db <- load_reference_db("transportdb", "Mus_musculus", toolkit_dir = ".")
       lapply(db, function(x) utils::head(x, 50)) })

# ---- GSEA on real symbols ---------------------------------------------------
case("run_gsea_hallmark", run_gsea(de_real, rank_metric = "t",
                                   species = "Mus musculus", collection = "H",
                                   nperm = 10000, seed = 123))
case("normalize_gsea_results", {
  g <- run_gsea(de_real, rank_metric = "t", species = "Mus musculus",
                collection = "H", nperm = 10000, seed = 123)
  normalize_gsea_results(g, database = "Hallmark", contrast = "KO_vs_WT") })

# ---- plots: DE --------------------------------------------------------------
case("create_standard_volcano_fdr",
     create_standard_volcano(tt, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1))
case("create_standard_volcano_p",
     create_standard_volcano(tt, decision_by = "p", p_cutoff = 0.05, fc_cutoff = 1))
case("create_standard_volcano_counts",
     create_standard_volcano(tt, decision_by = "fdr", p_cutoff = 0.05,
                             fc_cutoff = 1, annotate_counts = TRUE))
# create_MD_plot takes a limma fit + coef, so refit from the fixture intermediates.
ii   <- rd("de_intermediates.rds")
.fit <- limma::eBayes(limma::contrasts.fit(
          limma::lmFit(ii$logcpm, ii$design), ii$contrast))
case("create_MD_plot", create_MD_plot(.fit, coef = "KO_vs_WT", de_results = tt))
case("custom_minimal_theme_with_grid", custom_minimal_theme_with_grid())

# ---- plots: GSEA ------------------------------------------------------------
.gsea_obj <- NULL
case("gsea_dotplot",       gsea_dotplot(.gsea_obj, showCategory = 15, padj_cutoff = 0.05))
case("gsea_dotplot_facet", gsea_dotplot_facet(.gsea_obj, showCategory = 10))
case("gsea_barplot",       gsea_barplot(.gsea_obj, top_n = 15))
case("gsea_running_sum_plot",
     gsea_running_sum_plot(.gsea_obj,
       gene_set_ids = c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                        "HALLMARK_MYC_TARGETS_V1")))

# ---- run ---------------------------------------------------------------------
.gsea_obj <- run_gsea(de_real, rank_metric = "t", species = "Mus musculus",
                      collection = "H", nperm = 10000, seed = 123)

manifest <- data.frame(case = character(), status = character(),
                       note = character(), stringsAsFactors = FALSE)

for (nm in names(CASES)) {
  if (!is.null(ONLY) && !nm %in% ONLY) next
  set.seed(123)
  res <- try(suppressWarnings(suppressMessages(eval(CASES[[nm]], globalenv()))),
             silent = TRUE)
  if (inherits(res, "try-error")) {
    manifest <- rbind(manifest, data.frame(case = nm, status = "ERROR",
                        note = trimws(sub("\n.*", "", as.character(res)))))
    cat(sprintf("  %-34s ERROR\n", nm))
    next
  }
  val <- try(canon(normalize(res)), silent = TRUE)
  if (inherits(val, "try-error")) {
    manifest <- rbind(manifest, data.frame(case = nm, status = "NORMALIZE_FAIL",
                        note = trimws(sub("\n.*", "", as.character(val)))))
    cat(sprintf("  %-34s NORMALIZE_FAIL\n", nm)); next
  }
  path <- file.path(OUT, paste0(nm, ".rds"))
  if (MODE == "capture") {
    saveRDS(val, path)
    manifest <- rbind(manifest, data.frame(case = nm, status = "OK",
                        note = paste(class(res), collapse = "/")))
    cat(sprintf("  %-34s OK   (%s)\n", nm, paste(class(res), collapse = "/")))
  } else {
    if (!file.exists(path)) {
      manifest <- rbind(manifest, data.frame(case = nm, status = "NEW",
                          note = "no golden on record"))
      cat(sprintf("  %-34s NEW\n", nm)); next
    }
    cmp <- all.equal(canon(readRDS(path)), val, tolerance = 1e-6,
                     check.attributes = FALSE)
    if (isTRUE(cmp)) {
      manifest <- rbind(manifest, data.frame(case = nm, status = "PASS", note = ""))
      cat(sprintf("  %-34s PASS\n", nm))
    } else {
      manifest <- rbind(manifest, data.frame(case = nm, status = "FAIL",
                          note = paste(utils::head(cmp, 3), collapse = " | ")))
      cat(sprintf("  %-34s FAIL  %s\n", nm, utils::head(cmp, 1)))
    }
  }
}

# Functions deliberately not captured, with the reason.
skipped <- data.frame(
  case = c("download_gatom_references", "run_gsea_analysis", "plot_all_gsea_results",
           "save_gsea_log", "parse_gmx", "parse_mitoxplorer", "convert_human_to_mouse"),
  status = "SKIPPED",
  note = c("network download",
           "pipeline: writes files, wraps run_gsea + plotters already captured",
           "pipeline: writes files, wraps plotters already captured",
           "writes a log file; no return value to compare",
           "needs a raw .gmx excluded from the built package",
           "needs a raw mitoXplorer file excluded from the built package",
           "homologene lookup; covered indirectly, no fixture input"),
  stringsAsFactors = FALSE)
manifest <- rbind(manifest, skipped)

if (is.null(ONLY)) {
  write.csv(manifest, file.path("tests/golden",
            if (MODE == "capture") "manifest.csv" else "verify-report.csv"),
            row.names = FALSE)
} else {
  cat("\nSelective capture: manifest.csv left untouched.\n")
}
if (MODE == "capture") {
  cat("\n", sum(manifest$status == "OK"), " captured, ",
      sum(manifest$status == "ERROR"), " errors, ",
      sum(manifest$status == "SKIPPED"), " skipped\n", sep = "")
} else {
  nf <- sum(manifest$status %in% c("FAIL", "ERROR", "NORMALIZE_FAIL", "NEW"))
  cat("\n", sum(manifest$status == "PASS"), " pass, ", nf, " fail/new\n", sep = "")
  if (nf > 0) { cat("GOLDEN VERIFY FAILED\n"); quit(status = 1) }
  cat("GOLDEN VERIFY PASSED\n")
}
cat("R ", as.character(getRversion()), " | ggplot2 ",
    as.character(utils::packageVersion("ggplot2")), " | fgsea ",
    as.character(utils::packageVersion("fgsea")), "\n", sep = "")
