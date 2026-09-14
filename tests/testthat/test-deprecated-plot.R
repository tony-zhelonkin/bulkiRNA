# The frozen dotplot/barplot/running-sum names take a bare `gseaResult` S4
# object. Neither clusterProfiler nor DOSE is a dependency of this package
# (see the handback report), so we define a minimal stand-in class locally,
# with just the three slots `.dep_gsea_to_gs_result()` reads.
methods::setClass(
  "gseaResult",
  representation(
    result = "data.frame",
    geneList = "numeric",
    geneSets = "list"
  )
)

fake_gsea_result <- function(n = 6L) {
  nes <- seq(-2.5, 2.5, length.out = n)
  df <- data.frame(
    ID = paste0("SET_", seq_len(n)),
    Description = paste0("HALLMARK_SET_", seq_len(n)),
    setSize = 20L + seq_len(n),
    enrichmentScore = nes / 2,
    NES = nes,
    pvalue = seq(1e-6, 0.4, length.out = n),
    p.adjust = seq(1e-4, 0.6, length.out = n),
    qvalues = seq(1e-4, 0.6, length.out = n),
    rank = seq_len(n),
    core_enrichment = vapply(seq_len(n), function(i) {
      paste(paste0("G", seq_len(1 + (i %% 4))), collapse = "/")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$ID
  genes <- paste0("G", 1:10)
  gene_list <- stats::setNames(seq(3, -3, length.out = length(genes)), genes)
  gene_sets <- stats::setNames(
    lapply(seq_len(n), function(i) paste0("G", seq_len(1 + (i %% 4)))),
    df$ID
  )
  methods::new("gseaResult", result = df, geneList = gene_list,
              geneSets = gene_sets)
}

test_that(".dep_gsea_to_gs_result builds a valid gs_result", {
  g <- fake_gsea_result()
  x <- .dep_gsea_to_gs_result(g)
  expect_s3_class(x, "gs_result")
  expect_equal(nrow(x), 6L)
  expect_equal(unique(x$stat_type), "NES")
  expect_equal(unique(x$database), "gsea")
  expect_identical(attr(x, "ranks"), g@geneList)
  expect_identical(attr(x, "gene_sets"), g@geneSets)
})

# The composed path that mattered in practice: `run_gsea()` (a deprecated
# shim returning a `gs_result`, not an S4 `gseaResult`) piped straight into
# each of the four renderer shims. This is the exact object shape `run_gsea()`
# hands back -- built here via `gs_test()` directly (as `run_gsea()` does
# internally) against the shared `fake_ranks()`/`fake_gs_db()` fixtures,
# so no MSigDB network fetch is needed. `.dep_gsea_to_gs_result()` must pass
# such an object through untouched rather than rejecting it as "not a
# gseaResult".
fake_run_gsea_result <- function() {
  ranks <- fake_ranks()
  db <- fake_gs_db()
  res <- gs_test(ranks, db, min_size = 5, max_size = 50)
  attr(res, "ranks") <- ranks
  attr(res, "gene_sets") <- db
  res
}

test_that(".dep_gsea_to_gs_result passes a gs_result through untouched", {
  res <- fake_run_gsea_result()
  x <- .dep_gsea_to_gs_result(res)
  expect_identical(x, res)
})

test_that("gsea_dotplot accepts run_gsea()-shaped gs_result input", {
  res <- fake_run_gsea_result()
  expect_warning(p <- gsea_dotplot(res, showCategory = 3), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("gsea_dotplot_facet accepts run_gsea()-shaped gs_result input", {
  res <- fake_run_gsea_result()
  expect_warning(p <- gsea_dotplot_facet(res, showCategory = 3), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("gsea_barplot accepts run_gsea()-shaped gs_result input", {
  res <- fake_run_gsea_result()
  expect_warning(p <- gsea_barplot(res, padj_cutoff = 1, top_n = 3),
                "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("gsea_running_sum_plot accepts run_gsea()-shaped gs_result input, using its own ranks/gene_sets attributes", {
  res <- fake_run_gsea_result()
  expect_warning(
    p <- gsea_running_sum_plot(res, gene_set_ids = c("SET_UP", "SET_DOWN")),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
})

test_that(".dep_gsea_to_gs_result rejects genuinely unsupported input", {
  expect_error(
    .dep_gsea_to_gs_result(data.frame(x = 1)),
    "gseaResult.*gs_result|gs_result.*gseaResult"
  )
})

test_that("gsea_dotplot warns and forwards to gs_plot_dot", {
  g <- fake_gsea_result()
  expect_warning(p <- gsea_dotplot(g, showCategory = 4), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("gsea_dotplot filterBy maps direction correctly", {
  g <- fake_gsea_result(n = 8L)
  p_pos <- suppressWarnings(gsea_dotplot(g, filterBy = "NES_positive",
                                         showCategory = 20))
  expect_true(all(attr(p_pos, "gs_source")$stat > 0))

  p_neg <- suppressWarnings(gsea_dotplot(g, filterBy = "NES_negative",
                                         showCategory = 20))
  expect_true(all(attr(p_neg, "gs_source")$stat < 0))
})

test_that("gsea_dotplot_facet warns and forwards with facet = direction", {
  g <- fake_gsea_result(n = 8L)
  expect_warning(
    p <- gsea_dotplot_facet(g, showCategory = 3),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  src <- attr(p, "gs_source")
  expect_setequal(as.character(unique(src$.facet_row)), c("Up", "Down"))
})

test_that("gsea_barplot warns, hard-filters by padj_cutoff, forwards", {
  g <- fake_gsea_result(n = 8L)
  expect_warning(
    p <- gsea_barplot(g, padj_cutoff = 0.6, top_n = 20),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  expect_true(all(attr(p, "gs_source")$padj < 0.6))
})

test_that("gsea_running_sum_plot warns and forwards to gs_plot_running", {
  g <- fake_gsea_result()
  expect_warning(
    p <- gsea_running_sum_plot(g, gene_set_ids = c("SET_1", "SET_2")),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
})

test_that("custom_minimal_theme_with_grid warns and forwards to theme_bulki", {
  expect_warning(th <- custom_minimal_theme_with_grid(), "deprecated")
  expect_s3_class(th, "theme")
  # The old theme had no base-size floor: its documented 12 pt default really
  # rendered at 12 pt. C2 originally forwarded to the public `theme_bulki()`,
  # whose deliberate 14 pt floor silently enlarged every legacy figure, and
  # documented that as a divergence. The integrator removed the divergence
  # instead by routing this shim through `.theme_bulki(floor = NULL)`.
  expect_equal(th$text$size, 12)

  th2 <- suppressWarnings(custom_minimal_theme_with_grid(base_size = 20))
  expect_equal(th2$text$size, 20)

  # A value under the floor is honoured on the deprecated path only; the
  # public theme keeps the floor. Without this pair, a future "simplification"
  # back to theme_bulki() would pass every other assertion here.
  expect_equal(suppressWarnings(
    custom_minimal_theme_with_grid(base_size = 9))$text$size, 9)
  expect_equal(theme_bulki(base_size = 9)$text$size, 14)
})

test_that("plot_all_gsea_results warns, writes files, forwards to .gs_plot_all", {
  g1 <- fake_gsea_result(n = 6L)
  g2 <- fake_gsea_result(n = 6L)
  out <- file.path(tempdir(), "pagr")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_warning(
    written <- plot_all_gsea_results(
      list(dbA = g1, dbB = g2), analysis_name = "an1", out_root = out,
      n_pathways = 5
    ),
    "deprecated"
  )
  expect_true(length(written) > 0L)
  expect_true(all(file.exists(written)))
})

test_that("plot_all_gsea_results returns empty invisibly for an empty list", {
  expect_warning(
    written <- plot_all_gsea_results(list(), "an1", tempdir()),
    "deprecated"
  )
  expect_equal(written, character(0))
})

test_that("save_gsea_log warns, writes a log, forwards to .gs_write_log", {
  g <- fake_gsea_result(n = 6L)
  d <- file.path(tempdir(), "sgl")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  expect_warning(
    path <- save_gsea_log(g, "log.txt", padj_cutoff = 0.6, dir = d),
    "deprecated"
  )
  expect_true(file.exists(path))
  expect_true(any(grepl("Gene-set results log", readLines(path))))
})

test_that("create_standard_volcano warns and forwards to de_volcano", {
  de <- fake_de_table()
  expect_warning(p <- create_standard_volcano(de), "deprecated")
  expect_s3_class(p, "ggplot")

  p_ref <- de_volcano(de)
  expect_equal(p$labels, p_ref$labels)
})

test_that("create_MD_plot warns and forwards to de_md_plot", {
  de <- fake_de_table()
  fit <- list(Amean = de$AveExpr)
  expect_warning(
    p <- create_MD_plot(fit, coef = "KO_vs_WT", de_results = de),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "MD plot: KO_vs_WT")
})

# --- regression: the commonest legacy call on an unremarkable result -----------
# `gs_plot_bar()` hard-filters by `padj_max` and returns an empty plot when
# nothing survives; the shim then asked `nrow(p$data)` to decide whether to
# re-sort. With `$data` a waiver that threw "argument is of length zero", so
# `gsea_barplot()` errored outright for any run with no significant pathway --
# the normal outcome at the default cutoff of 0.05. No shim test covered it,
# which is why it shipped.
test_that("gsea_barplot survives a result where nothing is significant", {
  g <- fake_gsea_result(n = 4L)
  g@result$p.adjust <- rep(0.9, 4)
  expect_warning(p <- gsea_barplot(g, padj_cutoff = 0.05, top_n = 10),
                 "deprecated")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 0L)
})

# --- regressions from the architecture review --------------------------------

# Every shim test chose a cutoff that keeps rows (padj_cutoff = 1, or 0.6), so
# the default-cutoff / nothing-significant path was untested for all four plot
# shims -- which is why gsea_barplot() shipped erroring on it. The new renderers
# had that case covered (test-gs-plot-bar.R), the shims did not: each side
# tested its own half of the seam.

test_that("the plot shims return a figure when nothing is significant", {
  g <- fake_gsea_result(n = 6L)   # padj runs 1e-4 .. 0.6

  shims <- list(
    gsea_dotplot       = function() gsea_dotplot(g, padj_cutoff = 1e-12),
    gsea_dotplot_facet = function() gsea_dotplot_facet(g, padj_cutoff = 1e-12),
    gsea_barplot       = function() gsea_barplot(g, padj_cutoff = 1e-12)
  )
  for (nm in names(shims)) {
    p <- suppressWarnings(shims[[nm]]())
    expect_s3_class(p, "ggplot")
  }
})

test_that("gsea_barplot on an empty result is a figure, not an error", {
  empty <- suppressWarnings(empty_gsea_tibble())
  p <- suppressWarnings(gsea_barplot(empty))
  expect_s3_class(p, "ggplot")
})
