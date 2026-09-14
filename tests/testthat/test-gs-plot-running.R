# Contract tests for gs_plot_running(). Plot assertions read
# ggplot_build(p)$data / $layout, never pixels.
#
# The renderer returns a patchwork of three panels, so most assertions address
# a panel by index: p[[1]] is the ES curve, p[[2]] the gene ticks, p[[3]] the
# ranked metric.

rs_ranks <- function(n = 60L) {
  stats::setNames(seq(3, -3, length.out = n), paste0("G", seq_len(n)))
}

rs_sets <- function() {
  list(
    SET_A = paste0("G", c(1:6, 30)),
    SET_B = paste0("G", c(55:60, 20)),
    SET_C = paste0("G", seq(2, 60, by = 6))
  )
}

rs_db <- function() {
  sets <- rs_sets()
  structure(
    sets,
    pathway_names = c(SET_A = "Alpha response", SET_B = "Beta response",
                      SET_C = "Gamma response"),
    database = "testdb", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
}

rs_result <- function() {
  bulkiRNA:::gs_result(
    data.frame(
      pathway_id = c("SET_A", "SET_B", "SET_C"),
      pathway_name = c("Alpha response", "Beta response", "Gamma response"),
      n_genes = c(7L, 7L, 10L),
      n_genes_tested = c(7L, 7L, 10L),
      stat = c(2.1, -1.8, 0.4),
      p_value = c(0.001, 0.01, 0.4),
      padj = c(0.01, 0.05, 0.6),
      stringsAsFactors = FALSE
    ),
    database = "testdb", contrast = "KO-WT", method = "fgsea",
    stat_type = "NES"
  )
}

# The legend lives on the ES panel, which owns the colour scale.
rs_guide <- function(p) ggplot2::get_guide_data(p[[1]], "colour")

# Panel heights as laid out, in null units.
rs_heights <- function(p) {
  g <- patchwork::patchworkGrob(p)
  as.numeric(g$heights[grepl("null", as.character(g$heights))])
}

rs_yrange <- function(p, i) {
  ggplot2::ggplot_build(p[[i]])$layout$panel_params[[1]]$y.range
}

# --- the contract whose loss broke every consumer figure --------------------
# gs_plot_running() used to return a single faceted ggplot. scio's
# style_series() branches on the object's class, so returning the wrong type
# silently routed every figure into a path that applied coord_cartesian() to a
# synthetic y scale: the 2.4:0.7:0.9 height ratio flattened to 1:1:1 and the
# gene ticks fell below 1% of panel height. Nothing failed; the figures were
# just wrong. These two assertions are the guard.

test_that("returns a patchwork carrying a restyle closure", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  expect_s3_class(p, "patchwork")
  expect_true(is.function(attr(p, "grs_restyle")))
  expect_length(rs_heights(p), 3L)
})

test_that("the three panels share one x scale", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  xr <- lapply(1:3, function(i) {
    ggplot2::ggplot_build(p[[i]])$layout$panel_params[[1]]$x.range
  })
  expect_equal(xr[[1]], xr[[2]])
  expect_equal(xr[[1]], xr[[3]])
})

test_that("panel heights honour panel_heights exactly", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       panel_heights = c(3, 1, 2))
  expect_equal(rs_heights(p), c(3, 1, 2), tolerance = 1e-6)
})

test_that("a y clamp on the ES panel leaves the layout and other panels alone", {
  # The exact operation that used to destroy the figure, applied through the
  # documented argument and through the restyle closure.
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(), es_ylim = c(-1, 1))
  expect_equal(rs_heights(p), c(2.4, 0.7, 0.9), tolerance = 1e-6)
  expect_equal(rs_yrange(p, 2), c(0, 3), tolerance = 1e-6)

  bare <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  styled <- attr(bare, "grs_restyle")(es_ylim = c(-1, 1))
  expect_s3_class(styled, "patchwork")
  expect_equal(rs_heights(styled), c(2.4, 0.7, 0.9), tolerance = 1e-6)
  # The ES panel is clamped (plus scale expansion) and the tick panel is not.
  expect_true(all(abs(rs_yrange(styled, 1)) <= 1.1 + 1e-9))
  expect_equal(rs_yrange(styled, 2), c(0, 3), tolerance = 1e-6)
})

test_that("the restyle closure honours every layout argument", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  restyle <- attr(p, "grs_restyle")

  st <- restyle(panel_heights = c(4, 1, 1), legend_position = "right",
                xticks = "all", rug_ylabels = TRUE)
  expect_equal(rs_heights(st), c(4, 1, 1), tolerance = 1e-6)
  # xticks = "all" leaves the ES panel's x axis text in place.
  expect_false(inherits(st[[1]]$theme$axis.text.x, "element_blank"))
  # rug_ylabels = TRUE keeps the lane indices.
  expect_false(inherits(st[[2]]$theme$axis.text.y, "element_blank"))

  st2 <- restyle(xticks = "bottom", rug_ylabels = FALSE)
  expect_s3_class(st2[[1]]$theme$axis.text.x, "element_blank")
  expect_s3_class(st2[[2]]$theme$axis.text.y, "element_blank")
})

test_that("gene ticks span the whole lane, whatever top_n is", {
  # Each of n lanes must be exactly 1/n of the tick panel. The previous
  # renderer drew 0.6/n inside a rescaled window, which at top_n = 5 left a
  # tick at roughly 2% of panel height.
  for (n in c(1L, 3L)) {
    p <- gs_plot_running(rs_sets(), ranks = rs_ranks(), top_n = n)
    d <- ggplot2::ggplot_build(p[[2]])$data[[1]]
    expect_equal(sort(unique(d$yend - d$y)), 1, tolerance = 1e-9)
    expect_equal(rs_yrange(p, 2), c(0, n), tolerance = 1e-9)
  }
})

test_that("panel_heights and es_ylim are validated", {
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               panel_heights = c(1, 2)),
               "three positive finite numbers")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               panel_heights = c(1, 0, 2)),
               "three positive finite numbers")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               es_ylim = c(1, -1)),
               "two finite increasing numbers")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               es_ylim = 1),
               "two finite increasing numbers")
})

test_that("colours are keyed by pathway id, not by position or label", {
  # Palette given in reverse-alphabetical order, and labels that would sort
  # the other way: a positional bug would swap these.
  pal <- c(SET_C = "#000001", SET_B = "#000002", SET_A = "#000003")
  p <- gs_plot_running(rs_db(), ranks = rs_ranks(),
                       pathways = c("SET_A", "SET_B", "SET_C"),
                       palette = pal,
                       labels = c(SET_A = "zzz", SET_C = "aaa"))
  gd <- rs_guide(p)
  expect_equal(as.character(gd$.value), c("SET_A", "SET_B", "SET_C"))
  expect_equal(gd$colour, c("#000003", "#000002", "#000001"))
  expect_equal(gd$.label, c("zzz", "Beta response", "aaa"))
})

test_that("an unnamed palette zips to the declared pathway order", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       pathways = c("SET_C", "SET_A"),
                       palette = c("#111111", "#222222"))
  gd <- rs_guide(p)
  expect_equal(gd$colour[as.character(gd$.value) == "SET_C"], "#111111")
  expect_equal(gd$colour[as.character(gd$.value) == "SET_A"], "#222222")
})

test_that("the default palette does not depend on the order ids arrive in", {
  a <- rs_guide(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                                pathways = c("SET_A", "SET_C")))
  b <- rs_guide(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                                pathways = c("SET_C", "SET_A")))
  key <- function(gd) stats::setNames(gd$colour, as.character(gd$.value))
  expect_equal(key(a)[["SET_A"]], key(b)[["SET_A"]])
  expect_equal(key(a)[["SET_C"]], key(b)[["SET_C"]])
})

test_that("a partially matching named palette warns and falls back", {
  expect_warning(
    gs_plot_running(rs_sets(), ranks = rs_ranks(),
                    pathways = c("SET_A", "SET_B"),
                    palette = c(SET_A = "#111111")),
    "do not cover every plotted pathway"
  )
})

test_that("the ES curve comes from fgsea, unaltered", {
  ranks <- rs_ranks()
  sets <- rs_sets()
  p <- gs_plot_running(sets, ranks = ranks, pathways = "SET_A")
  # Real panels mean real units: the drawn y IS the enrichment score now, not
  # a value rescaled into a synthetic window.
  drawn <- ggplot2::ggplot_build(p[[1]])$data[[2]]
  ref <- fgsea::plotEnrichmentData(pathway = sets$SET_A, stats = ranks)$curve
  expect_equal(nrow(drawn), nrow(ref))
  expect_equal(drawn$x, as.numeric(ref$rank))
  expect_equal(drawn$y, as.numeric(ref$ES), tolerance = 1e-12)
})

test_that("each panel's y axis reads in its own units", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  es <- rs_yrange(p, 1)
  expect_true(all(abs(es) <= 1.2))
  met <- rs_yrange(p, 3)
  expect_true(max(met) >= 2 && min(met) <= -2)
  # The tick panel's lane index carries no meaning, so it is hidden by default.
  expect_s3_class(p[[2]]$theme$axis.text.y, "element_blank")
})

test_that("a gs_result selects its top pathways by abs(stat) and needs a db", {
  res <- rs_result()
  ranks <- rs_ranks()
  expect_error(gs_plot_running(res, ranks = ranks), "`db` is required")
  p <- gs_plot_running(res, ranks = ranks, db = rs_db(), top_n = 2)
  gd <- rs_guide(p)
  expect_equal(sort(as.character(gd$.value)), c("SET_A", "SET_B"))
  expect_equal(gd$.label[as.character(gd$.value) == "SET_A"], "Alpha response")
})

test_that("integer pathways index gs_result rows, but not a bare set list", {
  res <- rs_result()
  p <- gs_plot_running(res, ranks = rs_ranks(), db = rs_db(),
                       pathways = c(3L, 1L))
  gd <- rs_guide(p)
  expect_equal(as.character(gd$.value), c("SET_C", "SET_A"))
  expect_error(gs_plot_running(res, ranks = rs_ranks(), db = rs_db(),
                               pathways = 99L),
               "index outside the 3 rows")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(), pathways = 1L),
               "index the rows of a `gs_result`")
})

test_that("ranks and gene sets can arrive as attributes of x", {
  x <- rs_result()
  attr(x, "ranks") <- rs_ranks()
  attr(x, "gene_sets") <- rs_sets()
  expect_s3_class(gs_plot_running(x, top_n = 1), "patchwork")
})

test_that("missing or malformed ranks error clearly", {
  expect_error(gs_plot_running(rs_sets()), "`ranks` is required")
  expect_error(gs_plot_running(rs_sets(), ranks = 1:5),
               "must be a \\*named\\* numeric vector")
})

test_that("unknown pathways and empty overlaps error", {
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               pathways = "NOPE"),
               "not found in `db`")
  expect_error(gs_plot_running(list(SET_X = c("zz1", "zz2")),
                              ranks = rs_ranks()),
               "no genes in `ranks`")
})

test_that("an empty gs_result is an error, not an empty plot", {
  res <- rs_result()[0, , drop = FALSE]
  expect_error(gs_plot_running(res, ranks = rs_ranks(), db = rs_db()),
               "no rows")
})

test_that("legend_position and metric_label reach the plot", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       legend_position = "none",
                       metric_label = "t statistic")
  expect_equal(p[[1]]$theme$legend.position, "none")
  # metric_label is the bottom panel's y axis title, since real panels have
  # axis titles rather than facet strips.
  expect_equal(p[[3]]$labels$y, "t statistic")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               metric_label = c("a", "b")),
               "single string")
})

test_that("long labels are wrapped, never truncated", {
  long <- paste(rep("verylongword", 6), collapse = " ")
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(), pathways = "SET_A",
                       labels = c(SET_A = long), max_name_length = 20)
  lab <- rs_guide(p)$.label[[1]]
  expect_true(grepl("\n", lab))
  expect_equal(gsub("\n", " ", lab), long)
})

# --- regression: raw MSigDB ids in the legend -----------------------------------
# Every other gs_plot_* renderer formats labels before display; this one did not,
# so a db with no `pathway_names` put snake_case ids straight into the legend.
test_that("legend labels are formatted, not raw MSigDB ids", {
  sets <- list(HALLMARK_P53_PATHWAY = paste0("G", 1:8),
               KEGG_APOPTOSIS       = paste0("G", 50:58))
  db <- structure(
    sets,
    pathway_names = stats::setNames(names(sets), names(sets)),
    database = "testdb", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
  p <- gs_plot_running(db, ranks = rs_ranks(),
                       pathways = c("HALLMARK_P53_PATHWAY", "KEGG_APOPTOSIS"))
  gd <- rs_guide(p)
  expect_false(any(grepl("_", gd$.label)))
  expect_false(any(grepl("HALLMARK|KEGG", gd$.label)))
})

# The formatting must be selective. format_pathway_name() is built for
# ALL_CAPS_SNAKE ids and is NOT idempotent on prose -- it turns "Beta response"
# into "beta Response" -- so applying it unconditionally would swap the raw-id
# bug for mangled capitalisation. A label the caller already made readable must
# survive untouched.
test_that("a caller-supplied label is never re-formatted", {
  p <- gs_plot_running(rs_db(), ranks = rs_ranks(),
                       pathways = c("SET_A", "SET_B"),
                       labels = c(SET_A = "Beta response", SET_B = "zzz"))
  gd <- rs_guide(p)
  expect_true("Beta response" %in% gd$.label)
  expect_false("beta Response" %in% gd$.label)
})
