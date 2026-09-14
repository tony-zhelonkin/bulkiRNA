test_that("de_volcano returns a ggplot and validates its input", {
  de <- fake_de_table()
  expect_s3_class(de_volcano(de), "ggplot")
  expect_error(de_volcano(de[, c("logFC", "B")]), "missing column")
})

test_that("the FDR dashed line sits at the raw p of the FDR boundary", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1)

  boundary <- max(de$P.Value[de$adj.P.Val <= 0.05])
  expect_equal(hline_yintercepts(p), -log10(boundary))
  # ... and NOT at -log10(p_cutoff), which is the bug this guards against.
  expect_false(isTRUE(all.equal(hline_yintercepts(p), -log10(0.05))))
})

test_that("decision_by = 'p' puts the line at -log10(p_cutoff)", {
  de <- fake_de_table()
  p  <- de_volcano(de, decision_by = "p", p_cutoff = 0.05, fc_cutoff = 1)
  expect_equal(hline_yintercepts(p), -log10(0.05))
})

test_that("fixed_p_boundary overrides the derived boundary", {
  de <- fake_de_table()
  p  <- de_volcano(de, fixed_p_boundary = 1e-3)
  expect_equal(hline_yintercepts(p), 3)
})

test_that("no significant genes means no line and an annotation instead", {
  de <- fake_de_table()
  de$adj.P.Val <- 0.9
  p <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05)

  expect_length(hline_yintercepts(p), 0L)
  labels <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) d$label))
  expect_true(any(grepl("No genes pass FDR", labels)))
  expect_match(p$labels$caption, "No genes pass FDR")
})

test_that("the caption reports the realised p, not the FDR cutoff", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1)
  boundary <- signif(max(de$P.Value[de$adj.P.Val <= 0.05]), 2)
  expect_match(p$labels$caption, sprintf("p \u2264 %.2g", boundary), fixed = TRUE)
})

test_that("annotate_counts appends counts to the populated legend entry", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, fc_cutoff = 1, annotate_counts = TRUE)
  labs <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$labels
  expect_true(any(grepl("\u2191", labs)))
})

test_that(".volcano_sig_counts walks the priority order", {
  df <- data.frame(logFC = c(2, -2, 0.1), sig_dec = c(TRUE, TRUE, FALSE),
                   sig_fc = c(TRUE, TRUE, FALSE))
  expect_equal(.volcano_sig_counts(df)$category, "both")
  expect_equal(.volcano_sig_counts(df)$n_up, 1L)

  df$sig_fc <- FALSE
  expect_equal(.volcano_sig_counts(df)$category, "fdr")

  df$sig_dec <- FALSE
  none <- .volcano_sig_counts(df)
  expect_equal(none$category, "none")
  expect_true(is.na(none$legend_key))
})

test_that("vertical orientation swaps the axes and the boundary line", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, fc_cutoff = 1, orientation = "vertical")
  boundary <- max(de$P.Value[de$adj.P.Val <= 0.05])

  expect_equal(vline_xintercepts(p), -log10(boundary))
  expect_setequal(hline_yintercepts(p), c(-1, 1))
  expect_equal(p$labels$y, "log2(FC)")
})

test_that("highlight_gene gets its own always-on label layer", {
  de <- fake_de_table()
  p  <- de_volcano(de, highlight_gene = "Gene12")
  repel <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(repel, 2L)
  bold <- Filter(function(l) identical(l$aes_params$fontface, "bold"), repel)
  expect_length(bold, 1L)
  # the highlight layer carries only the highlighted gene, and never drops it
  expect_equal(rownames(bold[[1]]$data), "Gene12")
  expect_identical(bold[[1]]$geom_params$max.overlaps, Inf)
})

test_that("de_volcano_grid combines panels on shared limits", {
  skip_if_not_installed("patchwork")
  de <- fake_de_table()
  p  <- de_volcano(de, fc_cutoff = 1, orientation = "vertical")
  g  <- de_volcano_grid(list(A = p, B = p))
  expect_s3_class(g, "patchwork")
  expect_error(de_volcano_grid(list()), "non-empty list")
})

test_that("de_volcano_grid applies shared limits without ggplot2 chatter", {
  skip_if_not_installed("patchwork")
  # Two panels with DIFFERENT fold-change ranges. The narrow panel must end up
  # on the wide panel's limits, so this fails if the shared coord is not
  # actually applied -- an equal-range fixture would pass either way.
  narrow <- fake_de_table()
  narrow$logFC <- narrow$logFC / 4          # range -1 .. 1
  wide   <- fake_de_table()                 # range -4 .. 4

  pn <- de_volcano(narrow, fc_cutoff = 1, orientation = "vertical")
  pw <- de_volcano(wide,   fc_cutoff = 1, orientation = "vertical")

  expect_silent(g <- de_volcano_grid(list(narrow = pn, wide = pw)))

  own_y    <- max(abs(ggplot2::layer_scales(pn)$y$range$range))
  global_y <- max(abs(ggplot2::layer_scales(pw)$y$range$range))
  expect_gt(global_y, own_y)                # the fixture is discriminating

  # read the BUILT plot: a silently ignored assignment leaves own_y here
  panel <- ggplot2::ggplot_build(g[[1]])$layout$panel_params[[1]]
  expect_equal(panel$y.range, c(-global_y, global_y) * 1.1, tolerance = 1e-8)
})

# --- regression: the grid's threshold caption ---------------------------------
# The caption was read as `ggplot_build(p)$layout$plot$labels$caption`, but
# ggplot2 4.0.3 dropped `$layout$plot`, so the read was always NULL: every panel
# then had `caption = NULL` set and the promotion branch was dead code. Both
# documented modes lost the caption, which is the only place the realised raw-p
# boundary behind the dashed lines is reported.
test_that("de_volcano_grid keeps the threshold caption in both modes", {
  skip_if_not_installed("patchwork")
  de <- fake_de_table(n = 60, n_sig = 20)
  pa <- de_volcano(de, fc_cutoff = 1)
  expect_true(nzchar(pa$labels$caption))

  # Default: promoted to a single figure-level annotation.
  g <- de_volcano_grid(list(A = pa, B = pa))
  expect_equal(g$patches$annotation$caption, pa$labels$caption)

  # keep_first_caption: kept on panel 1 instead.
  g2 <- de_volcano_grid(list(A = pa, B = pa), keep_first_caption = TRUE)
  expect_equal(g2[[1]]$labels$caption, pa$labels$caption)
})

# --- regressions from the DE review ------------------------------------------

test_that("an all-NA adj.P.Val takes the documented no-genes-pass path", {
  # `de_results$P.Value[sig_logic]` kept the NA positions, so the vector was
  # non-empty with nothing significant; max(na.rm = TRUE) returned -Inf and the
  # line became NaN with draw_line = TRUE, captioning the figure "p <= NaN".
  de <- fake_de_table()
  de$adj.P.Val <- NA_real_
  expect_silent(p <- de_volcano(de, fc_cutoff = 1))
  expect_length(vline_xintercepts(p), 2L)      # the fold-change lines only
  expect_length(hline_yintercepts(p), 0L)      # no p boundary line
  expect_false(grepl("NaN", p$labels$caption, fixed = TRUE))
  ann <- Filter(function(l) inherits(l$geom, "GeomText") &&
                  !inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_true(any(grepl("No genes pass FDR",
                        vapply(ann, function(l) l$aes_params$label %||% "",
                               character(1L)))))
})

test_that("a p-value of zero does not make the axis limit infinite", {
  de <- fake_de_table()
  de$P.Value[1] <- 0
  p <- de_volcano(de, fc_cutoff = 1)
  yr <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
  expect_true(all(is.finite(yr)))
})

test_that("x_breaks follows the fold-change axis in both orientations", {
  # `x_breaks` is documented as fold-change tick spacing, but in vertical mode
  # fold change is on y: it used to retune the -log10(p) axis instead.
  de <- fake_de_table()
  h <- de_volcano(de, fc_cutoff = 1, x_breaks = 0.5)
  v <- de_volcano(de, fc_cutoff = 1, x_breaks = 0.5, orientation = "vertical")

  fc_breaks <- function(p, axis) {
    sc <- ggplot2::layer_scales(p)[[axis]]
    b <- sc$get_breaks()
    b <- b[is.finite(b)]
    unique(round(diff(sort(b)), 8))
  }
  expect_identical(fc_breaks(h, "x"), 0.5)
  expect_identical(fc_breaks(v, "y"), 0.5)
})
