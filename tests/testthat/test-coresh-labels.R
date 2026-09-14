# Contract tests for coresh_labels() and the CoReSh naming fixes.
#
# The defect these guard: a CoReSh bar was labelled with its own set id,
# "Coresh Q De Wt T16 Vs T0 Gse233217". That repeats the database name on a
# CoReSh panel, exposes the internal query-type token, and shows the caller's
# own contrast name against an unrelated external GEO accession.

cl_hits <- function() {
  tibble::tibble(
    set_name = c("CORESH_Q_de_WT_t16_vs_t0_GSE233217",
                 "CORESH_Q_curated_ifn_response_GSE233216"),
    query_name = c("Q_de_WT_t16_vs_t0", "Q_curated_ifn_response"),
    gse = c("GSE233217", "GSE233216"),
    gpl = c("GPL21103", "GPL24247"),
    pct_var = c(40.12, 35.5),
    query_size = c(20L, 7L),
    rank_in_coresh = c(1L, 2L)
  )
}

test_that("coresh_labels() composes the web-UI metadata line", {
  out <- coresh_labels(cl_hits())
  expect_named(out, cl_hits()$set_name)
  expect_identical(out[[1L]], "GSE233217 · GPL21103 · n = 20 · 40.1% var")
  # The label must not carry the database name or the query-type token.
  expect_false(any(grepl("CORESH|Coresh", out)))
  expect_false(any(grepl("Q_de|Q_curated", out)))
  # Nor the caller's contrast name, which is the confusing part.
  expect_false(any(grepl("WT_t16_vs_t0", out)))
})

test_that("coresh_labels() prepends a GEO title when one is supplied", {
  titles <- c(GSE233217 = "Interferon response in bone marrow macrophages")
  out <- coresh_labels(cl_hits(), titles = titles)
  expect_match(out[[1L]], "^Interferon response in bone marrow macrophages — ")
  expect_match(out[[1L]], "GSE233217")
  # A dataset with no title falls back to metadata alone, not to an empty gap.
  expect_identical(out[[2L]], "GSE233216 · GPL24247 · n = 7 · 35.5% var")
})

test_that("coresh_labels() accepts a data frame title lookup", {
  titles <- data.frame(gse = "GSE233216", title = "Curated ISG panel",
                       stringsAsFactors = FALSE)
  out <- coresh_labels(cl_hits(), titles = titles)
  expect_match(out[[2L]], "^Curated ISG panel — ")
  expect_error(
    coresh_labels(cl_hits(), titles = data.frame(gse = "x")),
    "must have `gse` and `title` columns"
  )
  expect_error(coresh_labels(cl_hits(), titles = "no names"),
               "must be a named character vector")
})

test_that("coresh_labels() truncates a long title but never a metadata field", {
  long <- paste(rep("word", 60L), collapse = " ")
  out <- coresh_labels(cl_hits(), titles = c(GSE233217 = long),
                       title_width = 20L)
  expect_match(out[[1L]], "…")
  expect_match(out[[1L]], "GSE233217 · GPL21103")
  whole <- coresh_labels(cl_hits(), titles = c(GSE233217 = long),
                         title_width = Inf)
  expect_false(grepl("…", whole[[1L]]))
  expect_error(coresh_labels(cl_hits(), title_width = 2),
               "at least 8")
})

test_that("coresh_labels() styles and partial records both stay readable", {
  expect_identical(
    unname(coresh_labels(cl_hits(), style = "compact")[[1L]]),
    "GSE233217 · GPL21103 · 40.1% var"
  )
  expect_identical(
    unname(coresh_labels(cl_hits(), style = "title")[[1L]]),
    "GSE233217"
  )
  # A record missing gpl and pct_var must not leave dangling separators.
  bare <- tibble::tibble(gse = "GSE1")
  expect_identical(unname(coresh_labels(bare)[[1L]]), "GSE1")
})

test_that("coresh_labels() validates its input and handles zero rows", {
  expect_error(coresh_labels(list(gse = "GSE1")), "must be a data frame")
  expect_error(coresh_labels(tibble::tibble(x = 1)),
               "must have a `gse` column")
  expect_length(coresh_labels(cl_hits()[0, , drop = FALSE]), 0L)
})

test_that("format_pathway_name() strips the package's own CORESH prefix", {
  out <- format_pathway_name("CORESH_Q_curated_ifn_response_GSE233216")
  expect_false(grepl("Coresh", out))
  # The query name is the caller's, so it survives -- coresh_labels() is the
  # way to stop showing ids at all.
  expect_match(out, "IFN")
})

test_that("a placeholder GMT description falls back to the set id", {
  # coresh_derived_sets.gmt ships "-" in every description field, which used to
  # render as a literal dash on every bar.
  gmt <- c(
    "SET_A\t-\tG1\tG2\tG3",
    "SET_B\t.\tG4\tG5",
    "SET_C\tA real description\tG6\tG7"
  )
  parsed <- bulkiRNA:::.gsdb_parse_gmt(gmt)
  expect_identical(unname(parsed$labels[["SET_A"]]), "SET_A")
  expect_identical(unname(parsed$labels[["SET_B"]]), "SET_B")
  expect_identical(unname(parsed$labels[["SET_C"]]), "A real description")
})

# --- regression: prose labels mangled by the bar/dot/heatmap renderers ------
# .gs_plot_frame() applied format_pathway_name() unconditionally. That function
# is built for ALL_CAPS_SNAKE ids and is not idempotent on prose: it maps "." to
# a space and title-cases every word. A coresh_labels() label therefore arrived
# as "Gse174808 . 40 1% Var" -- the decimal point eaten, the accession
# lower-cased. gs_plot_running() had been selective since 1.0.0; the other three
# renderers had not.

cl_result <- function(names_are_prose) {
  ids <- c("CORESH_Q_de_A_GSE1", "CORESH_Q_de_B_GSE2")
  nm <- if (names_are_prose) {
    c("GSE1 · GPL570 · n = 20 · 40.1% var",
      "GSE2 · GPL1261 · n = 7 · 3.5% var")
  } else {
    ids
  }
  bulkiRNA:::gs_result(
    data.frame(
      pathway_id = ids, pathway_name = nm,
      n_genes = c(20L, 7L), n_genes_tested = c(20L, 7L),
      stat = c(2.1, -1.8), p_value = c(0.001, 0.01), padj = c(0.01, 0.04),
      stringsAsFactors = FALSE
    ),
    database = "CoReSh", contrast = "KO-WT", method = "fgsea",
    stat_type = "NES"
  )
}

test_that("gs_plot_bar() leaves a real display name exactly as given", {
  p <- gs_plot_bar(cl_result(TRUE), top_n = 2)
  labs <- as.character(p$data$label)
  expect_true(any(grepl("GSE1", labs, fixed = TRUE)))
  # The three specific corruptions, each asserted absent.
  expect_false(any(grepl("Gse1", labs, fixed = TRUE)))
  expect_false(any(grepl("40 1", labs, fixed = TRUE)))
  expect_false(any(grepl("% Var", labs, fixed = TRUE)))
  expect_true(any(grepl("40.1% var", labs, fixed = TRUE)))
})

test_that("gs_plot_bar() still formats a name that is only the machine id", {
  # The selective rule must not stop ids being made readable -- that was the
  # point of formatting them in the first place.
  res <- cl_result(FALSE)
  res$pathway_id <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "KEGG_APOPTOSIS")
  res$pathway_name <- res$pathway_id
  labs <- as.character(gs_plot_bar(res, top_n = 2)$data$label)
  expect_false(any(grepl("_", labs, fixed = TRUE)))
  expect_false(any(grepl("HALLMARK|KEGG", labs)))
})

test_that("gs_plot_dot() and gs_plot_heatmap() share the selective rule", {
  # One helper serves all three renderers, so the fix must reach all three.
  dot_labs <- as.character(gs_plot_dot(cl_result(TRUE), top_n = 2)$data$label)
  expect_true(any(grepl("40.1% var", dot_labs, fixed = TRUE)))
  hm_labs <- as.character(
    gs_plot_heatmap(cl_result(TRUE), top_n = 2)$data$label
  )
  expect_true(any(grepl("40.1% var", hm_labs, fixed = TRUE)))
})

# --- the figure's sidecar keeps the provenance ------------------------------
# The chosen resolution for CoReSh labelling: a compact metadata axis, with the
# full GEO title and every other provenance column in the .tsv beside the
# figure. That only works if extra columns on the result survive into the
# plotted frame, which is what gs_save() writes.

test_that("extra result columns reach the figure's source table", {
  res <- cl_result(TRUE)
  res$gse <- c("GSE1", "GSE2")
  res$gpl <- c("GPL570", "GPL1261")
  res$pct_var <- c(40.12, 3.5)
  res$geo_title <- c(
    "The role of EGR1 in hypoxia/reoxygenation (H/R) injuries",
    "Glycolytic shift during West Nile virus infection"
  )

  p <- gs_plot_bar(res, top_n = 2)
  src <- attr(p, "gs_source")
  expect_true(all(c("gse", "gpl", "pct_var", "geo_title") %in% names(src)))
  # Aligned by pathway_id, not by position: .gs_select_top() reorders rows.
  i <- match("CORESH_Q_de_A_GSE1", src$pathway_id)
  expect_identical(src$gse[[i]], "GSE1")
  expect_identical(src$pct_var[[i]], 40.12)
  expect_match(src$geo_title[[i]], "^The role of EGR1")

  # And the title is written in full, never truncated, unlike an axis label.
  dir <- withr::local_tempdir()
  gs_save(p, file.path(dir, "bar"), formats = character(0L))
  tsv <- utils::read.delim(file.path(dir, "bar.tsv"), stringsAsFactors = FALSE)
  expect_true("geo_title" %in% names(tsv))
  expect_true(any(grepl("reoxygenation", tsv$geo_title, fixed = TRUE)))
})

test_that("carrying extra columns does not disturb the plotted aesthetics", {
  bare <- gs_plot_bar(cl_result(TRUE), top_n = 2)
  rich <- cl_result(TRUE)
  rich$anything <- c("a", "b")
  rich$a_list <- list(c("x", "y"), "z")
  with_extra <- gs_plot_bar(rich, top_n = 2)

  b1 <- ggplot2::ggplot_build(bare)$data
  b2 <- ggplot2::ggplot_build(with_extra)$data
  expect_equal(b1, b2)
})
