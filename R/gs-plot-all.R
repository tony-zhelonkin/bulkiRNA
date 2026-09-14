#' Render the standard plot set for a gene-set result
#'
#' The body behind the deprecated `plot_all_gsea_results()`, rebuilt on the new
#' renderers: for each database in `x` it writes an up dotplot, a down dotplot,
#' a direction-faceted dotplot, a barplot and a text log under
#' `<out_dir>/<database>/`. Every figure goes through [gs_save()], so each one
#' arrives with its source table beside it.
#'
#' Internal on purpose -- the new golden path is to call the renderer you want
#' and save it. This exists so the deprecation shim has something to delegate
#' to.
#'
#' @param x A [gs_result].
#' @param out_dir Output directory; created if missing.
#' @param name Stem prepended to every file name.
#' @param top_n Number of pathways per plot.
#' @param padj_cutoff FDR threshold for highlighting and for the log.
#' @param width Figure width in inches, or `NULL` to take a per-database
#'   suggestion from [gs_plot_size()].
#' @param height Figure height in inches, or `NULL` for the same.
#' @param verbose Logical. Report progress with [message()].
#' @return A character vector of every written path, invisibly.
#' @keywords internal
.gs_plot_all <- function(x, out_dir, name = "gsea", top_n = 20,
                         padj_cutoff = 0.05, width = NULL, height = NULL,
                         verbose = FALSE) {
  .gs_plot_check_result(x)
  if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    stop("`out_dir` must be a single non-empty directory path.",
         call. = FALSE)
  }
  ensure_dir(out_dir)
  written <- character(0)

  # One fill scale for the whole set. Each renderer otherwise rescales to its
  # own maximum, so a pathway at NES 1.5 comes out pale in the Hallmark panel
  # and saturated in the GO panel of the same figure -- and a reader compares
  # the colours. The old fixed +/-3.5 scale had this property by construction.
  shared_limits <- .gs_symmetric_limits(x[["stat"]])

  for (db in unique(x[["database"]])) {
    part <- x[x[["database"]] == db, , drop = FALSE]
    if (nrow(part) == 0L) next
    db_dir <- file.path(out_dir, db)
    ensure_dir(db_dir)
    stem <- file.path(db_dir, paste0(name, "_", db))
    if (verbose) {
      message("Rendering ", db, " (", nrow(part), " pathways) into ", db_dir)
    }

    # Sizing is per-database when the caller did not fix it: a Reactome panel
    # needs more room than a Hallmark one, and one flat canvas cramps whichever
    # database has the most sets.
    size_for <- function(type) {
      s <- gs_plot_size(type, database = db)
      list(width = width %||% s$width, height = height %||% s$height)
    }

    specs <- list(
      list(
        suffix = "_up_dot",
        plot = function() {
          gs_plot_dot(part, top_n = top_n, direction = "up",
                      highlight = padj_cutoff, limits = shared_limits,
                      title = paste0(db, ": up"))
        },
        size = size_for("dot")
      ),
      list(
        suffix = "_down_dot",
        plot = function() {
          gs_plot_dot(part, top_n = top_n, direction = "down",
                      highlight = padj_cutoff, limits = shared_limits,
                      title = paste0(db, ": down"))
        },
        size = size_for("dot")
      ),
      list(
        suffix = "_facet_dot",
        plot = function() {
          gs_plot_dot(part, top_n = top_n, facet = "direction",
                      highlight = padj_cutoff, limits = shared_limits,
                      title = db)
        },
        size = size_for("facet")
      ),
      list(
        suffix = "_bar",
        plot = function() {
          gs_plot_bar(part, top_n = top_n, highlight = padj_cutoff,
                      limits = shared_limits, title = db)
        },
        size = size_for("bar")
      )
    )

    for (spec in specs) {
      written <- c(written, gs_save(
        spec$plot(), paste0(stem, spec$suffix),
        width = spec$size$width, height = spec$size$height
      ))
    }
    written <- c(
      written,
      .gs_write_log(part, paste0(stem, "_log.txt"),
                    padj_cutoff = padj_cutoff)
    )
  }
  invisible(written)
}
