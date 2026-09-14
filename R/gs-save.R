#' Save a figure and its source table
#'
#' The craft rule "a figure's source table is its same-stem neighbour" is this
#' function's job, not something to remember. From one plot object `gs_save()`
#' writes, at the same stem:
#'
#' ```
#' <stem>.pdf   <stem>.png   <stem>.tsv
#' ```
#'
#' The table is the frame the plot was drawn from -- the `gs_source` attribute
#' every `gs_plot_*` renderer attaches -- or whatever you pass as `data`. List
#' columns are collapsed with `/`, matching [gs_write()].
#'
#' @param plot A ggplot object.
#' @param path Output path stem. A recognised extension (`.pdf`, `.png`,
#'   `.tsv`) is stripped; the parent directory is created if missing.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi Resolution for the raster output.
#' @param formats Character vector of image formats to write. Default
#'   `c("pdf", "png")`.
#' @param data Source table to write instead of the plot's own. Pass `NULL`
#'   (default) to use the plot's `gs_source` attribute, falling back to
#'   `plot$data`.
#' @param table Logical. Write the source table. `FALSE` writes only images.
#' @param rescale_fonts Logical. Scale text with the canvas via
#'   [gs_scale_fonts()], so a figure asked for at 9 x 10 in does not carry the
#'   same point sizes as one at 7 x 5. `TRUE` by default.
#' @param base_font_size Numeric. Reference point size at the 7 x 5 in
#'   reference canvas, used only when `rescale_fonts = TRUE`.
#' @return The written paths, invisibly.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "SET_A", SET_B = "SET_B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' p <- gs_plot_dot(res, top_n = 5)
#' gs_save(p, file.path(tempdir(), "figures", "demo_dot"))
#' @export
gs_save <- function(plot, path, width = 8, height = 6, dpi = 300,
                    formats = c("pdf", "png"), data = NULL, table = TRUE,
                    rescale_fonts = TRUE, base_font_size = 10) {
  if (!inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    stop("`plot` must be a ggplot or patchwork object; got ",
         paste(class(plot), collapse = "/"), ".", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single non-empty path stem.", call. = FALSE)
  }
  known <- c("pdf", "png", "tsv", "svg", "tiff", "jpeg")
  bad <- setdiff(formats, c("pdf", "png", "svg", "tiff", "jpeg"))
  if (length(bad)) {
    stop("Unsupported `formats`: ", paste(bad, collapse = ", "),
         ". Use any of pdf, png, svg, tiff, jpeg.", call. = FALSE)
  }

  ext <- tolower(sub("^.*\\.", "", basename(path)))
  stem <- if (ext %in% known) sub("\\.[^.]+$", "", path) else path
  ensure_parent_dir(stem)

  # Font scaling before writing, and never in the renderer: the renderer does
  # not know the canvas. Sizes are read off the object being drawn, so the
  # source table below is still taken from the unscaled plot.
  drawn <- if (isTRUE(rescale_fonts)) {
    gs_scale_fonts(plot, width = width, height = height,
                   base_font_size = base_font_size)
  } else {
    plot
  }

  written <- character(0)
  for (fmt in formats) {
    f <- paste0(stem, ".", fmt)
    ggsave(f, plot = drawn, width = width, height = height, dpi = dpi,
           units = "in")
    written <- c(written, f)
  }

  if (isTRUE(table)) {
    src <- data
    if (is.null(src)) src <- attr(plot, "gs_source")
    if (is.null(src)) src <- plot$data
    if (is.data.frame(src) && nrow(src) > 0L) {
      f <- paste0(stem, ".tsv")
      .gs_save_tsv(src, f)
      written <- c(written, f)
    }
  }
  invisible(written)
}

#' Write a plot's source table
#'
#' List columns are collapsed with `/` and factors written as their labels, so
#' the file round-trips through a spreadsheet unharmed.
#'
#' @param x A data frame.
#' @param path Output file path.
#' @return `path`, invisibly.
#' @keywords internal
.gs_save_tsv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    col <- x[[nm]]
    if (is.list(col)) {
      x[[nm]] <- vapply(col, function(v) paste(v, collapse = "/"),
                        character(1L))
    } else if (is.factor(col)) {
      x[[nm]] <- as.character(col)
    }
  }
  names(x) <- sub("^\\.", "", names(x))
  ensure_parent_dir(path)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "")
  invisible(path)
}

#' Write a plain-text log of a gene-set result
#'
#' The body behind the deprecated `save_gsea_log()`: a human-readable summary
#' of the significant pathways, split by direction.
#'
#' @param x A [gs_result].
#' @param path Output file path; the parent directory is created if missing.
#' @param padj_cutoff FDR threshold defining "significant".
#' @return `path`, invisibly.
#' @keywords internal
.gs_write_log <- function(x, path, padj_cutoff = 0.05) {
  .gs_plot_check_result(x)
  ensure_parent_dir(path)

  df <- as.data.frame(x, stringsAsFactors = FALSE)
  sig <- df[!is.na(df$padj) & df$padj < padj_cutoff, , drop = FALSE]
  up <- sig[sig$direction == "up", , drop = FALSE]
  down <- sig[sig$direction == "down", , drop = FALSE]

  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  say <- function(...) cat(..., "\n", sep = "", file = con)

  say("Gene-set results log")
  say("====================")
  say("")
  say("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  say("Databases: ", paste(unique(df$database), collapse = ", "))
  say("Contrasts: ", paste(unique(df$contrast), collapse = ", "))
  say("Method: ", paste(unique(df$method), collapse = ", "),
      "  Statistic: ", gs_stat_label(x))
  say("")
  say("Total pathways: ", nrow(df))
  say("Significant up (padj < ", padj_cutoff, "): ", nrow(up))
  say("Significant down (padj < ", padj_cutoff, "): ", nrow(down))
  say("")

  if (nrow(sig) == 0L) {
    say("No significant pathways at padj < ", padj_cutoff, ".")
    return(invisible(path))
  }

  block <- function(part, heading) {
    if (nrow(part) == 0L) return(invisible(NULL))
    part <- part[order(part$padj), , drop = FALSE]
    say(heading)
    say(strrep("-", nchar(heading)))
    for (i in seq_len(nrow(part))) {
      row <- part[i, ]
      say(i, ". ", row$pathway_name, "  [", row$pathway_id, "]")
      say("   ", gs_stat_label(x), ": ", sprintf("%.3f", row$stat))
      say("   padj: ", sprintf("%.3e", row$padj),
          "   genes tested: ", row$n_genes_tested, "/", row$n_genes)
      le <- df$leading_edge
      if (!is.null(le)) {
        say("   leading edge: ",
            paste(le[[which(df$pathway_id == row$pathway_id &
                              df$contrast == row$contrast)[1L]]],
                  collapse = "/"))
      }
      say("")
    }
  }
  block(up, "UPREGULATED PATHWAYS")
  block(down, "DOWNREGULATED PATHWAYS")
  say("End of log")
  invisible(path)
}
