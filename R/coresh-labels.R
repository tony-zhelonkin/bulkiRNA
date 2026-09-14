#' Web-UI-style labels for CoReSh-derived gene sets
#'
#' A CoReSh set id encodes its provenance -- `CORESH_<query_name>_<GSE>` --
#' which makes it a poor axis label twice over. It repeats the database name on
#' a figure that is already a CoReSh figure, and it shows the *query* name,
#' which for a DE-seeded query is the caller's own contrast name, against an
#' unrelated external GEO dataset. A reader sees a familiar contrast label and
#' reasonably concludes the bar is about that contrast.
#'
#' This composes the label the CoReSh web interface shows instead, from the
#' columns [coresh_search()] and [coresh_sets()] already return: GEO title (if
#' you supply one), accession, platform, query size and percent of variation.
#'
#' @param x A data frame carrying at least `gse`. Recognised optional columns:
#'   `gpl`, `pct_var`, `p_value`, `query_size` or `size`, `rank` or
#'   `rank_in_coresh`, and `set_name`. Both a [coresh_search()] result and a
#'   [coresh_sets()] `provenance` attribute satisfy this.
#' @param titles Optional GEO titles: a named character vector (`gse` ->
#'   title), or a data frame with `gse` and `title` columns. Titles are not
#'   fetched -- there is no network call here -- so pass the lookup your
#'   project already keeps.
#' @param style One of:
#'   * `"ui"` (default) -- title, then accession, platform, size and percent of
#'     variation, mirroring the web interface's columns.
#'   * `"compact"` -- accession, platform and percent of variation only.
#'   * `"title"` -- the GEO title alone, falling back to the accession.
#' @param title_width Integer. Titles longer than this are truncated with an
#'   ellipsis, because a GEO title can run to 200 characters. `Inf` keeps them
#'   whole; renderers wrap rather than truncate, so `Inf` is safe there.
#' @return A character vector of labels, named by `set_name` when `x` has that
#'   column and by `gse` otherwise. Feed it to `gsdb_register(pathway_names =)`
#'   or to a renderer's `labels =`.
#' @seealso [coresh_search()], [coresh_sets()], [gsdb_register()]
#' @examples
#' hits <- tibble::tibble(
#'   set_name = c("CORESH_Q_iron_GSE1", "CORESH_Q_iron_GSE2"),
#'   gse = c("GSE1", "GSE2"),
#'   gpl = c("GPL570", "GPL1261"),
#'   pct_var = c(0.4012, 0.3550),
#'   query_size = c(20L, 20L)
#' )
#' coresh_labels(hits)
#' coresh_labels(hits, titles = c(GSE1 = "Iron overload in macrophages"))
#' @export
coresh_labels <- function(x, titles = NULL,
                          style = c("ui", "compact", "title"),
                          title_width = 60L) {
  style <- match.arg(style)
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame from `coresh_search()` or the ",
         "`provenance` attribute of `coresh_sets()`.", call. = FALSE)
  }
  if (!"gse" %in% names(x)) {
    stop("`x` must have a `gse` column; got: ",
         paste(names(x), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(title_width) || length(title_width) != 1L ||
        is.na(title_width) || title_width < 8) {
    stop("`title_width` must be a single number of at least 8, or Inf.",
         call. = FALSE)
  }
  if (!nrow(x)) return(stats::setNames(character(0L), character(0L)))

  gse <- as.character(x$gse)
  lookup <- .coresh_title_lookup(titles)
  title <- .coresh_truncate(unname(lookup[gse]), title_width)

  gpl <- if ("gpl" %in% names(x)) as.character(x$gpl) else rep(NA_character_,
                                                               length(gse))
  size <- if ("query_size" %in% names(x)) {
    x$query_size
  } else if ("size" %in% names(x)) {
    x$size
  } else {
    rep(NA_integer_, length(gse))
  }
  pct <- if ("pct_var" %in% names(x)) as.numeric(x$pct_var) else
    rep(NA_real_, length(gse))

  out <- switch(
    style,
    title = ifelse(is.na(title) | !nzchar(title), gse, title),
    compact = .coresh_join_parts(list(gse, gpl, .coresh_fmt_pct(pct))),
    ui = {
      # The web interface reports percent of variation as a percentage; the
      # search returns it already multiplied by 100, so it is not rescaled
      # here. Reporting a fraction as a percentage was worth guarding against.
      meta <- .coresh_join_parts(list(
        gse, gpl,
        ifelse(is.na(size), NA_character_, paste0("n = ", size)),
        .coresh_fmt_pct(pct)
      ))
      ifelse(is.na(title) | !nzchar(title), meta, paste0(title, " \u2014 ", meta))
    }
  )

  keys <- if ("set_name" %in% names(x)) as.character(x$set_name) else gse
  stats::setNames(out, keys)
}

#' Normalise a GEO title lookup
#'
#' @param titles `NULL`, a named character vector, or a data frame with `gse`
#'   and `title` columns.
#' @return A named character vector, possibly empty.
#' @keywords internal
.coresh_title_lookup <- function(titles) {
  if (is.null(titles)) return(stats::setNames(character(0L), character(0L)))
  if (is.data.frame(titles)) {
    if (!all(c("gse", "title") %in% names(titles))) {
      stop("A data frame `titles` must have `gse` and `title` columns; got: ",
           paste(names(titles), collapse = ", "), ".", call. = FALSE)
    }
    return(stats::setNames(as.character(titles$title),
                           as.character(titles$gse)))
  }
  if (!is.character(titles) || is.null(names(titles))) {
    stop("`titles` must be a named character vector (gse -> title), a data ",
         "frame with `gse` and `title`, or NULL.", call. = FALSE)
  }
  titles
}

#' Truncate a title, marking that it was truncated
#'
#' @param x Character vector, possibly containing `NA`.
#' @param width Maximum width, or `Inf`.
#' @return A character vector of the same length.
#' @keywords internal
.coresh_truncate <- function(x, width) {
  if (!is.finite(width)) return(x)
  long <- !is.na(x) & nchar(x) > width
  x[long] <- paste0(substr(x[long], 1L, max(width - 1L, 1L)), "\u2026")
  x
}

#' Format percent of variation
#'
#' @param pct Numeric vector of percentages.
#' @return A character vector, `NA` where `pct` is missing.
#' @keywords internal
.coresh_fmt_pct <- function(pct) {
  ifelse(is.na(pct), NA_character_,
         paste0(formatC(pct, format = "f", digits = 1L), "% var"))
}

#' Join label parts, dropping the ones that are missing
#'
#' A missing platform or percentage leaves no empty separator behind, so a
#' partial record still reads as a label rather than as a formatting bug.
#'
#' @param parts List of equal-length character vectors.
#' @return A character vector.
#' @keywords internal
.coresh_join_parts <- function(parts) {
  n <- max(vapply(parts, length, integer(1L)))
  vapply(seq_len(n), function(i) {
    vals <- vapply(parts, function(p) as.character(p[[i]]), character(1L))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    paste(vals, collapse = " \u00b7 ")
  }, character(1L))
}
