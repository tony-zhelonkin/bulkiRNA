#' Write a `gs_result` to the standard table layout
#'
#' `gs_write()` and [gs_read()] are the only functions in the compute layer
#' that touch disk. The layout under `dir` is the project convention:
#'
#' ```
#' <dir>/
#'   by_contrast/<contrast>/<name>_<database>.tsv
#'   _overview/<name>_all.tsv
#'   _overview/<name>_summary.tsv
#' ```
#'
#' `dir` is normally `03_results/<stage>/tables`. List columns
#' (`leading_edge`) are collapsed to `"/"`-separated strings on the way out and
#' split again by [gs_read()], so the round trip is lossless.
#'
#' This layout is for `gs_result` only. A [gs_matrix] of per-sample scores has
#' no table format here -- persist it with [saveRDS()], which keeps its
#' `sample_data` and score-type attributes as they are.
#'
#' @param x A [gs_result].
#' @param dir Directory to write into; created if missing.
#' @param name Character stem for the file names. Default `"gsea"`.
#' @param overview Logical. Also write the pooled table and the
#'   [summary.gs_result()] table under `_overview/`.
#' @param by_contrast Logical. Also write one file per contrast x database.
#' @param prune Logical. Delete the existing `by_contrast/` tree before writing.
#'   `gs_write()` otherwise only ever adds files, so re-running a stage with a
#'   contrast dropped leaves the old directory in place for [gs_read()] to pick
#'   up. Off by default because it deletes data; a manifest is written either
#'   way, and `gs_read()` warns about files that are not in it.
#' @return The directory, invisibly, with the written paths in the `"files"`
#'   attribute.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' d <- file.path(tempdir(), "tables")
#' gs_write(res, d)
#' gs_read(d)
#' @export
gs_write <- function(x, dir, name = "gsea", overview = TRUE,
                     by_contrast = TRUE, prune = FALSE) {
  .gs_check_result(x)
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir)) {
    stop("`dir` must be a single non-empty directory path.", call. = FALSE)
  }
  ensure_dir(dir)
  files <- character(0)

  if (prune && by_contrast) {
    unlink(file.path(dir, "by_contrast"), recursive = TRUE)
  }

  if (by_contrast) {
    parts <- gs_split(x, by = c("contrast", "database"))
    for (nm in names(parts)) {
      part <- parts[[nm]]
      if (!nrow(part)) next
      cdir <- file.path(dir, "by_contrast", .gs_slug(part$contrast[1L]))
      ensure_dir(cdir)
      f <- file.path(
        cdir, paste0(name, "_", .gs_slug(part$database[1L]), ".tsv")
      )
      .gs_write_tsv(part, f)
      files <- c(files, f)
    }
  }

  if (overview) {
    odir <- file.path(dir, "_overview")
    ensure_dir(odir)
    f_all <- file.path(odir, paste0(name, "_all.tsv"))
    .gs_write_tsv(x, f_all)
    f_sum <- file.path(odir, paste0(name, "_summary.tsv"))
    utils::write.table(
      summary(x), f_sum,
      sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    files <- c(files, f_all, f_sum)
  }

  # The manifest is what makes a stale leftover visible. Without it `gs_read()`
  # globs `by_contrast/*/<name>_*.tsv` and row-binds whatever is there, so a
  # contrast dropped from the analysis reappears in the figure, dated from the
  # previous run, with nothing to indicate it.
  .gs_write_manifest(dir, name, files)

  invisible(structure(dir, files = files))
}

#' Record which files a `gs_write()` call produced
#'
#' @param dir The output directory.
#' @param name The file-name stem.
#' @param files Absolute paths just written.
#' @return The manifest path, invisibly.
#' @keywords internal
.gs_write_manifest <- function(dir, name, files) {
  path <- file.path(dir, paste0("_manifest_", .gs_slug(name), ".txt"))
  rel <- .gs_relative_to(files, dir)
  writeLines(rel, path)
  invisible(path)
}

#' Express written paths relative to the output directory
#'
#' Literal prefix removal, not a regex: a directory name may contain any of
#' `. ( ) [ ] + *`.
#'
#' @param paths Character vector of paths.
#' @param dir The output directory.
#' @return `paths` with a leading `dir/` removed where present.
#' @keywords internal
.gs_relative_to <- function(paths, dir) {
  pre <- paste0(sub("/+$", "", dir), "/")
  ifelse(startsWith(paths, pre), substring(paths, nchar(pre) + 1L), paths)
}

#' Read a `gs_result` back from the standard table layout
#'
#' Reads `_overview/<name>_all.tsv` when it exists, otherwise every
#' `by_contrast/*/<name>_*.tsv` and row-binds them. Also accepts a path to a
#' single `.tsv` file.
#'
#' @param dir Directory written by [gs_write()], or a single `.tsv` path.
#' @param name Character stem used when writing.
#' @return A [gs_result].
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' d <- file.path(tempdir(), "tables2")
#' gs_write(res, d)
#' gs_read(d)
#' @export
gs_read <- function(dir, name = "gsea") {
  if (!is.character(dir) || length(dir) != 1L) {
    stop("`dir` must be a single path.", call. = FALSE)
  }
  paths <- if (grepl("\\.tsv$", dir)) {
    dir
  } else {
    all_f <- file.path(dir, "_overview", paste0(name, "_all.tsv"))
    if (file.exists(all_f)) {
      all_f
    } else {
      found <- list.files(
        file.path(dir, "by_contrast"),
        pattern = paste0("^", name, "_.*\\.tsv$"),
        recursive = TRUE, full.names = TRUE
      )
      .gs_warn_unmanifested(dir, name, found)
      found
    }
  }
  missing <- paths[!file.exists(paths)]
  if (!length(paths) || length(missing)) {
    stop("No gs_result tables found at ", encodeString(dir, quote = "\""),
         ". Expected `_overview/", name, "_all.tsv` or ",
         "`by_contrast/*/", name, "_*.tsv`.", call. = FALSE)
  }
  parts <- lapply(paths, .gs_read_tsv)
  gs_result(dplyr::bind_rows(parts))
}

#' Warn about globbed files the last `gs_write()` did not produce
#'
#' @param dir The directory being read.
#' @param name The file-name stem.
#' @param found Absolute paths the glob returned.
#' @return `NULL`, invisibly.
#' @keywords internal
.gs_warn_unmanifested <- function(dir, name, found) {
  man <- file.path(dir, paste0("_manifest_", .gs_slug(name), ".txt"))
  if (!length(found) || !file.exists(man)) {
    return(invisible(NULL))
  }
  listed <- readLines(man, warn = FALSE)
  rel <- .gs_relative_to(found, dir)
  stale <- rel[!rel %in% listed]
  if (length(stale)) {
    warning(
      length(stale), " file(s) under ", encodeString(dir, quote = "\""),
      " were not written by the last `gs_write()` and are being read anyway: ",
      paste(utils::head(stale, 5L), collapse = ", "),
      ". Re-run `gs_write(prune = TRUE)` to clear them.", call. = FALSE
    )
  }
  invisible(NULL)
}

# ---- internals --------------------------------------------------------------

#' Write one `gs_result` as a TSV, flattening list columns
#'
#' @param x A [gs_result].
#' @param path File path.
#' @return `path`, invisibly.
#' @keywords internal
.gs_write_tsv <- function(x, path) {
  df <- as.data.frame(.as_plain_tibble(x))
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(
        df[[nm]], function(z) paste(as.character(unlist(z)), collapse = "/"),
        character(1L)
      )
    }
  }
  ensure_parent_dir(path)
  utils::write.table(
    df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  invisible(path)
}

#' Read one TSV back, restoring the `leading_edge` list column
#'
#' @param path File path.
#' @return A data frame.
#' @keywords internal
.gs_read_tsv <- function(path) {
  df <- utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "",
    stringsAsFactors = FALSE, na.strings = c("", "NA"),
    check.names = FALSE
  )
  if (!is.null(df[["leading_edge"]])) {
    df[["leading_edge"]] <- lapply(df[["leading_edge"]], function(s) {
      if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, "/")[[1L]]
    })
  }
  df
}

#' Make a path-safe slug from a label
#'
#' @param x Character vector.
#' @return Character vector with path-hostile characters replaced by `_`.
#' @keywords internal
.gs_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9._+-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}
