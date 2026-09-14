#' Read gene sets from a GMT or GMX file
#'
#' One reader for both MSigDB text formats; the layout is sniffed rather than
#' declared. A `.gmt`/`.gmx` extension decides it outright; otherwise the file
#' is inspected -- GMX has a fixed column count in its first two rows (row 1
#' descriptions, row 2 set names, rows 3+ one gene per set per row), while GMT
#' rows are independent sets of the form
#' `name<TAB>description<TAB>gene1<TAB>...`.
#'
#' @param path Character(1) path to a GMT or GMX file.
#' @param database Character(1) registry key for the result; defaults to the
#'   file's base name reduced to snake_case (the key lands in
#'   `gs_result$database`, where it is a join key).
#' @param species Character(1) species the file's symbols belong to.
#' @param prefix Character(1) or `NULL`; prepended to every set id as
#'   `paste0(prefix, "_", id)`.
#' @param min_size,max_size Integer(1) or `NULL`; drop sets outside these
#'   bounds.
#' @param database_label Character(1) display name for renderers, or `NULL`
#'   to default to the file's base name (independent of `database`, which is
#'   the machine-typeable key -- see the `gs_db()` contract).
#' @param verbose Logical(1); message what was parsed.
#' @return A [gs_db()].
#' @examples
#' gmt <- tempfile(fileext = ".gmt")
#' writeLines(c("SET_A\tfirst set\tGene1\tGene2",
#'              "SET_B\tsecond set\tGene2\tGene3"), gmt)
#' gsdb_from_file(gmt, species = "Mus musculus")
#' @export
gsdb_from_file <- function(path,
                           database = NULL,
                           species = "Mus musculus",
                           prefix = NULL,
                           min_size = NULL,
                           max_size = NULL,
                           verbose = FALSE,
                           database_label = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("`path` must be a single file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    stop("`path` is empty: ", path, call. = FALSE)
  }

  format <- .gsdb_sniff_format(path, lines)
  parsed <- switch(
    format,
    gmt = .gsdb_parse_gmt(lines),
    gmx = .gsdb_parse_gmx(lines)
  )

  ids <- names(parsed$sets)
  if (!is.null(prefix)) {
    ids <- paste0(prefix, "_", ids)
    names(parsed$sets) <- ids
    names(parsed$labels) <- ids
  }

  db <- gs_db(
    parsed$sets,
    database = database %||% .gsdb_key_from_path(path),
    species = species,
    pathway_names = parsed$labels,
    database_label = database_label %||% basename(path)
  )
  db <- .gs_filter_size(db, min_size, max_size, verbose = verbose)
  if (verbose) {
    message(sprintf("Parsed %d sets from %s (%s).", length(db),
                    basename(path), format))
  }
  db
}

#' Derive a machine-typeable database key from a file path
#'
#' `gs_result$database` is a join key, so a file-derived default has to be
#' typeable: the base name loses its extension and every run of
#' non-alphanumeric characters becomes a single underscore.
#'
#' @param path Character(1) file path.
#' @return Character(1) snake_case key.
#' @keywords internal
.gsdb_key_from_path <- function(path) {
  key <- tools::file_path_sans_ext(basename(path))
  key <- gsub("[^A-Za-z0-9]+", "_", key)
  key <- gsub("(^_+)|(_+$)", "", key)
  if (!nzchar(key)) "genesets" else key
}

#' Sniff GMT versus GMX
#'
#' @param path Character(1) file path (used for its extension).
#' @param lines Character vector of non-blank file lines.
#' @return `"gmt"` or `"gmx"`.
#' @keywords internal
.gsdb_sniff_format <- function(path, lines) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("gmt", "gmx")) {
    return(ext)
  }
  sample <- utils::head(lines, 6L)
  fields <- strsplit(sample, "\t", fixed = TRUE)
  widths <- vapply(fields, length, integer(1L))
  is_gmx <- length(widths) >= 3L && widths[1] == widths[2] && widths[1] >= 2L &&
    all(widths[-c(1, 2)] <= widths[1])
  if (is_gmx) {
    # GMX's row 1 holds per-set descriptions and rows 3+ hold genes; a
    # description that also shows up as a gene/id in the data rows means the
    # uniform-width heuristic above is matching a GMT by coincidence (every
    # set happening to have the same number of tab fields), not a real GMX
    # layout. Fall back to GMT rather than silently mis-parsing it.
    row1 <- trimws(fields[[1]])
    data_vals <- trimws(unlist(fields[-c(1, 2)]))
    if (any(nzchar(row1) & row1 %in% data_vals)) {
      is_gmx <- FALSE
    }
  }
  fmt <- if (is_gmx) "gmx" else "gmt"
  if (!nzchar(ext)) {
    message(sprintf(
      "`%s` has no .gmt/.gmx extension; inferred %s format from its contents.",
      basename(path), toupper(fmt)))
  }
  fmt
}

#' Parse GMT lines
#'
#' @param lines Character vector of non-blank GMT lines.
#' @return `list(sets = <named list>, labels = <named chr>)`.
#' @keywords internal
.gsdb_parse_gmt <- function(lines) {
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1L)) >= 3L]
  if (!length(fields)) {
    stop("No GMT records found: every line has fewer than three ",
         "tab-separated fields.", call. = FALSE)
  }
  ids <- trimws(vapply(fields, `[[`, character(1L), 1L))
  desc <- trimws(vapply(fields, `[[`, character(1L), 2L))
  sets <- lapply(fields, function(f) trimws(f[-c(1, 2)]))
  names(sets) <- ids
  # A GMT description field is mandatory, so writers that have nothing to say
  # put a placeholder there. "-" and "." are the common ones and both used to
  # pass through, which put a literal dash on every bar of a figure.
  blank <- is.na(desc) | !nzchar(desc) |
    tolower(desc) %in% c("na", "null", "-", ".", "--", "n/a")
  desc[blank] <- ids[blank]
  list(sets = sets, labels = stats::setNames(desc, ids))
}

#' Parse GMX lines
#'
#' @param lines Character vector of non-blank GMX lines.
#' @return `list(sets = <named list>, labels = <named chr>)`.
#' @keywords internal
.gsdb_parse_gmx <- function(lines) {
  if (length(lines) < 3L) {
    stop("A GMX file needs at least three rows (descriptions, set names, ",
         "genes); got ", length(lines), ".", call. = FALSE)
  }
  desc <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  ids <- trimws(strsplit(lines[2], "\t", fixed = TRUE)[[1]])
  rows <- strsplit(lines[-c(1, 2)], "\t", fixed = TRUE)

  keep <- !is.na(ids) & nzchar(ids)
  sets <- lapply(which(keep), function(i) {
    genes <- vapply(rows, function(r) if (i <= length(r)) r[i] else NA_character_,
                    character(1L))
    trimws(genes[!is.na(genes)])
  })
  names(sets) <- ids[keep]
  labels <- vapply(which(keep), function(i) {
    d <- if (i <= length(desc)) trimws(desc[i]) else NA_character_
    if (is.na(d) || !nzchar(d)) ids[i] else d
  }, character(1L))
  list(sets = sets, labels = stats::setNames(labels, ids[keep]))
}

#' Register a user-supplied list of gene sets as a database
#'
#' Turns a plain named list of gene symbols into a first-class [gs_db()], so
#' ad-hoc signatures go through `gs_test()` / `gs_score()` exactly like the
#' bundled databases.
#'
#' @param sets Named list of character vectors: set id -> gene symbols.
#' @param database Character(1) registry key; lands in `gs_result$database`,
#'   where it is a join and filter key, so prefer a stable snake_case value
#'   such as `"my_signatures"`.
#' @param species Character(1) species the symbols belong to.
#' @param pathway_names Named character of human-readable labels, or `NULL`
#'   to use the ids.
#' @param database_label Character(1) display name for renderers, or `NULL`
#'   to reuse `database`.
#' @param min_size,max_size Integer(1) or `NULL`; drop sets outside these
#'   bounds.
#' @return A `gs_db`.
#' @examples
#' gsdb_register(
#'   list(MY_SET = c("Actb", "Gapdh"), OTHER = c("Sdha", "Ndufa1")),
#'   database = "my_signatures",
#'   species = "Mus musculus",
#'   database_label = "My signatures"
#' )
#' @export
gsdb_register <- function(sets,
                          database,
                          species = "Mus musculus",
                          pathway_names = NULL,
                          database_label = NULL,
                          min_size = NULL,
                          max_size = NULL) {
  if (!is.list(sets)) {
    stop("`sets` must be a named list of character vectors, one per gene ",
         "set.", call. = FALSE)
  }
  bad <- !vapply(sets, function(g) is.character(g) || is.factor(g),
                 logical(1L))
  if (any(bad)) {
    stop("`sets` must contain character vectors; element ",
         which(bad)[1], " is ", class(sets[[which(bad)[1]]])[1], ".",
         call. = FALSE)
  }
  db <- gs_db(sets, database = database, species = species,
              pathway_names = pathway_names,
              database_label = database_label)
  .gs_filter_size(db, min_size, max_size)
}
