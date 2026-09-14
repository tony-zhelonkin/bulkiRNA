#' Guess the field separator of a delimited text file
#'
#' Reads the header line and picks whichever of tab, comma or semicolon appears
#' most often. Beats trusting the extension, which lies routinely ("`.txt`"
#' featureCounts output is tab separated, "`.csv`" exports are sometimes
#' semicolon separated).
#'
#' @param path Path to a text file.
#' @return A single-character separator.
#' @keywords internal
.io_guess_sep <- function(path) {
  line <- .io_first_noncomment_line(path)
  if (!length(line)) {
    stop("`", path, "` is empty.", call. = FALSE)
  }
  cand <- c("\t" = "\t", "," = ",", ";" = ";")
  counts <- vapply(cand, function(s) lengths(regmatches(line, gregexpr(s, line,
    fixed = TRUE))), integer(1))
  if (max(counts) == 0L) {
    stop("Could not find a tab, comma or semicolon separator in the header of `",
         path, "`.", call. = FALSE)
  }
  unname(cand[which.max(counts)])
}

#' Count and return the leading `#`-comment lines of a text file
#'
#' `featureCounts` always writes a `# Program:featureCounts ...` line before
#' its header. Sniffing or reading line 1 unconditionally treats that comment
#' as the header/data, which breaks both. This scans from the top and returns
#' how many leading lines start with `#`, plus the first line that does not.
#'
#' @param path Path to a text file.
#' @return A list with `skip` (integer count of leading comment lines) and
#'   `line` (the first non-comment line, or `character(0)` if the file is
#'   exhausted).
#' @keywords internal
.io_skip_leading_comments <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return(list(skip = skip, line = character(0)))
    if (!startsWith(line, "#")) return(list(skip = skip, line = line))
    skip <- skip + 1L
  }
}

#' First non-`#`-comment line of a text file
#'
#' @param path Path to a text file.
#' @return A single-character line, or `character(0)` if the file has none.
#' @keywords internal
.io_first_noncomment_line <- function(path) {
  .io_skip_leading_comments(path)$line
}

#' Read a wide counts table into a matrix
#'
#' Accepts the four shapes a bulk RNA-seq counts file actually arrives in:
#'
#' * featureCounts, post-processed -- `Geneid` plus sample columns.
#' * featureCounts, raw -- `Geneid, Chr, Start, End, Strand, Length` then
#'   samples from column 7 on.
#' * Salmon gene level -- `gene_id`, `gene_name`, then samples.
#' * generic -- `gene_id` plus sample columns.
#'
#' Sample column names are stripped of directory paths and of the
#' `.bam`/`.sam`/`.sorted`/`.markdup`/`.txt` suffixes aligners leave behind.
#'
#' For the Salmon shape the original `gene_name` column is not discarded: it is
#' returned as an `input_gene_name` attribute, a character vector named by
#' **version-stripped** stable ID, so [annotate_genes()] can fall back to the
#' quantifier's own symbol when Ensembl has none. The attribute is `NA` for the
#' other shapes.
#'
#' @param path Path to the counts file (tab, comma or semicolon delimited;
#'   detected from the header).
#' @return A numeric matrix, genes in rows, with an `input_gene_name`
#'   attribute.
#' @export
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' write.table(data.frame(Geneid = c("ENSG1", "ENSG2"), S1 = c(3, 4),
#'                        S2 = c(5, 6)),
#'             f, sep = "\t", row.names = FALSE, quote = FALSE)
#' read_counts_matrix(f)
read_counts_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  n_skip <- .io_skip_leading_comments(path)$skip
  dt <- utils::read.delim(path, sep = .io_guess_sep(path), skip = n_skip,
                          check.names = FALSE, stringsAsFactors = FALSE)
  nm <- names(dt)
  fc_meta <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

  input_gene_name <- NA
  if ("Geneid" %in% nm && !"Chr" %in% nm) {
    rn  <- dt[["Geneid"]]
    mat <- as.matrix(dt[, setdiff(nm, "Geneid"), drop = FALSE])
  } else if (all(fc_meta %in% nm)) {
    rn  <- dt[["Geneid"]]
    mat <- as.matrix(dt[, seq(7L, ncol(dt)), drop = FALSE])
  } else if (all(c("gene_id", "gene_name") %in% nm)) {
    rn  <- dt[["gene_id"]]
    mat <- as.matrix(dt[, setdiff(nm, c("gene_id", "gene_name")), drop = FALSE])
    stripped <- sub("\\..*$", "", rn)
    keep <- !duplicated(stripped)
    input_gene_name <- stats::setNames(dt[["gene_name"]][keep],
                                       stripped[keep])
  } else if ("gene_id" %in% nm) {
    rn  <- dt[["gene_id"]]
    mat <- as.matrix(dt[, setdiff(nm, "gene_id"), drop = FALSE])
  } else {
    stop("Unsupported counts file format: `", path, "` has no `Geneid` or ",
         "`gene_id` column.", call. = FALSE)
  }

  storage.mode(mat) <- "numeric"
  cn <- basename(colnames(mat))
  colnames(mat) <- sub("\\.bam$|\\.sam$|\\.sorted$|\\.markdup$|\\.txt$", "", cn)
  rownames(mat) <- rn
  attr(mat, "input_gene_name") <- input_gene_name
  mat
}

#' Read a sample metadata sheet
#'
#' Reads `.xlsx`/`.xls` through readxl and anything else as delimited text,
#' then renames whichever of `sample_col_candidates` it finds to the canonical
#' `Sample_ID`. Additional required columns are opt-in, so a sheet from another
#' project loads without a project-specific schema error.
#'
#' @param path Path to an Excel or delimited-text metadata file.
#' @param sample_col_candidates Column names tried in order for the sample-ID
#'   column.
#' @param required_cols Character vector of further columns that must be
#'   present; empty by default.
#' @return A data frame with a `Sample_ID` column.
#' @export
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write.csv(data.frame(Sample_ID = c("S1", "S2"), group = c("WT", "KO")),
#'           f, row.names = FALSE)
#' read_metadata(f)
read_metadata <- function(path,
                          sample_col_candidates = c("Sample_ID", "Sample ID"),
                          required_cols = character(0)) {
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    .require_pkg("readxl", "Reading an Excel metadata sheet")
    md <- as.data.frame(readxl::read_excel(path), check.names = FALSE)
  } else {
    md <- utils::read.delim(path, sep = .io_guess_sep(path),
                            check.names = FALSE, stringsAsFactors = FALSE)
  }

  hit <- intersect(sample_col_candidates, names(md))
  if (!length(hit)) {
    stop("Metadata must have a sample-ID column named one of: ",
         paste(sprintf("`%s`", sample_col_candidates), collapse = ", "),
         "; found: ", paste(names(md), collapse = ", "), ".", call. = FALSE)
  }
  names(md)[names(md) == hit[1]] <- "Sample_ID"

  if (length(required_cols)) {
    missing <- setdiff(required_cols, names(md))
    if (length(missing)) {
      stop("Metadata is missing required column(s): ",
           paste(sprintf("`%s`", missing), collapse = ", "), ".", call. = FALSE)
    }
  }
  as.data.frame(md, check.names = FALSE)
}

#' Collapse duplicate row IDs by summing
#'
#' Rows sharing an ID are collapsed with `colSums`. Output rownames come out in
#' `split()` order -- `sort(unique(ids))`, **not** input order -- so callers that
#' carry a parallel annotation table must re-match to `rownames()` afterwards.
#'
#' The `input_gene_name` attribute survives: for each group the first non-`NA`
#' gene name is kept and re-ordered to the collapsed row order.
#'
#' @param mat Numeric matrix.
#' @param ids Character vector of length `nrow(mat)`; defaults to
#'   `rownames(mat)`.
#' @return A numeric matrix with unique rownames.
#' @keywords internal
.aggregate_duplicate_ids <- function(mat, ids = rownames(mat)) {
  stopifnot(is.matrix(mat), length(ids) == nrow(mat))
  ign <- attr(mat, "input_gene_name")

  if (!anyDuplicated(ids)) {
    rownames(mat) <- ids
    attr(mat, "input_gene_name") <- ign
    return(mat)
  }

  groups <- split(seq_len(nrow(mat)), ids)
  collapsed <- do.call(rbind, lapply(groups, function(idx) {
    if (length(idx) == 1L) mat[idx, , drop = FALSE]
    else colSums(mat[idx, , drop = FALSE])
  }))

  if (length(ign) > 1L) {
    new_ign <- vapply(groups, function(idx) {
      vals <- ign[idx]
      vals <- vals[!is.na(vals)]
      if (length(vals)) vals[1L] else NA_character_
    }, character(1L))
    names(new_ign) <- rownames(collapsed)
  } else {
    new_ign <- ign
  }
  attr(collapsed, "input_gene_name") <- new_ign
  collapsed
}

#' Record session provenance to a text file
#'
#' Writes the genome build, the requested Ensembl release, the **resolved**
#' biomaRt archive release, the bulkiRNA and hard-dependency versions, resolved
#' reference-data snapshots, all three [RNGkind()] components, the stochastic
#' defaults declared by [bulkirna_stochastic()], and the full `sessionInfo()`.
#' The resolved archive matters because `sessionInfo()` records only the
#' biomaRt *package* version -- not which remote Ensembl release the annotation
#' actually came from, which is the thing that changes an analysis underneath
#' you.
#'
#' @param path Output file path. Parent directories are created.
#' @param genome_build Character genome build, e.g. `"mm10"`, or `NULL`.
#' @param ensembl_version Ensembl release passed to `useEnsembl()`, or `NULL`
#'   for the floating current release.
#' @return `path`, invisibly.
#' @export
#' @examples
#' f <- tempfile(fileext = ".txt")
#' write_session_provenance(f, genome_build = "mm10")
#' length(readLines(f)) > 0
write_session_provenance <- function(path, genome_build = NULL,
                                     ensembl_version = NULL) {
  ensure_parent_dir(path)
  lines <- paste0("Provenance recorded: ",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  if (!is.null(genome_build)) {
    lines <- c(lines, paste0("Genome build: ", genome_build))
  }
  if (!is.null(ensembl_version)) {
    lines <- c(lines, paste0("Ensembl version (requested): ", ensembl_version))
  }

  rng_kind <- RNGkind()
  stochastic <- .bulkirna_stochastic_registry()
  lines <- c(
    lines,
    paste0(
      "RNGkind(): kind = ", rng_kind[[1L]],
      ", normal.kind = ", rng_kind[[2L]],
      ", sample.kind = ", rng_kind[[3L]]
    ),
    "Stochastic defaults:",
    paste0("  ", stochastic$name, " = ", stochastic$seed_default)
  )

  if (requireNamespace("biomaRt", quietly = TRUE)) {
    arch <- tryCatch({
      archives <- biomaRt::listEnsemblArchives()
      cur <- archives[!is.na(archives$current_release) &
                        archives$current_release %in% c(TRUE, "1", "true"), ]
      if (nrow(cur) > 0) paste0(cur$version[1], " (", cur$date[1], ")")
      else paste0(archives$version[1], " (", archives$date[1], ")")
    }, error = function(e) paste0("unavailable \u2014 ", conditionMessage(e)))
    lines <- c(lines, paste0("Ensembl archive release (resolved): ", arch))
  }

  desc <- suppressWarnings(tryCatch(
    utils::packageDescription("bulkiRNA"),
    error = function(e) e
  ))
  if (inherits(desc, "error") || !is.list(desc)) {
    detail <- if (inherits(desc, "error")) conditionMessage(desc) else
      "package metadata could not be read"
    lines <- c(
      lines,
      paste0("bulkiRNA version: unavailable \u2014 ", detail),
      paste0("Imports: unavailable \u2014 ", detail)
    )
  } else {
    package_version <- desc[["Version"]]
    if (is.null(package_version) || is.na(package_version) ||
        !nzchar(package_version)) {
      package_version <- "unavailable \u2014 Version is missing from package metadata"
    }
    lines <- c(lines, paste0("bulkiRNA version: ", package_version))

    imports_field <- desc[["Imports"]]
    if (is.null(imports_field) || is.na(imports_field) ||
        !nzchar(imports_field)) {
      lines <- c(lines, "Imports: none declared")
    } else {
      imports <- trimws(strsplit(imports_field, ",", fixed = TRUE)[[1L]])
      imports <- sub("[[:space:]]*\\(.*\\)[[:space:]]*$", "", imports)
      imports <- unique(imports[nzchar(imports)])
      import_versions <- vapply(imports, function(package) {
        tryCatch(
          as.character(utils::packageVersion(package)),
          error = function(e) paste0("unavailable \u2014 ", conditionMessage(e))
        )
      }, character(1L))
      lines <- c(lines, "Imports:",
                 paste0("  ", imports, ": ", unname(import_versions)))
    }
  }

  registry_path <- system.file("extdata", "METADATA.yaml",
                               package = "bulkiRNA")
  registry_version <- if (!nzchar(registry_path)) {
    "unavailable \u2014 METADATA.yaml is missing"
  } else {
    registry_lines <- tryCatch(
      readLines(registry_path, warn = FALSE),
      error = function(e) e
    )
    if (inherits(registry_lines, "error")) {
      paste0("unavailable \u2014 ", conditionMessage(registry_lines))
    } else {
      version_line <- grep("^version:[[:space:]]*", registry_lines,
                           value = TRUE)
      if (!length(version_line)) {
        "unavailable \u2014 top-level version key is missing"
      } else {
        value <- sub("^version:[[:space:]]*", "", version_line[[1L]])
        value <- trimws(sub("[[:space:]]+#.*$", "", value))
        first <- substr(value, 1L, 1L)
        last <- substr(value, nchar(value), nchar(value))
        if (nchar(value) >= 2L && first %in% c("\"", "'") && first == last) {
          value <- substr(value, 2L, nchar(value) - 1L)
        }
        if (nzchar(value)) value else
          "unavailable \u2014 top-level version key is empty"
      }
    }
  }
  lines <- c(lines, paste0("Bundled reference registry version: ",
                           registry_version))

  refs <- .ref_resolutions()
  if (!nrow(refs)) {
    lines <- c(lines, "Reference data: none resolved this session")
  } else {
    snapshot <- ifelse(is.na(refs$snapshot), "caller-supplied", refs$snapshot)
    lines <- c(lines, "Reference data:",
               paste0("  ", refs$source, ": ", snapshot,
                      " (", refs$path, ")"))
  }

  lines <- c(lines, "", "--- sessionInfo ---",
             utils::capture.output(utils::sessionInfo()))
  writeLines(lines, con = path)
  invisible(path)
}
