#' Registry of bundled reference databases
#'
#' Maps a database key to its `inst/extdata` directory, processed file name
#' and display label. The registry is the fallback when the `yaml` package is
#' not installed; `METADATA.yaml` is the richer source used by [gsdb_info()].
#'
#' @return Named list of registry entries.
#' @keywords internal
.gsdb_registry <- function() {
  list(
    mitopathways = list(
      dir      = "mitocarta3.0",
      rds_file = "mito_mitopathways.rds",
      label    = "MitoPathways 3.0"
    ),
    mitoxplorer = list(
      dir      = "mitoxplorer3.0",
      rds_file = "mito_mitoxplorer.rds",
      label    = "mitoXplorer 3.0"
    ),
    mito_unified = list(
      dir      = "mitochondria_unified",
      rds_file = "unified_mito_pathways.rds",
      label    = "Unified Mitochondrial Pathways"
    ),
    transportdb = list(
      dir      = "transportdb",
      rds_file = "transportdb_genesets.rds",
      label    = "TransportDB 2.0"
    )
  )
}

#' Locate a file inside the installed `extdata` tree
#'
#' The only path-resolution mechanism in the package. Returns `""` when the
#' file is absent, matching [system.file()].
#'
#' @param ... Path components below `inst/extdata`.
#' @return Character(1) absolute path, or `""`.
#' @keywords internal
.gsdb_extdata <- function(...) {
  system.file("extdata", ..., package = "bulkiRNA")
}

#' Read `METADATA.yaml` from the installed reference tree
#'
#' @return The parsed YAML as a list.
#' @keywords internal
.gsdb_metadata <- function() {
  .require_pkg("yaml", "Reading the reference registry")
  path <- .gsdb_extdata("METADATA.yaml")
  if (!nzchar(path)) {
    stop("`METADATA.yaml` is missing from the installed package. ",
         "Reinstall bulkiRNA.", call. = FALSE)
  }
  yaml::read_yaml(path)
}

#' List the bundled gene-set databases
#'
#' Reports every reference database registered in the packaged
#' `METADATA.yaml`, together with the species for which a processed file
#' actually ships.
#'
#' @return A data frame with columns `database`, `name`, `bundled`,
#'   `species` and `description`. `species` is a comma-separated list of the
#'   species directories present in the installed package, or `""` when the
#'   database is not bundled.
#' @examples
#' gsdb_list()
#' @export
gsdb_list <- function() {
  reg <- .gsdb_registry()

  if (!requireNamespace("yaml", quietly = TRUE)) {
    out <- data.frame(
      database    = names(reg),
      name        = unname(vapply(reg, `[[`, character(1L), "label")),
      bundled     = TRUE,
      description = "",
      row.names   = NULL,
      stringsAsFactors = FALSE
    )
    out$species <- unname(vapply(names(reg), .gsdb_species_available,
                                 character(1L)))
    return(out[c("database", "name", "bundled", "species", "description")])
  }

  meta <- .gsdb_metadata()
  dbs <- meta$databases
  out <- data.frame(
    database    = names(dbs),
    name        = vapply(names(dbs),
                         function(k) dbs[[k]]$name %||% k, character(1L)),
    bundled     = vapply(names(dbs),
                         function(k) isTRUE(dbs[[k]]$bundled), logical(1L)),
    description = vapply(names(dbs),
                         function(k) dbs[[k]]$description %||% "",
                         character(1L)),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
  out$species <- vapply(seq_len(nrow(out)), function(i) {
    if (!out$bundled[i]) return("")
    .gsdb_species_available(out$database[i],
                            dbs[[out$database[i]]]$directory %||%
                              out$database[i])
  }, character(1L))
  out[c("database", "name", "bundled", "species", "description")]
}

#' Species with a processed file shipped for a database
#'
#' @param database Character(1) database key.
#' @param dir Character(1) `extdata` sub-directory; looked up in the registry
#'   when `NULL`.
#' @return Comma-separated character(1) of species directory names.
#' @keywords internal
.gsdb_species_available <- function(database, dir = NULL) {
  if (is.null(dir)) {
    dir <- .gsdb_registry()[[database]]$dir %||% database
  }
  proc <- .gsdb_extdata(dir, "processed")
  if (!nzchar(proc)) return("")
  paste(list.dirs(proc, full.names = FALSE, recursive = FALSE),
        collapse = ", ")
}

#' Metadata and citations for one reference database
#'
#' @param database Character(1) database key, e.g. `"mito_unified"`, or a
#'   [gs_db()] object. See [gsdb_list()].
#' @return A list whose first two fields are always `database`, the registry
#'   key, and `name`, the display name. For a registry key, the remaining
#'   fields come from `METADATA.yaml`, with `citations_path` and
#'   `citations_text` added when a `CITATIONS.bib` ships alongside the data.
#'   For a `gs_db`, the remaining fields identify the object and report its set
#'   count, database provenance, and set provenance.
#' @examples
#' info <- gsdb_info("mitopathways")
#' info$name
#' @export
gsdb_info <- function(database) {
  if (is_gs_db(database)) {
    return(list(
      database = attr(database, "database"),
      name = attr(database, "database_label"),
      database_label = attr(database, "database_label"),
      species = attr(database, "species"),
      gene_id_type = attr(database, "gene_id_type"),
      n_sets = length(database),
      provenance = attr(database, "provenance"),
      set_provenance = attr(database, "set_provenance")
    ))
  }
  meta <- .gsdb_metadata()
  if (!is.character(database) || length(database) != 1L) {
    stop("`database` must be a single database name; see `gsdb_list()`.",
         call. = FALSE)
  }
  if (!database %in% names(meta$databases)) {
    stop("`database` must be one of ",
         paste0("\"", names(meta$databases), "\"", collapse = ", "),
         "; got \"", database, "\".", call. = FALSE)
  }
  info <- meta$databases[[database]]
  info <- c(
    list(database = database, name = info$name %||% database),
    info[setdiff(names(info), "name")]
  )
  bib <- .gsdb_extdata(info$directory %||% database, "CITATIONS.bib")
  if (nzchar(bib)) {
    info$citations_path <- bib
    info$citations_text <- readLines(bib, warn = FALSE)
  }
  info
}

#' Load a bundled gene-set database
#'
#' Loads one or more of the processed reference databases that ship with the
#' package and returns them as [gs_db()] objects. Paths resolve only through
#' `system.file("extdata", ..., package = "bulkiRNA")`.
#'
#' @param database Character vector of database keys; see [gsdb_list()]. A
#'   vector of length > 1 returns a named list of `gs_db` objects.
#' @param species Character(1) species, `"Mus musculus"` or
#'   `"Homo sapiens"`. The underscore form is accepted too.
#' @param min_size,max_size Integer(1) or `NULL`; drop sets outside these
#'   bounds. `NULL` (the default) applies no bound.
#' @param rebuild Logical(1). Re-parsing raw source files requires a git
#'   checkout of the toolkit; `TRUE` therefore always errors in an installed
#'   package.
#' @param verbose Logical(1); message what was loaded.
#' @return A `gs_db`, or a named list of `gs_db` when `database` has length
#'   greater than one.
#' @examples
#' db <- gsdb_load("mitopathways")
#' length(db)
#' @export
gsdb_load <- function(database,
                      species = "Mus musculus",
                      min_size = NULL,
                      max_size = NULL,
                      rebuild = FALSE,
                      verbose = FALSE) {
  if (!is.character(database) || length(database) < 1L ||
        anyNA(database) || any(!nzchar(database))) {
    stop("`database` must be a non-empty character vector of database ",
         "names; see `gsdb_list()`.", call. = FALSE)
  }
  if (isTRUE(rebuild)) {
    .gsdb_rebuild_unavailable()
  }

  if (length(database) > 1L) {
    out <- lapply(database, gsdb_load, species = species,
                  min_size = min_size, max_size = max_size,
                  verbose = verbose)
    names(out) <- database
    return(out)
  }

  if (tolower(database) == "gatom") {
    stop("`database = \"gatom\"` names network files that are not bundled ",
         "(too large, ~24 MB) and are ",
         "not gene sets. Fetch them with ",
         "gatom_download_refs(dir = \"00_data/references/gatom\").",
         call. = FALSE)
  }

  reg <- .gsdb_registry()
  if (!database %in% names(reg)) {
    stop("`database` must be one of ",
         paste0("\"", names(reg), "\"", collapse = ", "),
         "; got \"", database, "\". See `gsdb_list()`.", call. = FALSE)
  }
  entry <- reg[[database]]
  sp_dir <- .gsdb_species_dir(species)
  path <- .gsdb_extdata(entry$dir, "processed", sp_dir, entry$rds_file)

  if (!nzchar(path)) {
    have <- .gsdb_species_available(database, entry$dir)
    stop("No processed \"", database, "\" data ships for `species = \"",
         .gsdb_species_label(species), "\"`",
         if (nzchar(have)) paste0("; available: ", have) else "", ".",
         call. = FALSE)
  }

  raw <- readRDS(path)
  if (!is.list(raw) || is.null(raw$T2G)) {
    stop("The processed file for \"", database, "\" is corrupted (no `T2G` ",
         "table): ", path, ". Reinstall bulkiRNA.", call. = FALSE)
  }

  db <- .gsdb_from_t2g(raw, database = database,
                       species = .gsdb_species_label(species),
                       database_label = entry$label)
  db <- .gs_filter_size(db, min_size, max_size, verbose = verbose)
  if (verbose) {
    message(sprintf("Loaded %s (%s): %d sets.", entry$label,
                    .gsdb_species_label(species), length(db)))
  }
  db
}
