#' Raw reference files are not shipped
#'
#' Rebuilding a reference database re-parses the raw source files under
#' `data/references/`, which live in the git checkout and are
#' `.Rbuildignore`d. Any code path that needs them fails through here with an
#' explicit message instead of a missing-file error.
#'
#' @param what Character(1) describing the caller, used in the message.
#' @return Never returns; always throws.
#' @keywords internal
.gsdb_rebuild_unavailable <- function(what = "`rebuild = TRUE`") {
  stop(what, " needs a source checkout of bulkiRNA; raw reference files are ",
       "not shipped with the installed package. Clone the repository and run ",
       "the rebuild from `data/references/`.", call. = FALSE)
}

#' Rebuild the bundled processed reference files from raw sources
#'
#' The internal successor to `scripts/GSEA/GSEA_processing/build_reference_databases.R`,
#' which executed at source time. This is a plain function: it parses the raw
#' files under `refs_dir`, applies the size filter and Jaccard
#' de-duplication, and writes the processed RDS files in the legacy
#' `list(T2G, T2N)` shape that [gsdb_load()] reads. It is a maintainer tool,
#' not part of the user surface, and errors via
#' [.gsdb_rebuild_unavailable()] when the raw tree is absent.
#'
#' @param species Character(1), `"Mus musculus"` or `"Homo sapiens"`.
#' @param refs_dir Character(1) path to the `data/references` directory of a
#'   source checkout.
#' @param min_size,max_size Integer(1) size bounds applied to every database.
#' @param write Logical(1); write the processed RDS files under `refs_dir`.
#' @return Named list of `gs_db` objects, invisibly.
#' @keywords internal
.gsdb_rebuild <- function(species = "Mus musculus",
                          refs_dir = file.path("data", "references"),
                          min_size = 5L,
                          max_size = 500L,
                          write = TRUE) {
  species <- .gsdb_species_label(species)
  sp_dir <- .gsdb_species_dir(species)
  if (!dir.exists(refs_dir)) {
    .gsdb_rebuild_unavailable("Rebuilding reference databases")
  }

  out <- list()

  gmx <- file.path(refs_dir, "mitocarta3.0", "raw", "MitoPathways3.0.gmx")
  if (file.exists(gmx)) {
    db <- gsdb_from_file(gmx, database = "mitopathways",
                         species = species, prefix = "MITOPATHWAYS",
                         database_label = "MitoPathways 3.0")
    if (grepl("musculus", species, ignore.case = TRUE)) {
      db <- .gsdb_human_to_mouse(db)
    }
    out$mitopathways <- .gs_filter_size(db, min_size, max_size)
  }

  mx <- file.path(refs_dir, "mitoxplorer3.0", "raw", "mouse_gene_function.txt")
  if (file.exists(mx) && !grepl("sapiens", species, ignore.case = TRUE)) {
    out$mitoxplorer <- .gs_filter_size(
      .gsdb_parse_mitoxplorer(mx, species = species), min_size, max_size
    )
  }

  if (!is.null(out$mitopathways)) {
    # `!duplicated(names(merged))` is first-wins, not a union of the two
    # sets' genes (the old builder unioned via `split()` on rbind'ed T2G
    # frames). This is safe today only because the two inputs are
    # prefix-disjoint ("MITOPATHWAYS_..." vs "MITOXPLORER_..." set ids) and
    # can never collide; revisit if a third, non-prefixed component is added
    # to this merge.
    merged <- c(unclass(out$mitopathways), unclass(out$mitoxplorer))
    labels <- c(attr(out$mitopathways, "pathway_names"),
                attr(out$mitoxplorer, "pathway_names"))
    unified <- gs_db(merged[!duplicated(names(merged))],
                     database = "mito_unified",
                     species = species,
                     pathway_names = labels,
                     database_label = "Unified Mitochondrial Pathways")
    out$mito_unified <- .gsdb_dedup_sets(unified, threshold = 0.99)
  }

  tdb <- file.path(refs_dir, "transportdb", "raw", "TransportDB2.0.csv")
  if (file.exists(tdb)) {
    out$transportdb <- .gs_filter_size(
      .gsdb_parse_transportdb(tdb, species = species), min_size, max_size
    )
  }

  if (write) {
    targets <- list(
      mitopathways = c("mitocarta3.0", "mito_mitopathways.rds"),
      mitoxplorer  = c("mitoxplorer3.0", "mito_mitoxplorer.rds"),
      mito_unified = c("mitochondria_unified", "unified_mito_pathways.rds"),
      transportdb  = c("transportdb", "transportdb_genesets.rds")
    )
    for (key in names(out)) {
      tgt <- targets[[key]]
      path <- file.path(refs_dir, tgt[1], "processed", sp_dir, tgt[2])
      ensure_parent_dir(path)
      legacy <- .gsdb_as_t2g(out[[key]])
      # Deliberately the display label ("MitoPathways 3.0"), not the old
      # builder's machine string ("MitoPathways3.0" / "TransportDB2.0"). No
      # code in this package reads `source`; a caller that string-matched
      # the old provenance value will need updating.
      legacy$source <- attr(out[[key]], "database_label")
      legacy$created <- Sys.time()
      saveRDS(legacy, path)
      message(sprintf("Wrote %s (%d sets).", path, length(out[[key]])))
    }
  }

  invisible(out)
}

#' Parse the mitoXplorer gene-function table
#'
#' @param path Character(1) path to `mouse_gene_function.txt`.
#' @param species Character(1) species label.
#' @param prefix Character(1) set-id prefix.
#' @param gene_col Character(1) column holding gene symbols.
#' @param process_col Character(1) column holding the mitochondrial process.
#' @return A [gs_db()].
#' @keywords internal
.gsdb_parse_mitoxplorer <- function(path,
                                    species = "Mus musculus",
                                    prefix = "MITOXPLORER",
                                    gene_col = "MGI_symbol",
                                    process_col = "mito_process") {
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  d <- utils::read.delim(path, stringsAsFactors = FALSE, header = TRUE)
  for (col in c(gene_col, process_col)) {
    if (!col %in% names(d)) {
      stop("Column `", col, "` is missing from ", basename(path),
           "; found ", paste0("`", names(d), "`", collapse = ", "), ".",
           call. = FALSE)
    }
  }
  by_proc <- split(as.character(d[[gene_col]]), as.character(d[[process_col]]))
  by_proc <- by_proc[nzchar(names(by_proc))]
  ids <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9]", "_",
                                          names(by_proc))))
  gs_db(stats::setNames(by_proc, ids),
        database = "mitoxplorer", species = species,
        pathway_names = stats::setNames(names(by_proc), ids),
        database_label = "mitoXplorer 3.0")
}

#' Parse the TransportDB CSV
#'
#' Gene sets are the transporter families. Column names vary between
#' TransportDB exports, so the gene, family and description columns are
#' detected from a list of known aliases.
#'
#' @param path Character(1) path to the TransportDB CSV.
#' @param species Character(1) species label.
#' @param prefix Character(1) set-id prefix.
#' @return A [gs_db()].
#' @keywords internal
.gsdb_parse_transportdb <- function(path,
                                    species = "Mus musculus",
                                    prefix = "TRANSPORTDB") {
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  pick <- function(aliases, label) {
    hit <- intersect(aliases, names(d))
    if (!length(hit)) {
      stop("Could not find the ", label, " column in ", basename(path),
           "; found ", paste0("`", names(d), "`", collapse = ", "), ".",
           call. = FALSE)
    }
    hit[1]
  }
  gene_col <- pick(c("Symbol", "Gene_Symbol", "gene_symbol", "SYMBOL",
                     "Gene"), "gene symbol")
  fam_col <- pick(c("Family", "family", "TC_Family"), "transporter family")
  desc_hit <- intersect(c("Description", "description",
                          "Family_Description"), names(d))

  fams <- unique(d[[fam_col]])
  fams <- fams[!is.na(fams) & nzchar(fams)]
  ids <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9]", "_", fams)))
  sets <- lapply(fams, function(f) as.character(d[[gene_col]][d[[fam_col]] == f]))
  labels <- vapply(fams, function(f) {
    if (!length(desc_hit)) return(f)
    val <- unique(d[[desc_hit[1]]][d[[fam_col]] == f])[1]
    if (is.na(val) || !nzchar(val)) f else val
  }, character(1L), USE.NAMES = FALSE)

  gs_db(stats::setNames(sets, ids),
        database = "transportdb", species = species,
        pathway_names = stats::setNames(labels, ids),
        database_label = "TransportDB 2.0")
}

#' Drop near-duplicate gene sets by Jaccard similarity
#'
#' When two sets overlap at or above `threshold`, one is dropped; a
#' `MITOPATHWAYS_` set wins over any other, matching the merge rule the
#' shipped unified mitochondrial database was built with.
#'
#' @param db A [gs_db()].
#' @param threshold Numeric(1) Jaccard cutoff in `[0, 1]`.
#' @param verbose Logical(1); message how many sets were dropped.
#' @return A `gs_db` with redundant sets removed.
#' @keywords internal
.gsdb_dedup_sets <- function(db, threshold = 0.99, verbose = FALSE) {
  ids <- names(db)
  n <- length(ids)
  if (n <= 1L) return(db)

  drop <- character()
  for (i in seq_len(n - 1L)) {
    if (ids[i] %in% drop) next
    for (j in seq(i + 1L, n)) {
      if (ids[j] %in% drop) next
      a <- db[[ids[i]]]
      b <- db[[ids[j]]]
      jac <- length(intersect(a, b)) / length(union(a, b))
      if (jac >= threshold) {
        drop <- c(drop, if (grepl("^MITOPATHWAYS_", ids[j])) ids[i] else ids[j])
      }
    }
  }
  if (verbose && length(drop)) {
    message(sprintf("Dropped %d redundant sets (Jaccard >= %.2f).",
                    length(drop), threshold))
  }
  db[setdiff(ids, drop)]
}

#' Map human gene symbols to mouse orthologs
#'
#' @param db A [gs_db()] holding human symbols.
#' @param verbose Logical(1); message mapping statistics.
#' @return A `gs_db` with mouse symbols and `species = "Mus musculus"`.
#' @keywords internal
.gsdb_human_to_mouse <- function(db, verbose = FALSE) {
  .require_pkg("homologene", "Converting human symbols to mouse orthologs")
  human <- unique(unlist(db, use.names = FALSE))
  map <- homologene::homologene(human, inTax = 9606, outTax = 10090)
  if (!nrow(map)) {
    stop("`db` symbols produced no human-to-mouse orthologs; supply a gs_db ",
         "whose genes are human symbols.", call. = FALSE)
  }
  map <- map[!duplicated(map[[1]]), ]
  lookup <- stats::setNames(as.character(map[[2]]), as.character(map[[1]]))
  if (verbose) {
    message(sprintf("Mapped %d / %d human symbols to mouse orthologs.",
                    sum(human %in% names(lookup)), length(human)))
  }
  sets <- lapply(unclass(db), function(g) {
    m <- unname(lookup[g])
    m[!is.na(m)]
  })
  retained <- names(sets)[vapply(sets, length, integer(1L)) > 0L]
  set_provenance <- attr(db, "set_provenance")
  if (!is.null(set_provenance)) {
    set_provenance <- set_provenance[
      match(retained, set_provenance$set_name),
      ,
      drop = FALSE
    ]
  }
  gs_db(sets, database = attr(db, "database"), species = "Mus musculus",
        pathway_names = attr(db, "pathway_names"),
        database_label = attr(db, "database_label"),
        pathway_descriptions = attr(db, "pathway_descriptions"),
        set_provenance = set_provenance,
        provenance = attr(db, "provenance"))
}
