#' Drop msigdbr's cross-collection ortholog cache
#'
#' Works around a correctness bug in msigdbr (confirmed 26.1.0) that silently
#' truncates every collection queried after the first one in a session, but
#' only when ortholog mapping is active (`db_species = "HS"` with a non-human
#' `species`).
#'
#' `msigdbr()` builds its ortholog table from *the genes of the collection
#' currently being queried* --
#' `babelgene::orthologs(genes = unique(mdb$db_ensembl_gene), ...)` -- and then
#' caches it in the package environment under `paste0("orthologs", taxon_id)`.
#' That key does not mention the collection. The next call for a *different*
#' collection finds the key, reuses the previous collection's ortholog table,
#' and `inner_join`s against it -- dropping every gene absent from the first
#' collection. Measured for `species = "Mus musculus"`, `db_species = "HS"`:
#'
#' | call order            | Reactome rows | unique symbols |
#' | --------------------- | ------------- | -------------- |
#' | Reactome alone        | 105806        | 10762          |
#' | Hallmark, do Reactome |  44635        |  3688          |
#'
#' GO:BP is hit harder still (4313 of 15988 symbols survive). The result is not
#' an error but quietly under-tested gene sets, truncated set memberships, and
#' a correspondingly wrong Benjamini-Hochberg family -- so it corrupts `padj`
#' for every collection after the first.
#'
#' Clearing the key before each query makes each call independent, which is the
#' documented contract of `gsdb_msigdb()`. The cost is that `babelgene`
#' recomputes the ortholog mapping per call; that is a few seconds against a
#' warm cache, and it buys correctness. Priming the cache with the union of all
#' collections' genes would avoid the recompute but requires loading the whole
#' database, which is far more expensive than the mapping it saves.
#'
#' This reaches into another package's internals, so it is deliberately *not*
#' the thing that guarantees correctness. If msigdbr renames `pkg_env` or its
#' cache key, this silently becomes a no-op -- and a silent no-op restores
#' exactly the silent truncation it exists to prevent. The guarantee therefore
#' comes from `.msigdbr_assert_ortholog_coverage()`, which checks the *result*
#' and does not care how the mapping was produced. This function is the fix;
#' that one is the seatbelt.
#'
#' @return `TRUE` if a cache entry was dropped, otherwise `FALSE`, invisibly.
#' @keywords internal
#' @noRd
.msigdbr_drop_ortholog_cache <- function() {
  dropped <- tryCatch({
    ns <- asNamespace("msigdbr")
    if (!exists("pkg_env", envir = ns, inherits = FALSE)) return(invisible(FALSE))
    pe <- get("pkg_env", envir = ns, inherits = FALSE)
    if (!is.environment(pe)) return(invisible(FALSE))
    keys <- grep("^orthologs", ls(pe, all.names = TRUE), value = TRUE)
    if (!length(keys)) return(invisible(FALSE))
    rm(list = keys, envir = pe)
    TRUE
  }, error = function(e) FALSE)
  invisible(isTRUE(dropped))
}


#' Assert that ortholog mapping did not silently drop most of a collection
#'
#' The mechanism-independent half of the msigdbr ortholog-cache guard (see
#' `.msigdbr_drop_ortholog_cache()`). Truncation is invisible in the returned
#' table -- the `inner_join` removes whole genes, so nothing looks anomalous
#' from the inside -- but it *is* visible against the source collection, which
#' can be fetched with `species = "Homo sapiens"`. That path never touches the
#' ortholog cache (msigdbr skips the branch entirely for human), so it can
#' neither be poisoned nor poison, and it costs one cached read, no network and
#' no `babelgene` call.
#'
#' The check is a coverage floor: what fraction of the collection's human genes
#' survived mapping to the target species. Measured with the fix in place,
#' `species = "Mus musculus"`, `db_species = "HS"`:
#'
#' | collection        | human genes | mapped | ratio |
#' | ----------------- | ----------- | ------ | ----- |
#' | H                 |  4384       |  4393  | 1.002 |
#' | C2 CP:REACTOME    | 11448       | 10762  | 0.940 |
#' | C2 CP:KEGG_LEGACY |  5246       |  5032  | 0.959 |
#' | C5 GO:MF          | 15971       | 14444  | 0.904 |
#' | C5 GO:CC          | 14953       | 13407  | 0.897 |
#' | C5 GO:BP          | 18215       | 15988  | 0.878 |
#' | C7                | 21385       | 17517  | 0.819 |
#' | C3 TFT:GTRD       | 27022       | 16143  | 0.597 |
#'
#' Truncated by the bug, the same collections land at 0.32 (Reactome) and 0.24
#' (GO:BP). `0.45` sits in the gap: below the worst legitimate collection
#' (0.597, C3's very wide TF-target sets) and above the mildest truncation seen.
#' It is a floor against catastrophic loss, not a quality metric -- a threshold
#' tight enough to grade ortholog quality would fire on legitimate collections.
#'
#' Failure is an error, not a warning. The whole problem with the upstream bug
#' is that it produces plausible numbers, so a warning in a long pipeline log
#' is indistinguishable from success. Anything that gets past this is a
#' complete collection or a stopped script.
#'
#' Skipped when no ortholog mapping happens (human target, or `db_species =
#' "MM"`), and skipped -- with a warning rather than an error -- if the
#' reference fetch itself fails, since a broken cache is not evidence of
#' truncation.
#'
#' @param tbl The tibble returned by [msigdbr::msigdbr()].
#' @param args The argument list that produced `tbl`.
#' @param min_coverage Numeric(1) coverage floor.
#' @return `tbl`, invisibly and unchanged.
#' @keywords internal
#' @noRd
.msigdbr_assert_ortholog_coverage <- function(tbl, args, min_coverage = 0.45) {
  if (!identical(args$db_species, "HS")) return(invisible(tbl))
  if (args$species %in% c("Homo sapiens", "human")) return(invisible(tbl))
  if (!all(c("gene_symbol", "db_gene_symbol") %in% names(tbl))) {
    return(invisible(tbl))
  }

  ref_args <- args
  ref_args$species <- "Homo sapiens"
  ref <- tryCatch(do.call(msigdbr::msigdbr, ref_args), error = function(e) NULL)
  if (is.null(ref) || !"db_gene_symbol" %in% names(ref)) {
    warning("Could not verify ortholog coverage for MSigDB `collection = \"",
            args$collection, "\"`; proceeding unchecked. Gene sets may be ",
            "silently truncated -- see `?gsdb_msigdb`.", call. = FALSE)
    return(invisible(tbl))
  }

  n_source <- length(unique(as.character(ref$db_gene_symbol)))
  n_mapped <- length(unique(as.character(tbl$db_gene_symbol)))
  if (n_source == 0L) return(invisible(tbl))
  coverage <- n_mapped / n_source

  if (coverage < min_coverage) {
    stop("MSigDB ortholog mapping dropped ", format(n_source - n_mapped,
         big.mark = ","), " of ", format(n_source, big.mark = ","),
         " genes for `collection = \"", args$collection, "\"`",
         if (is.null(args$subcollection)) "" else
           paste0(", `subcollection = \"", args$subcollection, "\"`"),
         " (", round(100 * coverage), "% retained, expected >",
         round(100 * min_coverage), "%).\n",
         "This is the msigdbr cross-collection ortholog-cache bug: the gene ",
         "sets would be truncated to another collection's gene space, giving ",
         "under-tested sets and a wrong `padj`. bulkiRNA's workaround did not ",
         "hold, most likely because msigdbr changed its internals.\n",
         "Workaround: query this collection in a fresh R session, before any ",
         "other collection. Please report this so the guard can be updated.",
         call. = FALSE)
  }
  invisible(tbl)
}


#' Load an MSigDB collection
#'
#' Thin provider over [msigdbr::msigdbr()]. MSigDB is downloaded at first use
#' and cached under `$HOME/.cache/R/msigdbr`, so a cold cache needs network
#' access.
#'
#' @param species Character(1) target species, e.g. `"Mus musculus"` or
#'   `"Homo sapiens"`. Human sets are mapped to orthologs when `db_species`
#'   is `"HS"`.
#' @param collection Character(1) MSigDB collection, e.g. `"H"` for the
#'   hallmarks or `"C2"`. See [msigdbr::msigdbr_collections()].
#' @param subcollection Character(1) sub-collection such as `"CP:REACTOME"`,
#'   or `NULL` for the whole collection.
#' @param db_species Character(1) source database, `"HS"` (human MSigDB, the
#'   default, with ortholog mapping) or `"MM"` (mouse-native sets).
#' @param min_size,max_size Integer(1) or `NULL`; drop sets outside these
#'   bounds.
#' @param verbose Logical(1); message what was loaded.
#' @return A [gs_db()] whose `database` key is derived from the collection --
#'   `"msigdb_H"`, `"msigdb_C2_CP_KEGG"` -- with `"MSigDB H"` kept as the
#'   `database_label` attribute for display.
#' @examplesIf requireNamespace("msigdbr", quietly = TRUE) && interactive()
#' db <- gsdb_msigdb("Mus musculus", collection = "H")
#' summary(db)[1:3, ]
#' @export
gsdb_msigdb <- function(species = "Mus musculus",
                        collection = "H",
                        subcollection = NULL,
                        db_species = c("HS", "MM"),
                        min_size = NULL,
                        max_size = NULL,
                        verbose = FALSE) {
  # `msigdbr` is a hard dependency (DESCRIPTION `Imports`), unlike the other
  # optional providers in this package (homologene, yaml, gatom, ...), which
  # are in `Suggests` and genuinely may be missing. It is therefore always
  # installed alongside bulkiRNA and a `requireNamespace()` guard here would
  # be unreachable dead code; call it directly instead.
  species <- .gsdb_species_label(species)
  if (!is.character(collection) || length(collection) != 1L ||
        is.na(collection) || !nzchar(collection)) {
    stop("`collection` must be a single MSigDB collection name such as ",
         "\"H\" or \"C2\"; see `msigdbr::msigdbr_collections()`.",
         call. = FALSE)
  }
  db_species <- match.arg(db_species)

  args <- list(db_species = db_species, species = species,
               collection = collection)
  if (!is.null(subcollection)) {
    args$subcollection <- subcollection
  }
  # Two-part guard against an msigdbr ortholog-cache bug that silently truncates
  # collections: clear the stale cache, then verify the result independently of
  # whether that clear worked. See both helpers above.
  .msigdbr_drop_ortholog_cache()
  tbl <- do.call(msigdbr::msigdbr, args)
  .msigdbr_assert_ortholog_coverage(tbl, args)

  if (!nrow(tbl)) {
    stop("MSigDB returned no gene sets for `collection = \"", collection,
         "\"`", if (is.null(subcollection)) "" else
           paste0(", `subcollection = \"", subcollection, "\"`"),
         " and `species = \"", species, "\"`. Check the names against ",
         "`msigdbr::msigdbr_collections()`.", call. = FALSE)
  }

  sets <- split(as.character(tbl$gene_symbol), as.character(tbl$gs_name))

  # `pathway_names` is the axis/legend label (see `gs_result`), so it is the
  # formatted set name -- "Myc Targets V1" -- not MSigDB's `gs_description`,
  # which is a full sentence ("A subgroup of genes regulated by MYC - version 1
  # (v1).") and unusable as a label. The old toolkit formatted the id at every
  # call site instead; doing it once here keeps every renderer consistent.
  # The description is genuinely useful, so it is kept alongside rather than
  # dropped.
  labels <- stats::setNames(format_pathway_name(names(sets)), names(sets))
  descriptions <- NULL
  if ("gs_description" %in% names(tbl)) {
    uniq <- tbl[!duplicated(tbl$gs_name), c("gs_name", "gs_description")]
    descriptions <- stats::setNames(as.character(uniq$gs_description),
                                    as.character(uniq$gs_name))
  }

  label <- paste(c("MSigDB", collection, subcollection), collapse = " ")
  key <- paste(c("msigdb", collection,
                 if (is.null(subcollection)) NULL
                 else strsplit(subcollection, "[:.]")[[1]]),
               collapse = "_")
  db <- gs_db(sets, database = key, species = species,
              pathway_names = labels, database_label = label,
              pathway_descriptions = descriptions)
  db <- .gs_filter_size(db, min_size, max_size, verbose = verbose)
  if (verbose) {
    message(sprintf("Loaded %s (%s): %d sets.", label, species, length(db)))
  }
  db
}
