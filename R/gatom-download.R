#' GATOM reference file names, by network
#'
#' The one table of upstream file names, taken from the index at
#' `artyomovlab.wustl.edu/publications/supp_materials/GATOM/`. Both
#' [gatom_download_refs()] and [gatom_refs()] read it, so the downloader and
#' the loader cannot disagree about a name. The vignette's prose names the
#' lipid metabolite file `met.rhea.lipids.db.rds`; the server serves
#' `met.lipids.db.rds`, which is also what the vignette's code reads.
#'
#' @param network Character(1), one of `"kegg"`, `"rhea"`, `"combined"`,
#'   `"lipids"`.
#' @param species A species record from `.species()`.
#' @return A named character vector: `network`, `met_db`, `org_anno` and,
#'   for every network except KEGG, `gene2reaction_extra`.
#' @keywords internal
.gatom_ref_files <- function(network, species) {
  files <- switch(
    network,
    kegg = c(network = "network.kegg.rds", met_db = "met.kegg.db.rds"),
    rhea = c(network = "network.rhea.rds", met_db = "met.rhea.db.rds",
             gene2reaction_extra = "gene2reaction.rhea.%s.eg.tsv"),
    combined = c(network = "network.combined.rds",
                 met_db = "met.combined.db.rds",
                 gene2reaction_extra = "gene2reaction.combined.%s.eg.tsv"),
    lipids = c(network = "network.rhea.lipids.rds",
               met_db = "met.lipids.db.rds",
               gene2reaction_extra = "gene2reaction.rhea.%s.eg.tsv")
  )
  if (is.null(files)) {
    stop("`network` must be one of ",
         paste0("\"", .gatom_networks, "\"", collapse = ", "), "; got \"",
         network, "\".", call. = FALSE)
  }
  if ("gene2reaction_extra" %in% names(files)) {
    files[["gene2reaction_extra"]] <- sprintf(files[["gene2reaction_extra"]],
                                              species$code)
  }
  c(files, org_anno = sprintf("org.%s.eg.gatom.anno.rds",
                              species$gatom_short))
}

.gatom_networks <- c("kegg", "rhea", "combined", "lipids")

.gatom_check_network <- function(network, arg = "network", several = FALSE) {
  ok <- is.character(network) && length(network) >= 1L && !anyNA(network) &&
    (several || length(network) == 1L) && all(network %in% .gatom_networks)
  if (!ok) {
    stop("`", arg, "` must be ", if (several) "a subset of " else "one of ",
         paste0("\"", .gatom_networks, "\"", collapse = ", "), "; got ",
         paste0("\"", paste(network, collapse = "\", \""), "\""), ".",
         call. = FALSE)
  }
  invisible(network)
}

#' Download GATOM reference network files
#'
#' GATOM's atom-transition metabolic networks are too large to bundle
#' (~24 MB), so they are fetched on demand from the Artyomov Lab server. This
#' is a downloader, not a gene-set provider: the files it writes are consumed
#' by GATOM itself through [gatom_refs()], and it does not return a [gs_db()].
#'
#' Each network fetches its network file, its metabolite database and, for
#' Rhea, combined and lipids, its `gene2reaction` supplement; the species
#' annotation is fetched once. Files already present are skipped unless
#' `overwrite = TRUE`. A file that fails to download warns and is left out of
#' the return value, so a partial run is visible rather than silent.
#'
#' The directory's `SHA256SUMS` is rewritten to cover every file present in
#' it, so [gatom_refs()] can verify what it loads later.
#'
#' @param dir Character(1) destination directory; created if missing.
#' @param species Character(1) human or mouse species alias.
#' @param networks Character vector of networks to fetch: any of `"kegg"`,
#'   `"rhea"`, `"combined"`, `"lipids"`.
#' @param overwrite Logical(1); re-download files that already exist.
#' @return Character vector of downloaded (or already present) file paths,
#'   invisibly.
#' @examples
#' \dontrun{
#' gatom_download_refs(dir = "00_data/references/gatom")
#' }
#' @export
gatom_download_refs <- function(
    dir = "00_data/references/gatom",
    species = "Mus_musculus",
    networks = c("kegg", "combined"),
    overwrite = FALSE
) {
  base_url <-
    "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM"

  sp <- .species(species)
  .gatom_check_network(networks, "networks", several = TRUE)

  files <- unique(unlist(lapply(networks, function(net) {
    unname(.gatom_ref_files(net, sp))
  })))

  ensure_dir(dir)

  downloaded <- character()
  for (fname in files) {
    dest_file <- file.path(dir, fname)

    # An interrupted transfer used to leave a partial file in place that every
    # later run reported as `[skip] ... (exists)`, so `gatom_refs()` then died
    # in `readRDS()` naming neither the file nor the fix. A zero-byte file is
    # treated as absent, and the transfer itself goes to `<dest>.part` and is
    # renamed only on success, so an interruption cannot produce one.
    if (file.exists(dest_file) && !isTRUE(overwrite)) {
      if (isTRUE(file.info(dest_file)$size > 0)) {
        message(sprintf("  [skip] %s (exists; use overwrite = TRUE to replace)",
                        fname))
        downloaded <- c(downloaded, dest_file)
        next
      }
      message(sprintf("  [redo] %s (present but empty -- refetching)", fname))
    }

    message(sprintf("  Downloading %s ...", fname))
    part_file <- paste0(dest_file, ".part")
    ok <- tryCatch({
      utils::download.file(file.path(base_url, fname), part_file,
                           mode = "wb", quiet = TRUE)
      if (!isTRUE(file.info(part_file)$size > 0)) {
        stop("the download produced an empty file", call. = FALSE)
      }
      if (!file.rename(part_file, dest_file)) {
        stop("could not move the download into place at ", dest_file,
             call. = FALSE)
      }
      TRUE
    }, error = function(e) {
      unlink(part_file)
      warning(sprintf("  [FAIL] %s: %s", fname, conditionMessage(e)),
              call. = FALSE)
      FALSE
    })
    if (ok) {
      downloaded <- c(downloaded, dest_file)
      message(sprintf("  [ok] %s (%.1f MB)", fname,
                      file.info(dest_file)$size / 1e6))
    }
  }

  .gatom_write_sums(dir)
  message(sprintf("Downloaded %d / %d files to: %s",
                  length(downloaded), length(files), dir))
  invisible(downloaded)
}

#' Write `SHA256SUMS` for a GATOM reference directory
#'
#' Covers every reference file in `dir`, in the format `sha256sum -c` reads.
#'
#' @param dir Character(1) reference directory.
#' @return The path to `SHA256SUMS`, invisibly; `NA` when nothing was
#'   written.
#' @keywords internal
.gatom_write_sums <- function(dir) {
  present <- list.files(dir, pattern = "\\.(rds|tsv)$")
  if (!length(present)) return(invisible(NA_character_))
  sums <- .gatom_sha256(file.path(dir, present))
  if (anyNA(sums)) {
    message("  SHA256SUMS not written: tools::sha256sum() needs R >= 4.5.")
    return(invisible(NA_character_))
  }
  path <- file.path(dir, "SHA256SUMS")
  writeLines(paste0(unname(sums), "  ", present), path)
  invisible(path)
}

#' SHA-256 of files, where this R can compute it
#'
#' `tools::sha256sum()` arrived in R 4.5.0 and the package supports R 4.2.
#'
#' @param paths Character vector of file paths.
#' @return A character vector of hashes named by path, `NA` throughout on
#'   R < 4.5.
#' @keywords internal
.gatom_sha256 <- function(paths) {
  tools_ns <- asNamespace("tools")
  if (!exists("sha256sum", envir = tools_ns, inherits = FALSE)) {
    return(stats::setNames(rep(NA_character_, length(paths)), paths))
  }
  get("sha256sum", envir = tools_ns)(paths)
}

#' Deprecated GATOM reference downloader name
#'
#' `download_gatom_references()` is the frozen pre-package name. Use
#' [gatom_download_refs()] for the layer-prefixed API and its consistent `dir`
#' argument.
#'
#' @param dest_dir Character(1) destination directory; created if missing.
#' @param species Character(1) human or mouse species alias.
#' @param networks Character vector of networks to fetch.
#' @param overwrite Logical(1); re-download files that already exist.
#' @return Character vector returned invisibly by [gatom_download_refs()].
#' @examples
#' \dontrun{
#' download_gatom_references(dest_dir = "00_data/references/gatom")
#' }
#' @keywords internal
download_gatom_references <- function(
    dest_dir = "00_data/references/gatom",
    species = "Mus_musculus",
    networks = c("kegg", "combined"),
    overwrite = FALSE
) {
  .Deprecated(new = "gatom_download_refs")
  gatom_download_refs(
    dir = dest_dir,
    species = species,
    networks = networks,
    overwrite = overwrite
  )
}
