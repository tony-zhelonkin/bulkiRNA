#' Download GATOM reference network files
#'
#' GATOM's atom-transition metabolic networks are too large to bundle
#' (~24 MB), so they are fetched on demand from the Artyomov Lab server. This
#' is a downloader, not a gene-set provider: the files it writes are consumed
#' by GATOM itself, and it does not return a [gs_db()].
#'
#' Files already present are skipped unless `overwrite = TRUE`. A file that
#' fails to download warns and is left out of the return value, so a partial
#' run is visible rather than silent.
#'
#' @param dir Character(1) destination directory; created if missing.
#' @param species Character(1) human or mouse species alias.
#' @param networks Character vector of networks to fetch, e.g.
#'   `c("kegg", "combined")`.
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
  if (!is.character(networks) || !length(networks) || anyNA(networks) ||
        any(!nzchar(networks))) {
    stop("`networks` must be a non-empty character vector, e.g. ",
         "c(\"kegg\", \"combined\").", call. = FALSE)
  }

  files <- unique(c(
    unlist(lapply(networks, function(net) {
      c(sprintf("network.%s.rds", net),
        sprintf("met.%s.db.rds", net),
        sprintf("gene2reaction.%s.%s.eg.tsv", net, sp$code))
    })),
    sprintf("org.%s.eg.gatom.anno.rds", sp$gatom_short)
  ))

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

  message(sprintf("Downloaded %d / %d files to: %s",
                  length(downloaded), length(files), dir))
  invisible(downloaded)
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
