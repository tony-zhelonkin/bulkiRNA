#' Resolve a human or mouse species label
#'
#' This is the single owner of the package's species aliases. It retains the
#' partial scientific-name matching historically provided by
#' `annotate_genes()` and can optionally return a formatted custom label for
#' `gs_db` providers that support species beyond human and mouse.
#'
#' @param species Character(1) species label.
#' @param allow_custom Logical(1). Permit an unknown non-empty label and
#'   return it with underscores replaced by spaces.
#' @return A list containing the common name, scientific name, three-letter
#'   code, organism annotation package, biomaRt dataset, and GATOM labels.
#'   Fields unavailable for a custom species are `NA_character_`.
#' @keywords internal
.species <- function(species, allow_custom = FALSE) {
  if (!is.logical(allow_custom) || length(allow_custom) != 1L ||
      is.na(allow_custom)) {
    stop("`allow_custom` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.character(species) || length(species) != 1L || is.na(species) ||
      !nzchar(species)) {
    if (allow_custom) {
      stop("`species` must be a single non-empty string.", call. = FALSE)
    }
    stop(
      "`species` must be one of \"human\", \"Homo sapiens\", \"hsa\", ",
      "\"mouse\", \"Mus musculus\", or \"mmu\".",
      call. = FALSE
    )
  }

  key <- gsub("[ _]+", "_", tolower(trimws(species)))
  scientific <- c(human = "Homo sapiens", mouse = "Mus musculus")
  common <- if (key %in% c("human", "homo_sapiens", "hsa")) {
    "human"
  } else if (key %in% c("mouse", "mus_musculus", "mmu")) {
    "mouse"
  } else {
    # match.arg() accepted any unambiguous, case-sensitive prefix of its two
    # scientific choices. Preserve that part of annotate_genes()'s strict
    # contract without rewriting a custom provider label.
    if (!allow_custom) {
      partial <- names(scientific)[startsWith(scientific, species)]
      if (length(partial) == 1L) partial else NA_character_
    } else {
      NA_character_
    }
  }

  if (!is.na(common)) {
    return(switch(
      common,
      human = list(
        common = "human",
        scientific = "Homo sapiens",
        code = "hsa",
        orgdb = "org.Hs.eg.db",
        biomart_dataset = "hsapiens_gene_ensembl",
        gatom_short = "Hs",
        gatom_download = "Homo_sapiens"
      ),
      mouse = list(
        common = "mouse",
        scientific = "Mus musculus",
        code = "mmu",
        orgdb = "org.Mm.eg.db",
        biomart_dataset = "mmusculus_gene_ensembl",
        gatom_short = "Mm",
        gatom_download = "Mus_musculus"
      )
    ))
  }

  if (allow_custom) {
    return(list(
      common = NA_character_,
      scientific = gsub("_", " ", species, fixed = TRUE),
      code = NA_character_,
      orgdb = NA_character_,
      biomart_dataset = NA_character_,
      gatom_short = NA_character_,
      gatom_download = NA_character_
    ))
  }

  stop(
    "`species` must be one of \"human\", \"Homo sapiens\", \"hsa\", ",
    "\"mouse\", \"Mus musculus\", or \"mmu\"; got \"", species, "\".",
    call. = FALSE
  )
}
