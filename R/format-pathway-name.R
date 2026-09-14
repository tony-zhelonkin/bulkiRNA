#' Format pathway names with smart capitalization
#'
#' Turns a database pathway identifier (`HALLMARK_TNFA_SIGNALING_VIA_NFKB`)
#' into a label fit for an axis (`TNF-alpha Signaling via NF-kappaB`), while
#' preserving biological abbreviations, roman numerals and chemical
#' nomenclature. Function words are capitalized when sentence-initial and kept
#' lowercase elsewhere.
#'
#' Underscores **and** dots become spaces -- dots are hierarchy separators in
#' some custom databases (MitoPathways
#' `Metabolism.Lipid_metabolism.Fatty_acid_oxidation`), so each segment reads
#' cleanly.
#'
#' @param text Character vector of pathway names to format.
#' @param use_formatting Logical. `TRUE` (default) applies the smart
#'   capitalization dictionary; `FALSE` returns plain title case.
#' @param strip_prefix Logical. `TRUE` (default) removes common database
#'   prefixes such as `HALLMARK_`, `GOBP_` or `MITOPATHWAYS_`.
#' @return A character vector of formatted names, the same length as `text`.
#' @examples
#' format_pathway_name("HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING")
#' format_pathway_name("GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY")
#' format_pathway_name("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING")
#' @importFrom stringr str_replace str_replace_all str_to_title
#' @export
format_pathway_name <- function(text, use_formatting = TRUE,
                                strip_prefix = TRUE) {
  if (length(text) == 0L || all(is.na(text))) {
    return(text)
  }
  if (!is.character(text)) text <- as.character(text)

  # Underscores and dots to spaces, then collapse runs of spaces.
  text <- stringr::str_replace_all(text, "[_.]", " ")
  text <- stringr::str_replace_all(text, " +", " ")

  if (strip_prefix) {
    for (prefix in .gs_name_prefixes()) {
      text <- stringr::str_replace(text, paste0("^", prefix), "")
    }
  }

  if (!use_formatting) {
    return(stringr::str_to_title(text))
  }

  .gs_smart_capitalization(text)
}

#' Database prefixes stripped from pathway names
#'
#' Matched in their space form (underscores are already spaces) and anchored at
#' the start of the string.
#'
#' @return A character vector of prefix patterns.
#' @keywords internal
.gs_name_prefixes <- function() {
  c(
    # MSigDB collections
    "HALLMARK ", "KEGG ", "REACTOME ", "BIOCARTA ", "MEDICUS ",
    "GOBP ", "GOCC ", "GOMF ", "PID ", "WIKIPATHWAY ",
    "WP ", "^GO ", "GTRD ", "NABA ", "HP ",
    # Custom metabolic / transport databases
    "MITOXPLORER ", "MITOPATHWAYS ", "MITOCARTA ", "TRANSPORTDB ",
    # CoReSh set ids are built as CORESH_<query_name>_<GSE>[_<GPL>] by
    # coresh_sets(), so the prefix is this package's own and always redundant on
    # a figure whose panel is already a CoReSh panel. What remains -- the query
    # name -- is the caller's, and is left alone. Prefer coresh_labels() over
    # displaying the id at all.
    "CORESH "
  )
}

#' Apply smart capitalization with biological exceptions
#'
#' @param text Character vector to format.
#' @return A formatted character vector of the same length.
#' @keywords internal
#' @importFrom stringr str_to_title
.gs_smart_capitalization <- function(text) {
  exceptions <- .gs_name_exceptions()
  multiword <- .gs_name_multiword()
  # These exceptions alone change display form at the start of a label. Same
  # source as the exceptions map, so the two cannot disagree.
  function_words <- names(.gs_function_words())
  greek <- c(
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "iota", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma",
    "tau", "upsilon", "phi", "chi", "psi", "omega"
  )

  vapply(text, function(s) {
    if (is.na(s)) return(NA_character_)
    s_lower <- tolower(s)

    # Protect multi-word terms first, marked with << >>.
    for (info in multiword) {
      s_lower <- gsub(info$pattern, info$replacement, s_lower, fixed = TRUE)
    }

    words <- unlist(strsplit(s_lower, " "))
    formatted <- vapply(seq_along(words), function(i) {
      word <- words[[i]]
      if (nchar(word) == 0L) return("")

      if (grepl("^<<", word)) {
        return(gsub("<<|>>", "", word))
      }

      hit <- exceptions[tolower(names(exceptions)) == word]
      if (length(hit) > 0L) {
        if (i == 1L && word %in% function_words) {
          return(stringr::str_to_title(hit[[1L]]))
        }
        return(hit[[1L]])
      }

      # Roman numerals
      if (grepl("^(i|ii|iii|iv|v|vi|vii|viii|ix|x)$", word)) {
        return(toupper(word))
      }

      # Numbered items: il2 -> IL-2 when the prefix is a known abbreviation
      if (grepl("^([a-z]+)(\\d+[a-z]*)$", word)) {
        prefix <- sub("^([a-z]+)(\\d+[a-z]*)$", "\\1", word)
        number <- sub("^([a-z]+)(\\d+[a-z]*)$", "\\2", word)
        prefix_hit <- exceptions[tolower(names(exceptions)) == prefix]
        if (length(prefix_hit) > 0L) {
          return(paste0(prefix_hit[[1L]], "-", number))
        }
      }

      if (word %in% greek) return(word)

      # Chemical prefixes: n-, o-, cis-, trans- ...
      if (grepl("^(n|o|p|m|s|cis|trans|alpha|beta|gamma|delta)-", word)) {
        parts <- unlist(strsplit(word, "-"))
        return(paste(parts[1L], stringr::str_to_title(parts[-1L]),
                     sep = "-", collapse = "-"))
      }

      # Ordinals
      if (grepl("^\\d+(st|nd|rd|th)$", word)) return(word)

      stringr::str_to_title(word)
    }, character(1L), USE.NAMES = FALSE)

    gsub("<<|>>", "", paste(formatted, collapse = " "))
  }, character(1L), USE.NAMES = FALSE)
}

#' Multi-word biological terms protected before word-by-word formatting
#'
#' Order matters: longest and most specific first.
#'
#' @return A list of `pattern`/`replacement` pairs.
#' @keywords internal
.gs_name_multiword <- function() {
  list(
    # NF-kappaB family
    list(pattern = "nf kappa b", replacement = "<<NF-kappaB>>"),
    list(pattern = "nf-kappa-b", replacement = "<<NF-kappaB>>"),
    list(pattern = "nfkappab", replacement = "<<NF-kappaB>>"),
    list(pattern = "nfkb", replacement = "<<NF-kappaB>>"),

    # TGF-beta
    list(pattern = "tgf beta", replacement = "<<TGF-beta>>"),
    list(pattern = "tgfbeta", replacement = "<<TGF-beta>>"),

    # IFN variants
    list(pattern = "ifn alpha", replacement = "<<IFN-alpha>>"),
    list(pattern = "ifnalpha", replacement = "<<IFN-alpha>>"),
    list(pattern = "ifn beta", replacement = "<<IFN-beta>>"),
    list(pattern = "ifnbeta", replacement = "<<IFN-beta>>"),
    list(pattern = "ifn gamma", replacement = "<<IFN-gamma>>"),
    list(pattern = "ifngamma", replacement = "<<IFN-gamma>>"),

    # Interferon spelled out
    list(pattern = "interferon alpha", replacement = "<<Interferon-alpha>>"),
    list(pattern = "interferon beta", replacement = "<<Interferon-beta>>"),
    list(pattern = "interferon gamma", replacement = "<<Interferon-gamma>>"),

    # TNF-alpha
    list(pattern = "tnf alpha", replacement = "<<TNF-alpha>>"),
    list(pattern = "tnfa", replacement = "<<TNF-alpha>>"),
    list(pattern = "tumor necrosis factor alpha",
         replacement = "<<Tumor Necrosis Factor-alpha>>"),

    # HIF-1alpha
    list(pattern = "hif 1 alpha", replacement = "<<HIF-1alpha>>"),
    list(pattern = "hif1alpha", replacement = "<<HIF-1alpha>>"),

    # IL / STAT combinations
    list(pattern = "il2 stat5", replacement = "<<IL-2/STAT5>>"),
    list(pattern = "il6 jak stat3", replacement = "<<IL-6/JAK/STAT3>>"),

    # Common multi-word process terms
    list(pattern = "cell cycle", replacement = "<<Cell Cycle>>"),
    list(pattern = "cell death", replacement = "<<Cell Death>>"),
    list(pattern = "apoptotic process", replacement = "<<Apoptotic Process>>"),
    list(pattern = "immune response", replacement = "<<Immune Response>>"),
    list(pattern = "inflammatory response",
         replacement = "<<Inflammatory Response>>"),
    list(pattern = "innate immune", replacement = "<<Innate Immune>>"),
    list(pattern = "adaptive immune", replacement = "<<Adaptive Immune>>"),
    list(pattern = "antigen processing", replacement = "<<Antigen Processing>>"),
    list(pattern = "antigen presentation",
         replacement = "<<Antigen Presentation>>")
  )
}

#' Biological abbreviation dictionary
#'
#' Keys are lowercase for case-insensitive matching; values are the display
#' form.
#'
#' @return A named list.
#' @keywords internal
.gs_name_exceptions <- function() {
  c(
    list(
    # NF-kappaB family
    "nf" = "NF", "nfkb" = "NF-kappaB", "nfkappab" = "NF-kappaB",

    # Interleukins
    "il" = "IL", "il1" = "IL-1", "il2" = "IL-2", "il4" = "IL-4",
    "il6" = "IL-6", "il8" = "IL-8", "il10" = "IL-10", "il12" = "IL-12",
    "il17" = "IL-17",

    # Interferons
    "ifn" = "IFN", "ifna" = "IFN-alpha", "ifnb" = "IFN-beta",
    "ifng" = "IFN-gamma", "interferon" = "Interferon",

    # Tumor necrosis factors
    "tnf" = "TNF", "tnfa" = "TNF-alpha", "tnfr" = "TNFR",
    "tnfr1" = "TNFR1", "tnfr2" = "TNFR2",

    # Immunology
    "mhc" = "MHC", "hla" = "HLA", "tcr" = "TCR", "bcr" = "BCR",
    "tlr" = "TLR", "tlr4" = "TLR4", "cd" = "CD", "cd4" = "CD4",
    "cd8" = "CD8", "cd28" = "CD28", "cd80" = "CD80", "cd86" = "CD86",

    # Signaling molecules
    "stat" = "STAT", "stat1" = "STAT1", "stat3" = "STAT3",
    "stat5" = "STAT5", "jak" = "JAK", "mapk" = "MAPK", "erk" = "ERK",
    "jnk" = "JNK", "pi3k" = "PI3K", "akt" = "AKT", "mtor" = "mTOR",
    "ampk" = "AMPK",

    # Transcription factors
    "ap" = "AP", "ap1" = "AP-1", "creb" = "CREB", "irf" = "IRF",
    "irf1" = "IRF1", "irf3" = "IRF3", "irf7" = "IRF7",

    # Cell cycle and nucleic acids
    "cdk" = "CDK", "apc" = "APC", "dna" = "DNA", "rna" = "RNA",
    "mrna" = "mRNA", "rrna" = "rRNA", "trna" = "tRNA", "mirna" = "miRNA",

    # Metabolism
    "atp" = "ATP", "adp" = "ADP", "amp" = "AMP", "nadh" = "NADH",
    "nadph" = "NADPH", "fadh2" = "FADH2", "coa" = "CoA",
    "acetyl" = "Acetyl", "tca" = "TCA", "oxphos" = "OXPHOS",

    # Lipids
    "ldl" = "LDL", "hdl" = "HDL", "vldl" = "VLDL", "lps" = "LPS",
    "pge2" = "PGE2",

    # Reactive species
    "ros" = "ROS", "no" = "NO", "nos" = "NOS",

    # Enzymes
    "cox" = "COX", "cox2" = "COX-2", "lox" = "LOX", "lt" = "LT",
    "ex" = "EX", "pla2" = "PLA2",

    # Amino acids
    "gln" = "Gln", "glu" = "Glu", "gls" = "GLS", "glul" = "GLUL",
    "asn" = "Asn", "asp" = "Asp", "arg" = "Arg",

    # Growth factors
    "egf" = "EGF", "fgf" = "FGF", "pdgf" = "PDGF", "vegf" = "VEGF",
    "tgf" = "TGF", "tgfb" = "TGF-beta", "igf" = "IGF", "ngf" = "NGF",
    "gdnf" = "GDNF", "bdnf" = "BDNF",

    # Receptor types
    "gpcr" = "GPCR", "rtk" = "RTK",

    # Membrane-transporter families (TransportDB classification codes)
    "abc" = "ABC", "mfs" = "MFS", "mc" = "MC", "vic" = "VIC",
    "dmt" = "DMT", "mpt" = "MPT", "p-atpase" = "P-ATPase",
    "f-atpase" = "F-ATPase",

    # Organelles
    "er" = "ER", "ecm" = "ECM",

    # Misc
    "uv" = "UV", "hif" = "HIF", "hif1a" = "HIF-1alpha", "p53" = "p53",
    "rb" = "Rb", "brca1" = "BRCA1", "brca2" = "BRCA2", "parp" = "PARP",
    "sumo" = "SUMO", "ubiquitin" = "Ubiquitin",

    # Processes
    "autophagy" = "Autophagy", "apoptosis" = "Apoptosis",
    "necrosis" = "Necrosis", "ferroptosis" = "Ferroptosis",
    "pyroptosis" = "Pyroptosis"
    ),
    # Function words come from .gs_function_words() rather than being listed
    # here, because .gs_smart_capitalization() needs the same set to decide what
    # to capitalise sentence-initially. Two hand-kept copies of one list drift
    # the moment somebody adds "with". as.list() matters: this is a list(), so
    # splicing the character vector directly would bury all nine words in one
    # unnamed element and every lookup for them would miss.
    as.list(.gs_function_words())
  )
}

#' Function words kept lowercase mid-label
#'
#' The single source for both the exceptions map and the sentence-initial rule in
#' `.gs_smart_capitalization()`.
#'
#' @return Named character vector mapping each function word to itself.
#' @keywords internal
.gs_function_words <- function() {
  words <- c("via", "and", "or", "of", "in", "to", "by", "from", "the")
  stats::setNames(words, words)
}

#' Wrap a label onto two or three lines
#'
#' Word-aware wrapping that splits a long label in half, or into thirds when it
#' is very long. Ported from the old `smart_wrap()`.
#'
#' @param text Character vector of labels.
#' @param width Soft character width. Labels longer than `width` split in two;
#'   longer than `1.5 * width` split in three.
#' @return A character vector with embedded newlines.
#' @keywords internal
.gs_wrap_label <- function(text, width = 50) {
  vapply(text, function(s) {
    if (is.na(s)) return(NA_character_)
    words <- unlist(strsplit(s, " "))
    if (length(words) <= 1L) return(s)

    total <- nchar(s)
    if (total > width * 1.5) {
      third <- ceiling(length(words) / 3)
      two_thirds <- min(third * 2, length(words))
      parts <- list(
        words[seq_len(third)],
        words[seq(third + 1L, two_thirds)],
        if (two_thirds < length(words)) {
          words[seq(two_thirds + 1L, length(words))]
        } else {
          character(0)
        }
      )
      parts <- Filter(function(p) length(p) > 0L, parts)
      return(paste(vapply(parts, paste, character(1L), collapse = " "),
                   collapse = "\n"))
    }
    if (total > width) {
      mid <- ceiling(length(words) / 2)
      return(paste(
        paste(words[seq_len(mid)], collapse = " "),
        paste(words[seq(mid + 1L, length(words))], collapse = " "),
        sep = "\n"
      ))
    }
    s
  }, character(1L), USE.NAMES = FALSE)
}
