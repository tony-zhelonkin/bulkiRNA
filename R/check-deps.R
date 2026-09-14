#' Registry of optional dependencies
#'
#' Keeps the feature and repository metadata used by [bulkirna_check_deps()]
#' in one place. `babelgene` supports `msigdbr` ortholog mapping and the
#' corresponding cache regression tests.
#'
#' @return A data frame with one row per non-development `Suggests` package.
#' @keywords internal
.bulkirna_optional_deps <- function() {
  package <- c(
    "edgeR", "limma",
    "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
    # msigdbr uses babelgene internally for ortholog mapping;
    # test-gsdb-msigdb.R skips on it.
    "babelgene",
    "biomaRt", "homologene",
    "GSVA",
    "gatom", "mwcsr", "igraph",
    # qs2 reads the CoReSh chunk tree; BiocParallel spreads the search over
    # chunk files. Both are reached only through `coresh_*()`.
    "qs2", "BiocParallel",
    # patchwork left this registry at 1.1.0: gs_plot_running() composes three
    # real panels with it, so it is a hard Import rather than optional.
    "plotly",
    "readxl", "yaml"
  )
  feature <- c(
    rep("de", 2L),
    rep("annotation", 6L),
    "scoring",
    rep("network", 3L),
    rep("coresh", 2L),
    "plots",
    rep("io", 2L)
  )
  repository <- c(
    "Bioconductor", "Bioconductor",
    "Bioconductor", "Bioconductor", "Bioconductor", "CRAN",
    "Bioconductor", "CRAN",
    "Bioconductor",
    "Bioconductor", "CRAN", "CRAN",
    "CRAN", "Bioconductor",
    "CRAN",
    "CRAN", "CRAN"
  )
  install <- ifelse(
    repository == "Bioconductor",
    sprintf('BiocManager::install("%s")', package),
    sprintf('install.packages("%s")', package)
  )

  data.frame(
    package = package,
    feature = feature,
    repository = repository,
    install = install,
    stringsAsFactors = FALSE
  )
}

#' Check package installation state
#'
#' @param package Character vector of package names.
#' @return A data frame with logical `installed` and character `version`
#'   columns, in the same order as `package`.
#' @keywords internal
.bulkirna_dependency_status <- function(package) {
  installed <- vapply(
    package,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1L)
  )
  version <- vapply(seq_along(package), function(i) {
    if (!installed[[i]]) return(NA_character_)
    as.character(utils::packageVersion(package[[i]]))
  }, character(1L))

  data.frame(
    installed = unname(installed),
    version = unname(version),
    stringsAsFactors = FALSE
  )
}

#' Check optional bulkiRNA dependencies
#'
#' Reports whether the optional packages needed by one or more feature areas
#' are installed, along with their versions and source-appropriate install
#' commands. This is a convenience preflight; feature functions still check
#' their own dependencies when called.
#'
#' @param features `"all"`, or one or more of `"de"`, `"annotation"`,
#'   `"scoring"`, `"network"`, `"coresh"`, `"plots"`, and `"io"`.
#' @param quiet Logical. If `FALSE`, print the report and return it invisibly.
#'   If `TRUE`, do not print and return it visibly.
#' @param error Logical. Stop with a non-zero script exit when any requested
#'   optional package is missing.
#' @return A tibble with columns `package`, `feature`, `repository`,
#'   `installed`, `version`, and `install`.
#' @examples
#' deps <- bulkirna_check_deps(c("de", "scoring"), quiet = TRUE)
#' deps[, c("package", "installed")]
#' @export
bulkirna_check_deps <- function(features = "all", quiet = FALSE,
                                error = FALSE) {
  choices <- c("all", "de", "annotation", "scoring", "network", "coresh",
               "plots", "io")
  features <- match.arg(features, choices = choices, several.ok = TRUE)

  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be a single non-missing logical value.", call. = FALSE)
  }
  if (!is.logical(error) || length(error) != 1L || is.na(error)) {
    stop("`error` must be a single non-missing logical value.", call. = FALSE)
  }

  deps <- .bulkirna_optional_deps()
  if (!"all" %in% features) {
    deps <- deps[deps$feature %in% features, , drop = FALSE]
  }
  status <- .bulkirna_dependency_status(deps$package)
  deps$installed <- status$installed
  deps$version <- status$version
  deps <- deps[c("package", "feature", "repository", "installed", "version",
                 "install")]
  out <- tibble::as_tibble(deps)

  # n = Inf: a truncated report is worse than none, since the rows it hides
  # are exactly the ones someone is checking for.
  if (!quiet) print(out, n = Inf)

  missing <- !out$installed
  if (error && any(missing)) {
    commands <- paste(unique(out$install[missing]), collapse = "; ")
    stop(
      "Missing optional package(s) selected by `features`: ",
      paste(out$package[missing], collapse = ", "),
      ". Install with: ", commands, ".",
      call. = FALSE
    )
  }

  if (!quiet) return(invisible(out))
  out
}
