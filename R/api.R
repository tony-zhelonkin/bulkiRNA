#' Public API stability registry
#'
#' `bulkirna_api()` is the machine-readable contract for the package's public
#' surface. Stability and the historical signature freeze are separate axes:
#'
#' * `stable` functions follow semantic versioning. Incompatible changes need
#'   a major release, apart from removals that completed a documented
#'   deprecation cycle.
#' * `experimental` functions may change or be removed without a major release.
#'   CoReSh dataset-taking functions use the `(gse, gpl)` pair as the dataset
#'   key; a GSE accession alone is not unique in that compendium.
#' * `deprecated` functions, when present, remain callable and warn until the
#'   version recorded in `removed_in`.
#'
#' `frozen` records the remaining exported signatures carried forward from the
#' script-library API. Most stable functions were introduced after the freeze
#' and therefore have `frozen = FALSE`.
#'
#' `stochastic` records whether calling the function consumes randomness, so
#' its result depends on a seed. [bulkirna_stochastic()] reports the seed
#' interface, default, and source of randomness for each such function.
#'
#' The `superseded_by` and `removed_in` fields retain a stable registry schema
#' even when the current public API has no deprecated entries.
#' `layer` records the function's functional module. Stable and experimental
#' functions retain their curated registry group.
#'
#' @param lifecycle `"all"`, or one or more of `"stable"`, `"experimental"`,
#'   and `"deprecated"`.
#' @param quiet Logical. If `FALSE`, print the registry and return it
#'   invisibly. If `TRUE`, do not print and return it visibly.
#' @return A tibble with one row per selected exported function and columns `name`,
#'   `layer`, `lifecycle`, `frozen`, `stochastic`, `superseded_by`, and
#'   `removed_in`.
#' @examples
#' experimental <- bulkirna_api("experimental", quiet = TRUE)
#' @export
bulkirna_api <- function(lifecycle = "all", quiet = FALSE) {
  lifecycle_requested <- match.arg(
    lifecycle,
    choices = c("all", "stable", "experimental", "deprecated"),
    several.ok = TRUE
  )
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be a single non-missing logical value.", call. = FALSE)
  }

  out <- .bulkirna_api_registry()
  if (!"all" %in% lifecycle_requested) {
    out <- out[out$lifecycle %in% lifecycle_requested, , drop = FALSE]
  }
  out <- tibble::as_tibble(out)

  # n = Inf: a truncated report is worse than none, since the rows it hides
  # are exactly the ones someone is checking for.
  if (!quiet) print(out, n = Inf)
  if (!quiet) return(invisible(out))
  out
}

#' Registry of public API metadata
#'
#' Keeps the layer, lifecycle, signature-freeze, stochasticity, and deprecation
#' metadata used by [bulkirna_api()] in one place.
#'
#' @return A data frame with one row per exported function.
#' @keywords internal
.bulkirna_api_registry <- function() {
  stable <- list(
    gsdb = c(
      "gsdb_from_file", "gsdb_info", "gsdb_list", "gsdb_load",
      "gsdb_msigdb", "gsdb_register"
    ),
    gs = c(
      "gs_filter", "gs_leading_edge", "gs_palette", "gs_plot_bar",
      "gs_plot_dot", "gs_plot_heatmap", "gs_plot_running", "gs_plot_size",
      "gs_ranks", "gs_read", "gs_master_columns", "gs_save",
      "gs_scale_fonts", "gs_score", "gs_split", "gs_stat_types", "gs_test",
      "gs_to_master", "gs_top", "gs_validate_master", "gs_write"
    ),
    de = c(
      "de_bfc_plot", "de_md_plot", "de_pca", "de_pca_3d", "de_volcano",
      "de_volcano_grid"
    ),
    gatom = c(
      "gatom_de", "gatom_download_refs", "gatom_genes", "gatom_module",
      "gatom_refs", "gatom_save_html"
    ),
    `top-level` = c(
      "annotate_genes", "build_dge", "bulki_palettes", "bulkirna_api",
      "bulkirna_check_deps", "bulkirna_stochastic", "ensure_dir",
      "format_pathway_name", "read_counts_matrix", "read_metadata",
      "theme_bulki", "write_session_provenance"
    )
  )

  experimental <- list(
    gsdb = "gsdb_coresh",
    gs = "gs_coregulation",
    # All dataset-level CoReSh APIs share the (gse, gpl) composite key.
    coresh = c(
      "coresh_chunks", "coresh_convergence", "coresh_labels",
      "coresh_loadings", "coresh_match", "coresh_search", "coresh_sets",
      "coresh_validate"
    ),
    `top-level` = c(
      "entrez_to_gene", "filter_confounder_genes", "gene_to_entrez"
    )
  )

  deprecated <- character(0L)

  stable_names <- unlist(stable, use.names = FALSE)
  experimental_names <- unlist(experimental, use.names = FALSE)
  names_all <- c(stable_names, experimental_names, deprecated)

  deprecated_layer <- character(0L)
  layer <- c(
    rep(names(stable), lengths(stable)),
    rep(names(experimental), lengths(experimental)),
    deprecated_layer
  )
  lifecycle <- c(
    rep("stable", length(stable_names)),
    rep("experimental", length(experimental_names)),
    rep("deprecated", length(deprecated))
  )

  frozen_names <- c("build_dge", "ensure_dir", "format_pathway_name")

  superseded_by <- rep(NA_character_, length(names_all))
  names(superseded_by) <- names_all

  removed_in <- rep(NA_character_, length(names_all))

  stochastic_names <- .bulkirna_stochastic_registry()$name

  out <- data.frame(
    name = names_all,
    layer = unname(layer),
    lifecycle = lifecycle,
    frozen = names_all %in% frozen_names,
    stochastic = names_all %in% stochastic_names,
    superseded_by = unname(superseded_by),
    removed_in = removed_in,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$name), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read a successor from a legacy fixture
#'
#' The first expression in every legacy fixture is its `.Deprecated()` call.
#' For fixtures that use `.Deprecated(new =)` (including an unnamed first
#' argument), evaluate that argument. Message-only fixtures return
#' `NA_character_`.
#'
#' @param name Name of a legacy fixture.
#' @return A character scalar or `NA_character_`.
#' @keywords internal
.bulkirna_deprecation_target <- function(name) {
  namespace <- environment(.bulkirna_deprecation_target)
  fun <- get(name, envir = namespace, mode = "function", inherits = FALSE)
  expr <- body(fun)
  if (!is.call(expr) || !identical(expr[[1L]], as.name("{")) ||
      length(expr) < 2L) {
    stop("Deprecated shim `", name, "` has an unexpected body.", call. = FALSE)
  }

  call <- expr[[2L]]
  if (!is.call(call) || !identical(call[[1L]], as.name(".Deprecated"))) {
    stop("Deprecated shim `", name, "` must call .Deprecated() first.",
         call. = FALSE)
  }

  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (!is.null(arg_names) && "msg" %in% arg_names &&
      !"new" %in% arg_names) {
    return(NA_character_)
  }

  new <- if (!is.null(arg_names) && "new" %in% arg_names) {
    args[[which(arg_names == "new")[[1L]]]]
  } else {
    args[[1L]]
  }
  value <- eval(new, envir = environment(fun))
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("Deprecated shim `", name, "` has an invalid successor.", call. = FALSE)
  }
  value
}
