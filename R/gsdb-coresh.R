#' Build a CoReSh-derived gene-set database
#'
#' Searches the CoReSh compendium by percentage of variance explained, keeps
#' the requested number of datasets per query, and turns their strongest gene
#' loadings into a [gs_db()]. This is a provider-layer composition of
#' [coresh_search()] and [coresh_sets()]; it does not score, extract, map, or
#' de-duplicate genes itself.
#'
#' @param queries A non-empty named list of integer Entrez vectors, each with
#'   at least three unique IDs.
#' @param chunk_dir Optional explicit path to an `hsa` or `mmu` chunk directory.
#'   When `NULL`, resolve the corresponding directory through the shared
#'   refcache.
#' @param species A human or mouse alias accepted by [.species()].
#' @param top_hits Positive whole number of top CoReSh datasets to retain per
#'   query.
#' @param top_n Positive whole number of absolute loadings to retain per hit.
#' @param min_size Positive whole-number minimum set size.
#' @param max_size Positive whole-number maximum set size.
#' @param jaccard_threshold Numeric threshold in `[0, 1]` passed to
#'   [coresh_sets()].
#' @param n_cores Positive whole number passed to [coresh_search()].
#' @param seed Whole-number RNG seed pinned around the compendium search. The
#'   default matches the upstream CoReSh reference implementation.
#' @return A [gs_db()] with `database = "coresh"`. Its database provenance
#'   records the resolved snapshot, species, query names and unique Entrez-ID
#'   counts, and every set-building parameter. It does not store the exact
#'   Entrez IDs; reproducing the database requires the caller's own query
#'   definitions. Set provenance identifies the query, dataset, chunk, loading
#'   cutoff, and CoReSh rank for each retained set.
#' @examples
#' \dontrun{
#' db <- gsdb_coresh(
#'   list(iron_uptake = c(7037L, 4891L, 55240L)),
#'   species = "human",
#'   n_cores = 1L
#' )
#' }
#' @export
gsdb_coresh <- function(queries, chunk_dir = NULL, species = "human",
                         top_hits = 5L, top_n = 50L,
                         min_size = 15L, max_size = 500L,
                         jaccard_threshold = 0.8, n_cores = 4L,
                         seed = 1L) {
  top_hits <- .coresh_positive_integer(top_hits, "top_hits")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("`seed` must be one finite whole number.", call. = FALSE)
  }

  hits <- .with_pinned_seed(
    seed,
    coresh_search(
      queries = queries,
      chunk_dir = chunk_dir,
      species = species,
      n_cores = n_cores,
      pvalues = FALSE,
      seed = seed
    )
  )
  selected <- hits[hits$rank <= top_hits, , drop = FALSE]

  coresh_sets(
    selected,
    queries = queries,
    chunk_dir = chunk_dir,
    species = species,
    top_n = top_n,
    min_size = min_size,
    max_size = max_size,
    jaccard_threshold = jaccard_threshold
  )
}
