#' Registry of stochastic public functions
#'
#' Keeps each stochastic export's seed interface, default, source of
#' randomness, and reason for retaining that default in one place.
#'
#' @return A data frame with one row per stochastic exported function.
#' @keywords internal
.bulkirna_stochastic_registry <- function() {
  data.frame(
    name = c(
      "coresh_match", "coresh_search", "gatom_module",
      "gs_coregulation", "gs_test", "gsdb_coresh"
    ),
    seed_arg = c(
      "seed", "seed", "seed", "seed", NA_character_, "seed"
    ),
    seed_default = c("1L", "1L", "42", "123L", "123L", "1L"),
    source_of_randomness = c(
      "GESECA multilevel permutations",
      "GESECA multilevel permutations",
      "BioNet BUM fit and the MWCS solver heuristic",
      "GESECA multilevel permutations",
      "fgsea multilevel permutations",
      "CoReSh compendium search worker scheduling"
    ),
    note = c(
      "upstream's reference implementation passes the literal value 1L",
      "upstream's reference implementation passes the literal value 1L",
      "the historical default is retained to preserve published results",
      paste0(
        "a new general gene-set verb, so it follows the package's fgsea ",
        "adapter default rather than CoReSh's historical 1L"
      ),
      paste0(
        "the legacy default, accepted through `...` and documented in ",
        "`gs_test_fgsea_params`, is retained"
      ),
      paste0(
        "provider composition pins the compendium search and retains ",
        "CoReSh's historical 1L"
      )
    ),
    stringsAsFactors = FALSE
  )
}

#' List stochastic bulkiRNA functions
#'
#' Reports the public functions whose results consume randomness, the route by
#' which each function accepts a seed, and the default preserved by this
#' package version. `seed_arg` is `NA` when the seed is accepted only through
#' `...`; see the row's `note` for the documented parameter source.
#'
#' The defaults deliberately differ. Existing functions retain historical
#' reproducibility commitments; new functions record their method-family
#' policy rather than imposing a package-wide seed.
#'
#' @param quiet Logical. If `FALSE`, print the registry and return it
#'   invisibly. If `TRUE`, do not print and return it visibly.
#' @return A tibble with columns `name`, `seed_arg`, `seed_default`,
#'   `source_of_randomness`, and `note`.
#' @examples
#' stochastic <- bulkirna_stochastic(quiet = TRUE)
#' stochastic[, c("name", "seed_default")]
#' @export
bulkirna_stochastic <- function(quiet = FALSE) {
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be a single non-missing logical value.", call. = FALSE)
  }

  out <- tibble::as_tibble(.bulkirna_stochastic_registry())

  # n = Inf: a truncated report could hide the function being audited.
  if (!quiet) print(out, n = Inf)
  if (!quiet) return(invisible(out))
  out
}

#' Evaluate an expression with a pinned RNG seed
#'
#' `set.seed()` alone is not enough because `BiocParallel::bplapply()` changes
#' the active generator. This helper is the only place in bulkiRNA allowed to
#' seed stochastic work.
#'
#' @param seed Whole-number RNG seed, or `NULL` to evaluate `expr` without
#'   touching RNG state.
#' @param expr An expression to evaluate exactly once after seeding.
#' @return The value of `expr`.
#' @examples
#' bulkiRNA:::.with_pinned_seed(NULL, 1 + 1)
#' @keywords internal
.with_pinned_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    # RNGkind() re-randomises the seed, so restore the kind before the vector.
    # A legacy-sampler warning describes the caller's choice, not ours.
    suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}
