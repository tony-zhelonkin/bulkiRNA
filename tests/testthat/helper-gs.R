# Small builders shared by the unit tests. Not a fixture — the committed
# fixture in tests/fixtures/ belongs to the golden harness and is off limits.

fake_gs_df <- function(n = 3L, stat = NULL) {
  data.frame(
    pathway_id = paste0("SET_", seq_len(n)),
    pathway_name = paste("Set", seq_len(n)),
    n_genes = seq_len(n) * 10L,
    n_genes_tested = seq_len(n) * 8L,
    stat = stat %||% seq(-1, 1, length.out = n),
    p_value = seq(0.001, 0.2, length.out = n),
    padj = seq(0.01, 0.4, length.out = n),
    stringsAsFactors = FALSE
  )
}

fake_gs_result <- function(n = 3L, contrast = "KO-WT", ...) {
  bulkiRNA:::gs_result(
    fake_gs_df(n, ...),
    database = "testdb", contrast = contrast,
    method = "fgsea", stat_type = "NES"
  )
}

fake_gs_matrix <- function(n_path = 3L, n_samp = 4L) {
  m <- matrix(
    seq_len(n_path * n_samp) / 10,
    nrow = n_path,
    dimnames = list(paste0("SET_", seq_len(n_path)),
                    paste0("s", seq_len(n_samp)))
  )
  bulkiRNA:::gs_matrix(
    m,
    database = "testdb", method = "gsva",
    sample_data = data.frame(
      group = rep(c("WT", "KO"), length.out = n_samp),
      row.names = colnames(m)
    )
  )
}

fake_coregulation_input <- function(n_genes = 600L, n_samples = 16L) {
  set.seed(8041)
  genes <- paste0("G", seq_len(n_genes))
  expr <- matrix(
    stats::rnorm(n_genes * n_samples),
    nrow = n_genes,
    dimnames = list(genes, paste0("S", seq_len(n_samples)))
  )
  latent <- as.numeric(scale(sin(seq(0, 3 * pi, length.out = n_samples))))
  coregulated <- genes[seq_len(25L)]
  expr[coregulated, ] <- matrix(
    latent, nrow = length(coregulated), ncol = n_samples, byrow = TRUE
  ) + matrix(stats::rnorm(length(coregulated) * n_samples, sd = 0.15),
             nrow = length(coregulated))

  scrambled <- genes[seq.int(1L, n_genes, length.out = 25L)]
  random_sets <- stats::setNames(
    lapply(seq_len(14L), function(i) sample(genes, 25L)),
    paste0("RANDOM_", seq_len(14L))
  )
  sets <- c(
    list(COREGULATED = coregulated, SCRAMBLED = scrambled),
    random_sets
  )
  list(
    expr = expr,
    db = gsdb_register(
      sets,
      database = "coreg_demo",
      species = "Homo sapiens",
      pathway_names = stats::setNames(
        paste("Pretty", names(sets)), names(sets)
      )
    ),
    coregulated = coregulated,
    scrambled = scrambled
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# `tests/fixtures/` is in .Rbuildignore, so it ships with the source tree but not
# with the built package. Tests that need it run under devtools::test() and the
# golden harness, and skip under R CMD check rather than erroring there.
coresh_micro_fixture <- function() {
  path <- testthat::test_path("..", "fixtures", "coresh-chunk-micro.rds")
  testthat::skip_if_not(file.exists(path),
                        "coresh-chunk-micro.rds is not in the built package")
  readRDS(path)
}

# Restore the process-wide RNG state after a test that deliberately mutates it.
local_pinned_rng <- function(.local_envir = parent.frame()) {
  original_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  original_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  withr::defer({
    suppressWarnings(RNGkind(
      original_kind[1L], original_kind[2L], original_kind[3L]
    ))
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, envir = .local_envir)
  invisible(NULL)
}
