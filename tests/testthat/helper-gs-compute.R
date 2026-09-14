# Builders for the compute-layer tests. B1 owns the real `gsdb_*` providers;
# these build the frozen `gs_db` shape (CONVENTIONS.md 11a) with structure().

fake_gs_db <- function(sets = NULL, database = "testdb",
                       species = "Mus musculus") {
  sets <- sets %||% list(
    SET_UP   = paste0("G", 1:20),
    SET_DOWN = paste0("G", 81:100),
    SET_MID  = paste0("G", 40:60)
  )
  structure(
    sets,
    pathway_names = stats::setNames(
      paste("Pretty", names(sets)), names(sets)
    ),
    database = database,
    species = species,
    gene_id_type = "symbol",
    class = "gs_db"
  )
}

# A rank vector over G1..G100 where G1..G20 are high and G81..G100 are low.
fake_ranks <- function(n = 100L) {
  stats::setNames(seq(n, 1) / 10 - n / 20, paste0("G", seq_len(n)))
}

fake_expr <- function(n_gene = 100L, n_samp = 8L, seed = 1L) {
  set.seed(seed)
  m <- matrix(
    stats::rnorm(n_gene * n_samp), nrow = n_gene,
    dimnames = list(paste0("G", seq_len(n_gene)),
                    paste0("s", seq_len(n_samp)))
  )
  # plant a group effect in SET_UP
  m[1:20, 1:4] <- m[1:20, 1:4] + 2
  m
}

fake_sample_data <- function(n_samp = 8L) {
  data.frame(
    group = rep(c("KO", "WT"), each = n_samp / 2),
    row.names = paste0("s", seq_len(n_samp)),
    stringsAsFactors = FALSE
  )
}
