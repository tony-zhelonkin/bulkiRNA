master_result <- function(stat_type = "NES", leading_edge = TRUE) {
  x <- data.frame(
    pathway_id = c("SET_A", "SET_B"),
    pathway_name = c("SET A raw", "SET B raw"),
    n_genes = c(4L, 3L),
    n_genes_tested = c(3L, 2L),
    stat = c(2.1, -1.4),
    p_value = c(1e-20, 0.03),
    padj = c(1e-20, 0.04),
    stringsAsFactors = FALSE
  )
  if (leading_edge) {
    x$leading_edge <- list(c("A", "B"), "D")
  }
  bulkiRNA:::gs_result(
    x,
    database = "testdb",
    contrast = "treated-control",
    method = "fgsea",
    stat_type = stat_type
  )
}

master_db <- function() {
  list(
    SET_A = c("A", "B", "C", "outside"),
    SET_B = c("D", "E", "outside")
  )
}

test_that("master schema v1 has the contracted rows and order", {
  path <- system.file(
    "extdata", "master-schema-v1.csv", package = "bulkiRNA"
  )
  schema <- utils::read.csv(path, stringsAsFactors = FALSE)

  expect_identical(
    schema$column,
    c(
      "pathway_id", "pathway_name", "database", "contrast", "nes",
      "pvalue", "padj", "set_size", "leading_edge_size", "gene_ratio",
      "core_enrichment", "genes_full_set", "direction", "neg_log_padj",
      "entity_type"
    )
  )
  expect_identical(
    schema$type,
    c(
      rep("character", 4L), "numeric", "numeric", "numeric", "integer",
      "integer", "numeric", "character", "character", "character",
      "numeric", "character"
    )
  )
  expect_identical(schema$required, c(rep("yes", 14L), "no"))
  expect_true(all(nzchar(schema$description)))
})

test_that("gs_master_columns matches what gs_to_master actually emits", {
  without <- gs_to_master(master_result())
  with_entity_type <- gs_to_master(master_result(), entity_type = "pathway")

  # Each mode is compared against a real table rather than against the
  # documented order, so a drifting docstring cannot pass this.
  expect_identical(gs_master_columns(), names(without))
  expect_identical(gs_master_columns(optional = TRUE), names(with_entity_type))

  # Selecting by the canonical order is the idiom the accessor exists for, so
  # assert it works rather than inferring it from the column names.
  expect_identical(names(without[, gs_master_columns()]), names(without))
  expect_identical(
    names(with_entity_type[, gs_master_columns(optional = TRUE)]),
    names(with_entity_type)
  )

  expect_identical(
    setdiff(gs_master_columns(optional = TRUE), gs_master_columns()),
    "entity_type"
  )
})

test_that("gs_master_columns rejects a non-logical optional", {
  expect_error(gs_master_columns(optional = NA), "must be TRUE or FALSE")
  expect_error(gs_master_columns(optional = "yes"), "must be TRUE or FALSE")
})

test_that("gs_to_master creates a valid master table", {
  master <- gs_to_master(
    master_result(),
    db = master_db(),
    universe = c("A", "B", "C", "D", "E")
  )

  expect_s3_class(master, "tbl_df")
  expect_identical(attr(master, "schema_version"), "1")
  expect_identical(
    names(master),
    c(
      "pathway_id", "pathway_name", "database", "contrast", "nes",
      "pvalue", "padj", "set_size", "leading_edge_size", "gene_ratio",
      "core_enrichment", "genes_full_set", "direction", "neg_log_padj"
    )
  )
  expect_identical(master$direction, c("Up", "Down"))
  expect_identical(master$leading_edge_size, c(2L, 1L))
  expect_equal(master$gene_ratio, c(2 / 3, 1 / 2))
  expect_identical(master$core_enrichment, c("A/B", "D"))
  expect_identical(master$genes_full_set, c("A/B/C", "D/E"))
  expect_equal(master$neg_log_padj, -log10(master$padj))
  expect_silent(gs_validate_master(master))
  expect_visible(gs_validate_master(master, error = FALSE))
  expect_equal(nrow(gs_validate_master(master, error = FALSE)), 0L)
})

test_that("a partial NA fill in neg_log_padj is rejected", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  master <- master[c(1L, 2L, 1L), ]
  master$pathway_id <- c("MSIGDB_SET_A", "MSIGDB_SET_B", "CORESH_SET_A")
  master$neg_log_padj[3L] <- NA_real_

  expect_error(
    gs_validate_master(master),
    "neg_log_padj.*1 row.*CORESH_SET_A"
  )
  problems <- gs_validate_master(master, error = FALSE)
  problem <- problems[problems$check == "neg_log_padj_missing", ]
  expect_equal(nrow(problem), 1L)
  expect_identical(problem$column, "neg_log_padj")
  expect_match(problem$message, "1 row")
  expect_match(problem$message, "CORESH_SET_A")
})

test_that("the retired cap-at-16 neg_log_padj convention is rejected", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  master$neg_log_padj <- pmin(master$neg_log_padj, 16)

  expect_error(gs_validate_master(master), "cap-at-16")
})

test_that("reordered master columns are rejected", {
  master <- gs_to_master(master_result())
  master <- master[c(2L, 1L, seq.int(3L, ncol(master)))]

  expect_error(gs_validate_master(master), "Column order")
})

test_that("non-NES statistics require a deliberate override", {
  res <- master_result(stat_type = "t")

  expect_error(gs_to_master(res), "got.*t")
  expect_error(gs_to_master(res), "stat_as_nes = TRUE", fixed = TRUE)

  master <- gs_to_master(res, stat_as_nes = TRUE)
  expect_identical(master$nes, res$stat)
  expect_silent(gs_validate_master(master))
})

test_that("db NULL produces deliberately unavailable full-set genes", {
  master <- gs_to_master(master_result(), db = NULL)

  expect_type(master$genes_full_set, "character")
  expect_true(all(is.na(master$genes_full_set)))
  expect_silent(gs_validate_master(master))
})

test_that("entity_type is emitted first", {
  master <- gs_to_master(master_result(), entity_type = "pathway")

  expect_identical(names(master)[1L], "entity_type")
  expect_identical(master$entity_type, rep("pathway", 2L))
  expect_silent(gs_validate_master(master))
})

test_that("a result without leading_edge still forms a valid master", {
  master <- gs_to_master(master_result(leading_edge = FALSE),
                         db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))

  expect_true(all(is.na(master$leading_edge_size)))
  expect_true(all(is.na(master$gene_ratio)))
  expect_true(all(is.na(master$core_enrichment)))
  expect_silent(gs_validate_master(master))
})

test_that("leading_edge_size follows non-empty core_enrichment", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  master$leading_edge_size[2L] <- NA_integer_

  problems <- gs_validate_master(master, error = FALSE)
  problem <- problems[problems$check == "leading_edge_size_missing", ]
  expect_equal(nrow(problem), 1L)
  expect_identical(problem$column, "leading_edge_size")
  expect_match(problem$message, "1 row")
  expect_match(problem$message, "SET_B")
})

test_that("an NA pathway_id is rejected as a missing identity", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  master$pathway_id[2L] <- NA_character_

  problems <- gs_validate_master(master, error = FALSE)
  problem <- problems[problems$check == "identity_missing" &
                        problems$column == "pathway_id", ]

  expect_equal(nrow(problem), 1L)
  expect_match(problem$message, "1 row")
  expect_match(problem$message, "<NA>", fixed = TRUE)
  expect_error(gs_validate_master(master), "pathway_id.*1 row")
})

test_that("gene_ratio is required only when its inputs are present", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  master$gene_ratio[1L] <- NA_real_

  problems <- gs_validate_master(master, error = FALSE)
  problem <- problems[problems$check == "gene_ratio_missing", ]
  expect_equal(nrow(problem), 1L)
  expect_identical(problem$column, "gene_ratio")
  expect_match(problem$message, "1 row")
  expect_match(problem$message, "SET_A")
  expect_error(gs_validate_master(master), "gene_ratio.*1 row.*SET_A")

  master$leading_edge_size[1L] <- NA_integer_
  master$core_enrichment[1L] <- NA_character_
  expect_silent(gs_validate_master(master))
})

test_that("validation reports every problem and accepts safe coercion", {
  master <- gs_to_master(master_result(), db = master_db(),
                         universe = c("A", "B", "C", "D", "E"))
  character_master <- master
  character_master$set_size <- as.character(character_master$set_size)
  character_master$padj <- as.character(character_master$padj)
  expect_silent(gs_validate_master(character_master))

  master$direction[1L] <- "up"
  master$neg_log_padj <- pmin(master$neg_log_padj, 16)
  problems <- gs_validate_master(master, error = FALSE)

  expect_true(all(c("direction_values", "neg_log_padj_values") %in%
                  problems$check))
  expect_error(gs_validate_master(master), "direction")
  expect_error(gs_validate_master(master), "cap-at-16")
})
