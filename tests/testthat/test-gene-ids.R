test_that("gene_to_entrez() maps human symbols in order and drops failures", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  expect_warning(
    mapped <- gene_to_entrez(c("EGFR", "NOT_A_REAL_GENE", "TP53", "EGFR")),
    "1/4 symbols failed to map: NOT_A_REAL_GENE"
  )
  expect_identical(mapped, c(1956L, 7157L, 1956L))
  expect_null(names(mapped))
})

test_that("gene_to_entrez() accepts all documented species aliases", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")

  human <- lapply(c("human", "Homo sapiens", "hsa"), function(species) {
    gene_to_entrez("TP53", species = species)
  })
  mouse <- lapply(c("mouse", "Mus musculus", "mmu"), function(species) {
    gene_to_entrez("Trp53", species = species)
  })

  expect_true(all(vapply(human, identical, logical(1L), 7157L)))
  expect_true(all(vapply(mouse, identical, logical(1L), 22059L)))
})

test_that("gene_to_entrez() names at most 20 unmapped symbols", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  warning_text <- NULL
  mapped <- withCallingHandlers(
    gene_to_entrez(c("TP53", paste0("NOT_A_REAL_GENE_", seq_len(21L)))),
    warning = function(w) {
      warning_text <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(mapped, 7157L)
  expect_match(warning_text, "21/22 symbols failed to map")
  expect_match(warning_text, "NOT_A_REAL_GENE_20")
  expect_false(grepl("NOT_A_REAL_GENE_21", warning_text, fixed = TRUE))
})

test_that("entrez_to_gene() maps in order, names output, and drops failures", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  expect_warning(
    mapped <- entrez_to_gene(c(1956L, 999999999L, 7157L, 1956L)),
    "1/4 Entrez IDs failed to map: 999999999"
  )
  expect_identical(
    mapped,
    c(`1956` = "EGFR", `7157` = "TP53", `1956` = "EGFR")
  )
})

test_that("entrez_to_gene() maps mouse identifiers", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Mm.eg.db")

  expect_identical(
    entrez_to_gene(c(22059L, 13649L), species = "mmu"),
    c(`22059` = "Trp53", `13649` = "Egfr")
  )
})

test_that("gene identifier converters validate dependency-free inputs", {
  expect_error(gene_to_entrez("TP53", species = "rat"),
               "`species` must be one of")
  expect_error(gene_to_entrez(c("TP53", NA_character_)),
               "`symbols` must be a character vector")
  expect_error(gene_to_entrez("TP53", multi_vals = "list"),
               "`multi_vals` must be one of")
  expect_error(entrez_to_gene(c(7157L, NA_integer_)),
               "`entrez` must be an integer")
  expect_error(entrez_to_gene(7157.5), "must be whole Entrez identifiers")
  expect_identical(gene_to_entrez(character()), integer())
  expect_identical(entrez_to_gene(integer()), character())
})

test_that("each human confounder category filters only its matches", {
  cases <- list(
    ribosomal = list(c("RPS3", "KEEP", "RPL10"), 2L),
    mito = list(c("MT-ND1", "KEEP", "MT-CO1"), 2L),
    hemoglobin = list(c("HBA1", "KEEP", "HBB", "HBG2"), 3L),
    cell_cycle = list(
      c("MKI67", "TOP2A", "CCNB1", "KEEP", "CCNB2", "CDK1", "CDC20",
        "BIRC5", "TYMS"),
      8L
    ),
    sex = list(c("XIST", "KEEP", "RPS4Y1", "DDX3Y", "KDM5D"), 4L)
  )

  for (category in names(cases)) {
    case <- cases[[category]]
    expect_message(
      filtered <- filter_confounder_genes(
        case[[1L]], drop = category
      ),
      paste0(category, " removed ", case[[2L]], " gene\\(s\\)\\."),
      info = category
    )
    expect_identical(filtered, "KEEP", info = category)
  }
})

test_that("each confounder category handles mouse symbol capitalization", {
  cases <- list(
    ribosomal = list(c("Rps3", "Keep", "Rpl10"), 2L),
    mito = list(c("mt-Nd1", "Keep", "mt-Co1"), 2L),
    hemoglobin = list(c("Hba1", "Keep", "Hbb", "Hbg2"), 3L),
    cell_cycle = list(
      c("Mki67", "Top2a", "Ccnb1", "Keep", "Ccnb2", "Cdk1", "Cdc20",
        "Birc5", "Tyms"),
      8L
    ),
    sex = list(c("Xist", "Keep", "Rps4y1", "Ddx3y", "Kdm5d"), 4L)
  )

  for (category in names(cases)) {
    case <- cases[[category]]
    expect_message(
      filtered <- filter_confounder_genes(
        case[[1L]], drop = category
      ),
      paste0(category, " removed ", case[[2L]], " gene\\(s\\)\\."),
      info = category
    )
    expect_identical(filtered, "Keep", info = category)
  }
})

test_that("filter_confounder_genes() lists valid categories on error", {
  expect_error(
    filter_confounder_genes("TP53", drop = "pseudogenes"),
    "Valid categories are: \"ribosomal\", \"mito\", \"hemoglobin\", ",
    fixed = TRUE
  )
  expect_error(
    filter_confounder_genes("TP53", drop = "pseudogenes"),
    "\"cell_cycle\", \"sex\"",
    fixed = TRUE
  )
})

test_that("filter_confounder_genes() preserves survivor order", {
  symbols <- c("KEEP3", "RPS3", "KEEP1", "MKI67", "KEEP2")
  filtered <- suppressMessages(filter_confounder_genes(
    symbols, drop = c("ribosomal", "cell_cycle")
  ))

  expect_identical(filtered, c("KEEP3", "KEEP1", "KEEP2"))
})

test_that("active categories all report counts after any removal", {
  messages <- character()
  filtered <- withCallingHandlers(
    filter_confounder_genes(
      c("RPS3", "KEEP"), drop = c("ribosomal", "mito")
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_identical(filtered, "KEEP")
  # message() appends the newline, so conditionMessage() carries it.
  expect_identical(
    messages,
    c(
      "filter_confounder_genes(): ribosomal removed 1 gene(s).\n",
      "filter_confounder_genes(): mito removed 0 gene(s).\n"
    )
  )
})

test_that("filter_confounder_genes() is silent when nothing matches", {
  symbols <- c("TP53", "EGFR", "SLC11A2")
  expect_silent(filtered <- filter_confounder_genes(symbols))
  expect_identical(filtered, symbols)
})

test_that("filter_confounder_genes() can remove every input", {
  symbols <- c("RPS3", "MT-ND1", "HBA1", "MKI67", "XIST")
  messages <- character()
  filtered <- withCallingHandlers(
    filter_confounder_genes(symbols),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_identical(filtered, character())
  expect_identical(
    sub(" removed.*$", "", sub("^filter_confounder_genes\\(\\): ", "", messages)),
    c("ribosomal", "mito", "hemoglobin", "cell_cycle", "sex")
  )
})

test_that("an empty drop selection returns every symbol without a message", {
  symbols <- c("RPS3", "MT-ND1", "MKI67")
  expect_silent(
    filtered <- filter_confounder_genes(symbols, drop = character())
  )
  expect_identical(filtered, symbols)
})

test_that("the hemoglobin rule keeps genes that merely start with HB", {
  # `^HB[ABDEG][0-9]*`, the rule as written in prose, removed all four of
  # these. HBEGF and HBS1L are not globins, and HBZ and HBM are.
  keep <- c("HBEGF", "HBS1L", "HBP1", "Hbegf")
  expect_silent(
    expect_identical(filter_confounder_genes(keep, drop = "hemoglobin"), keep)
  )
  globins <- c("HBA1", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBZ", "HBQ1",
               "HBM", "Hba-a1", "Hbb-bs")
  expect_message(
    out <- filter_confounder_genes(c(globins, "KEEP"), drop = "hemoglobin"),
    paste0("hemoglobin removed ", length(globins), " gene\\(s\\)\\.")
  )
  expect_identical(out, "KEEP")
})

# --- Defect found by the G5 cross-check ------------------------------------
# MSigDB ships retired symbols. HALLMARK_HYPOXIA lists ILVBL, the previous
# symbol for HACL2 (Entrez 10994); a SYMBOL-only lookup dropped it, producing a
# constant size offset of -3 against the reference implementation.

test_that("gene_to_entrez() resolves a retired symbol through the alias table", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  expect_message(
    mapped <- gene_to_entrez("ILVBL"),
    "alias table"
  )
  expect_identical(mapped, 10994L)
})

test_that("gene_to_entrez(alias_fallback = FALSE) keeps the strict behaviour", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  # The pre-1.1.0 behaviour stays reachable, so a caller who needs current
  # symbols only is not forced through the alias table.
  expect_warning(
    mapped <- gene_to_entrez("ILVBL", alias_fallback = FALSE),
    "1/1 symbols failed to map: ILVBL"
  )
  expect_identical(mapped, integer(0L))
})

test_that("a current symbol wins over the same string as another gene's alias", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  # TP53 is a current symbol; the alias pass must never override a live hit.
  expect_identical(gene_to_entrez("TP53"), 7157L)
})

test_that("the unmapped warning still fires for genuinely unknown symbols", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  expect_warning(
    mapped <- gene_to_entrez(c("TP53", "NOT_A_REAL_GENE")),
    "1/2 symbols failed to map: NOT_A_REAL_GENE"
  )
  expect_identical(mapped, 7157L)
})

test_that("gene_to_entrez() validates alias_fallback", {
  expect_error(gene_to_entrez("TP53", alias_fallback = NA),
               "`alias_fallback` must be TRUE or FALSE")
  expect_error(gene_to_entrez("TP53", alias_fallback = "yes"),
               "`alias_fallback` must be TRUE or FALSE")
})
