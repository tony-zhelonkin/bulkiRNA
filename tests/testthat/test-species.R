test_that("the shared species resolver accepts every historical alias", {
  aliases <- c(
    human = "human",
    scientific_human = "Homo sapiens",
    underscore_human = "Homo_sapiens",
    code_human = "hsa",
    upper_human = "HOMO SAPIENS",
    mouse = "mouse",
    scientific_mouse = "Mus musculus",
    underscore_mouse = "Mus_musculus",
    code_mouse = "mmu",
    upper_mouse = "MOUSE"
  )
  observed <- lapply(aliases, bulkiRNA:::.species)

  expect_identical(
    vapply(observed, `[[`, character(1L), "scientific"),
    c(
      human = "Homo sapiens",
      scientific_human = "Homo sapiens",
      underscore_human = "Homo sapiens",
      code_human = "Homo sapiens",
      upper_human = "Homo sapiens",
      mouse = "Mus musculus",
      scientific_mouse = "Mus musculus",
      underscore_mouse = "Mus musculus",
      code_mouse = "Mus musculus",
      upper_mouse = "Mus musculus"
    )
  )
  expect_identical(
    vapply(observed, `[[`, character(1L), "code"),
    c(
      human = "hsa", scientific_human = "hsa", underscore_human = "hsa",
      code_human = "hsa", upper_human = "hsa", mouse = "mmu",
      scientific_mouse = "mmu", underscore_mouse = "mmu",
      code_mouse = "mmu", upper_mouse = "mmu"
    )
  )
})

test_that("the shared species records pin every representation", {
  expect_identical(
    bulkiRNA:::.species("human"),
    list(
      common = "human",
      scientific = "Homo sapiens",
      code = "hsa",
      orgdb = "org.Hs.eg.db",
      biomart_dataset = "hsapiens_gene_ensembl",
      gatom_short = "Hs",
      gatom_download = "Homo_sapiens"
    )
  )
  expect_identical(
    bulkiRNA:::.species("mouse"),
    list(
      common = "mouse",
      scientific = "Mus musculus",
      code = "mmu",
      orgdb = "org.Mm.eg.db",
      biomart_dataset = "mmusculus_gene_ensembl",
      gatom_short = "Mm",
      gatom_download = "Mus_musculus"
    )
  )
})

test_that("scientific-name partial matching remains available", {
  expect_identical(
    bulkiRNA:::.species("Homo")$scientific,
    "Homo sapiens"
  )
  expect_identical(
    bulkiRNA:::.species("Mus muscul")$scientific,
    "Mus musculus"
  )
})

test_that("species errors name every accepted spelling", {
  err <- tryCatch(
    bulkiRNA:::.species("Rattus norvegicus"),
    error = conditionMessage
  )
  for (alias in c("human", "Homo sapiens", "hsa", "mouse",
                  "Mus musculus", "mmu")) {
    expect_match(err, alias, fixed = TRUE)
  }
})

test_that("custom species remain reachable only through the gs_db formatter", {
  expect_identical(
    bulkiRNA:::.gsdb_species_label("Danio_rerio"),
    "Danio rerio"
  )
  expect_identical(bulkiRNA:::.gsdb_species_label("Mus"), "Mus")
  expect_identical(bulkiRNA:::.species("Homo")$scientific, "Homo sapiens")
  expect_error(
    bulkiRNA:::.gsdb_species_label(9606L),
    "`species` must be a single non-empty string",
    fixed = TRUE
  )
  expect_error(
    bulkiRNA:::.species("Danio_rerio"),
    "`species` must be one of"
  )
})

test_that("strict species errors use stable quotes", {
  expect_error(
    bulkiRNA:::.species("rat"),
    'got "rat".',
    fixed = TRUE
  )
})
