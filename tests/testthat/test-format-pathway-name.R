test_that("format_pathway_name strips prefixes and applies exceptions", {
  # These are the values the pre-package implementation produces, verified
  # against scripts/GSEA/GSEA_plotting/format_pathway_names.R. The old
  # roxygen examples claimed hyphens ("TNFR1-Induced") that the code has
  # never produced; the golden baseline follows the code, so we do too.
  expect_equal(
    format_pathway_name("HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING"),
    "TNFR1 Induced NF-kappaB Signaling"
  )
  expect_equal(
    format_pathway_name("GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY"),
    "Type II Interferon Signaling Pathway"
  )
  expect_equal(
    format_pathway_name("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING"),
    "Class I MHC Mediated Antigen Processing"
  )
  expect_equal(
    format_pathway_name("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    "TNF-alpha Signaling via NF-kappaB"
  )
  expect_equal(
    format_pathway_name("HALLMARK_IL2_STAT5_SIGNALING"),
    "IL-2/STAT5 Signaling"
  )
})

test_that("dots are hierarchy separators, like underscores", {
  expect_equal(
    format_pathway_name("Metabolism.Lipid_metabolism.Fatty_acid_oxidation"),
    "Metabolism Lipid Metabolism Fatty Acid Oxidation"
  )
})

test_that("strip_prefix = FALSE keeps the collection prefix", {
  expect_equal(
    format_pathway_name("HALLMARK_APOPTOSIS", strip_prefix = FALSE),
    "Hallmark Apoptosis"
  )
})

test_that("use_formatting = FALSE falls back to title case", {
  expect_equal(
    format_pathway_name("HALLMARK_DNA_REPAIR", use_formatting = FALSE),
    "Dna Repair"
  )
})

test_that("TransportDB acronyms survive prefix stripping", {
  expect_equal(format_pathway_name("TRANSPORTDB_ABC"), "ABC")
  expect_equal(format_pathway_name("MITOPATHWAYS_OXPHOS"), "OXPHOS")
})

test_that("sentence-initial function words are capitalised", {
  input <- paste0(
    c("VIA", "AND", "OR", "OF", "IN", "TO", "BY", "FROM", "THE"),
    "_SIGNALING"
  )
  expected <- paste(
    c("Via", "And", "Or", "Of", "In", "To", "By", "From", "The"),
    "Signaling"
  )

  expect_equal(format_pathway_name(input), expected)
})

test_that("position does not change deliberate lowercase forms", {
  expect_equal(format_pathway_name("BETA_OXIDATION"), "beta Oxidation")
  expect_equal(
    format_pathway_name("CIS-REGULATORY_ELEMENT"),
    "cis-Regulatory Element"
  )
  expect_equal(format_pathway_name("MTOR_SIGNALING"), "mTOR Signaling")
})

test_that("function words remain lowercase within a label", {
  expect_equal(
    format_pathway_name("SIGNALING_VIA_THE_CELL"),
    "Signaling via the Cell"
  )
})

test_that("sentence-initial function words are capitalised vectorially", {
  input <- c(
    "THE_ATP-BINDING_CASSETTE_(ABC)_SUPERFAMILY",
    "FROM_BETA-OXIDATION_TO_MTOR_SIGNALING",
    "CELL_SIGNALING_BY_CIS-REGULATORY_ELEMENTS"
  )

  expect_equal(
    format_pathway_name(input),
    c(
      "The Atp-Binding Cassette (Abc) Superfamily",
      "From beta-Oxidation to mTOR Signaling",
      "Cell Signaling by cis-Regulatory Elements"
    )
  )
})

test_that("the function is vectorised and length-preserving", {
  x <- c("HALLMARK_APOPTOSIS", "KEGG_RIBOSOME", "GOBP_AUTOPHAGY")
  out <- format_pathway_name(x)
  expect_length(out, 3L)
  expect_type(out, "character")
  expect_null(names(out))
})

test_that("empty and all-NA input come back unchanged", {
  expect_equal(format_pathway_name(character(0)), character(0))
  expect_equal(format_pathway_name(NA_character_), NA_character_)
})

test_that("NA elements survive alongside real names", {
  out <- format_pathway_name(c("HALLMARK_APOPTOSIS", NA))
  expect_equal(out[[1L]], "Apoptosis")
  expect_true(is.na(out[[2L]]))
})

test_that(".gs_wrap_label splits long labels and leaves short ones alone", {
  short <- "Apoptosis"
  expect_equal(bulkiRNA:::.gs_wrap_label(short, width = 20), short)

  # 23 characters: over `width`, under 1.5 * `width`, so two lines.
  two <- bulkiRNA:::.gs_wrap_label("One two three four five", width = 20)
  expect_equal(lengths(strsplit(two, "\n")), 2L)

  three <- bulkiRNA:::.gs_wrap_label(
    paste(rep("word", 20), collapse = " "), width = 20
  )
  expect_equal(lengths(strsplit(three, "\n")), 3L)
})

test_that(".gs_wrap_label never drops or duplicates words", {
  txt <- paste(paste0("w", seq_len(17)), collapse = " ")
  out <- bulkiRNA:::.gs_wrap_label(txt, width = 10)
  expect_equal(
    unlist(strsplit(out, "[ \n]")),
    unlist(strsplit(txt, " "))
  )
})
