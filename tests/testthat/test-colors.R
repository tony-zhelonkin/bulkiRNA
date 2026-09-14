# Contract tests for the vendored palettes and gs_palette().

test_that("bulki_palettes() returns hex vectors and names its defaults", {
  pals <- bulki_palettes()
  expect_type(pals, "list")
  expect_true(all(c("hat", "heatmap0", "heatmap2", "heatmap3") %in%
                    names(pals)))
  # Every entry must be usable as a colour vector.
  for (nm in names(pals)) {
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pals[[nm]])), info = nm)
  }
  expect_length(bulki_palettes("hat"), 10L)
  expect_error(bulki_palettes("nope"), "`name` is not a known palette")
})

test_that("the vendored hexes match the upstream ltc values", {
  # Transcription guard: these are copied constants, so a typo would be
  # invisible everywhere except in a figure nobody would think to question.
  expect_identical(
    bulki_palettes("hat"),
    c("#efb306", "#eb990c", "#e8351e", "#cd023d", "#852f88",
      "#4e54ac", "#0f8096", "#7db954", "#17a769", "#000000")
  )
  expect_identical(
    bulki_palettes("heatmap2"),
    c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")
  )
})

test_that("the diverging default keeps a near-white midpoint and down-is-blue", {
  # scale_fill_gradient2() needs a neutral midpoint to read as zero, and the
  # published figures established blue as the down direction.
  div <- .gs_diverging_colours()
  expect_identical(names(div), c("low", "mid", "high"))
  expect_identical(toupper(unname(div[["mid"]])), "#F7F7F7")
  rgb_low <- grDevices::col2rgb(div[["low"]])
  rgb_high <- grDevices::col2rgb(div[["high"]])
  expect_gt(rgb_low["blue", 1], rgb_low["red", 1])
  expect_gt(rgb_high["red", 1], rgb_high["blue", 1])
})

test_that("the diverging default beats the red-ended ramps under CVD", {
  # The reason blue-orange is the default, encoded so it cannot be swapped for
  # a red-ended ramp on aesthetic grounds alone. A reader with protanopia has
  # to be able to tell the two ends apart; red loses the most lightness there.
  skip_if_not_installed("colorspace")

  lab <- function(hex) {
    grDevices::convertColor(t(grDevices::col2rgb(hex) / 255), "sRGB", "Lab")
  }
  ends_apart <- function(pal, sim) {
    l <- lab(sim(pal))
    sqrt(sum((l[1, ] - l[3, ])^2))
  }
  worst <- function(pal) {
    min(vapply(
      list(colorspace::deutan, colorspace::protan, colorspace::tritan),
      function(sim) ends_apart(pal, sim), numeric(1L)
    ))
  }

  default <- unname(.gs_diverging_colours())
  rdbu <- unname(bulki_palettes("heatmap2")[c(5L, 3L, 1L)])
  rdylbu <- unname(bulki_palettes("heatmap3")[c(5L, 3L, 1L)])

  expect_gt(worst(default), worst(rdbu))
  expect_gt(worst(default), worst(rdylbu))
})

test_that("gs_palette() is named, order-invariant and deduplicating", {
  a <- gs_palette(c("SET_B", "SET_A", "SET_C"))
  b <- gs_palette(c("SET_C", "SET_B", "SET_A"))
  expect_named(a, c("SET_B", "SET_A", "SET_C"))
  # Returned in the caller's order, but assigned by sorted id.
  expect_identical(a[["SET_A"]], b[["SET_A"]])
  expect_identical(a[["SET_C"]], b[["SET_C"]])
  expect_length(gs_palette(c("A", "A", "B")), 2L)
  expect_length(gs_palette(character(0L)), 0L)
})

test_that("gs_palette() interpolates past the palette length", {
  ids <- paste0("S", seq_len(25L))
  pal <- gs_palette(ids)
  expect_length(pal, 25L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}", pal)))
  # Distinct colours, since a duplicate would silently merge two curves.
  expect_gt(length(unique(pal)), 20L)
})

test_that("gs_palette() spreads picks across a qualitative palette", {
  # Not the first n: hat is ordered as a hue wheel, so a prefix is a run of
  # neighbouring hues. Asking for 5 must not return amber next to orange.
  five <- unname(gs_palette(paste0("S", 1:5)))
  expect_length(five, 5L)
  expect_false(identical(five, bulki_palettes("hat")[1:5]))
  expect_identical(five[[1L]], bulki_palettes("hat")[[1L]])
  expect_identical(five[[5L]], bulki_palettes("hat")[[10L]])

  # Pairwise separation in CIE Lab must beat what the prefix would give.
  skip_if_not_installed("colorspace")
  lab <- function(hex) {
    grDevices::convertColor(t(grDevices::col2rgb(hex) / 255), "sRGB", "Lab")
  }
  min_gap <- function(hex) min(stats::dist(lab(hex)))
  expect_gt(min_gap(five), min_gap(bulki_palettes("hat")[1:5]))
})

test_that("gs_palette() accepts a name or an explicit vector", {
  expect_identical(
    unname(gs_palette(c("A", "B"), palette = "minou")),
    bulki_palettes("minou")[c(1L, 6L)]
  )
  expect_identical(
    unname(gs_palette(c("A", "B"), palette = c("#123456", "#654321"))),
    c("#123456", "#654321")
  )
  expect_error(gs_palette(c("A", "B"), palette = 1),
               "`palette` must be a palette name")
  expect_error(gs_palette(1:3), "`ids` must be a character")
  expect_error(gs_palette(c("A", "")), "missing or empty identifiers")
})

test_that("a frozen gs_palette() gives one colour per id across figures", {
  # The documented way to get cross-figure stability: one lookup over the full
  # universe, reused. Two figures drawing different subsets must agree.
  universe <- paste0("SET_", LETTERS[1:8])
  frozen <- gs_palette(universe)
  ranks <- stats::setNames(seq(3, -3, length.out = 60L),
                           paste0("G", seq_len(60L)))
  sets <- stats::setNames(
    lapply(seq_along(universe), function(i) paste0("G", seq(i, 60L, by = 7L))),
    universe
  )

  guide <- function(ids) {
    p <- gs_plot_running(sets, ranks = ranks, pathways = ids,
                         palette = frozen)
    gd <- ggplot2::get_guide_data(p[[1]], "colour")
    stats::setNames(gd$colour, as.character(gd$.value))
  }
  one <- guide(c("SET_A", "SET_B", "SET_C"))
  two <- guide(c("SET_C", "SET_F", "SET_A"))
  expect_identical(one[["SET_A"]], two[["SET_A"]])
  expect_identical(one[["SET_C"]], two[["SET_C"]])
})
