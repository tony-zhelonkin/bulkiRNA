# Contract tests for gs_scale_fonts() and gs_plot_size().

size_plot <- function() {
  ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) +
    ggplot2::geom_point()
}

text_size <- function(p) {
  th <- if (inherits(p, "patchwork")) p[[1]]$theme else p$theme
  th$text$size
}

test_that("gs_scale_fonts() grows text with canvas area", {
  small <- gs_scale_fonts(size_plot(), width = 7, height = 5)
  large <- gs_scale_fonts(size_plot(), width = 14, height = 10)
  expect_gt(text_size(large), text_size(small))
  # The 7 x 5 reference canvas is scale factor 1 by construction.
  expect_equal(text_size(small), 10, tolerance = 1e-9)
  # Doubling both dimensions doubles the linear scale factor.
  expect_equal(text_size(large), 20, tolerance = 1e-9)
})

test_that("gs_scale_fonts() distributes across a patchwork's panels", {
  pw <- patchwork::wrap_plots(size_plot(), size_plot(), ncol = 1L)
  scaled <- gs_scale_fonts(pw, width = 14, height = 10)
  expect_s3_class(scaled, "patchwork")
  # `&`, not `+`: every panel must move, not just the composition.
  expect_equal(scaled[[1]]$theme$text$size, 20, tolerance = 1e-9)
  expect_equal(scaled[[2]]$theme$text$size, 20, tolerance = 1e-9)
})

test_that("gs_scale_fonts() supersedes theme_bulki()'s 14 pt floor", {
  # The floor is deliberate for an unscaled theme, but once a canvas is known
  # the absolute size has to win, or a small figure carries oversized type.
  p <- size_plot() + theme_bulki(base_size = 14)
  scaled <- gs_scale_fonts(p, width = 4, height = 3, base_font_size = 10)
  expect_lt(text_size(scaled), 14)
})

test_that("gs_scale_fonts() validates its canvas", {
  expect_error(gs_scale_fonts(size_plot(), width = 0, height = 5),
               "`width` must be one finite positive number")
  expect_error(gs_scale_fonts(size_plot(), width = 7, height = NA),
               "`height` must be one finite positive number")
  expect_error(gs_scale_fonts("not a plot", width = 7, height = 5),
               "must be a ggplot or patchwork")
})

test_that("gs_plot_size() suggests per-database sizes with a running-sum rule", {
  expect_identical(
    gs_plot_size("bar", "Hallmark"),
    list(width = 7, height = 5, base_font_size = 10)
  )
  # A running-sum stacks three panels, so it gets 20% more height for the same
  # database. This is the one structural rule in the table.
  bar <- gs_plot_size("bar", "GO_BP")
  run <- gs_plot_size("running", "GO_BP")
  expect_identical(bar$width, run$width)
  expect_equal(run$height / bar$height, 1.2, tolerance = 1e-9)
  expect_equal(gs_plot_size("facet", "GO_BP")$height / bar$height, 1.4,
               tolerance = 1e-9)
})

test_that("gs_plot_size() normalises the database spelling", {
  expect_identical(gs_plot_size("dot", "GO_BP"), gs_plot_size("dot", "gobp"))
  expect_identical(gs_plot_size("dot", "GO BP"), gs_plot_size("dot", "gobp"))
  # An unrecognised or absent database falls back rather than erroring: the
  # table is a suggestion, and a custom database is not a mistake.
  expect_identical(gs_plot_size("dot", "CoReSh"), gs_plot_size("dot"))
  expect_identical(gs_plot_size("dot", NULL), gs_plot_size("dot"))
  expect_error(gs_plot_size("nope"), "should be one of")
})

test_that("gs_save() writes a patchwork and can skip font rescaling", {
  dir <- withr::local_tempdir()
  pw <- patchwork::wrap_plots(size_plot(), size_plot(), ncol = 1L)
  written <- gs_save(pw, file.path(dir, "pw"), width = 6, height = 5,
                     formats = "png", table = FALSE)
  expect_true(file.exists(written[[1]]))

  # rescale_fonts must not mutate the caller's object either way.
  before <- pw[[1]]$theme$text
  gs_save(pw, file.path(dir, "pw2"), width = 12, height = 9,
          formats = "png", table = FALSE)
  expect_identical(pw[[1]]$theme$text, before)

  plain <- gs_save(size_plot(), file.path(dir, "plain"), width = 6, height = 5,
                   formats = "png", table = FALSE, rescale_fonts = FALSE)
  expect_true(file.exists(plain[[1]]))
})
