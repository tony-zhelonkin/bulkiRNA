test_that("theme_bulki returns a theme that composes with a plot", {
  th <- theme_bulki()
  expect_s3_class(th, "theme")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point() + th
  expect_s3_class(p, "ggplot")
})

test_that("the 14 pt base-size floor is enforced", {
  expect_equal(theme_bulki(base_size = 8)$text$size, 14)
  expect_equal(theme_bulki()$text$size, 14)
  expect_equal(theme_bulki(base_size = 18)$text$size, 18)
})

test_that("grid = FALSE removes the grid and grid = TRUE restores it", {
  expect_s3_class(theme_bulki()$panel.grid.major, "element_blank")
  expect_s3_class(theme_bulki(grid = TRUE)$panel.grid.major, "element_line")
  expect_s3_class(theme_bulki(grid = TRUE)$panel.grid.minor, "element_blank")
})

test_that("panel and plot backgrounds use transparent, never NA, colour", {
  th <- theme_bulki()
  expect_identical(th$panel.background$colour, "transparent")
  expect_identical(th$plot.background$colour, "transparent")
})

test_that("bad arguments error with the offending name", {
  expect_error(theme_bulki(base_size = "big"), "`base_size`")
  expect_error(theme_bulki(grid = NA), "`grid`")
})
