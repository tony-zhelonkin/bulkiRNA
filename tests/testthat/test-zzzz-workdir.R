test_that("the test suite leaves its working directory unchanged", {
  expect_identical(
    .bulkirna_test_workdir_entries(),
    .bulkirna_test_workdir_start
  )
})
