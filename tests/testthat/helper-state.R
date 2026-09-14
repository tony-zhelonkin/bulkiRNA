# Record paths rather than a file allowlist so this covers future tests.
#
# Two entries are excluded because testthat and the graphics device own them,
# not the tests: `_snaps/` is testthat's snapshot directory, and `Rplots.pdf`
# is what base graphics writes when a plot is drawn with no device open. Both
# appear during a normal run and neither is a test writing outside its own
# temporary directory, which is what this guard exists to catch. Any other new
# entry fails the suite.
.bulkirna_test_workdir_ignore <- c("_snaps", "Rplots.pdf")

.bulkirna_test_workdir_entries <- function() {
  entries <- list.files(
    ".",
    all.files = TRUE,
    no.. = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )
  owned <- entries %in% .bulkirna_test_workdir_ignore |
    startsWith(entries, "_snaps/")
  sort(entries[!owned])
}

.bulkirna_test_workdir_start <- .bulkirna_test_workdir_entries()

# Attribute a leak to the test that caused it. The final suite test below is
# the enforcing assertion.
testthat::set_state_inspector(.bulkirna_test_workdir_entries)
