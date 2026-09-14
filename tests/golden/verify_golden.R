# Thin wrapper: recompute every golden case and diff against tests/golden/data.
# The case definitions live in capture_golden.R so there is one source of truth.
#
# Run:  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
#         -v <cache>:/cache -v <repo>:/pkg -w /pkg \
#         scdock-r-dev:v0.5.11 Rscript tests/golden/verify_golden.R
#
# Exits non-zero on any FAIL/NEW. This is Step C's gate.
commandArgs <- function(trailingOnly = FALSE) "--verify"
source("tests/golden/capture_golden.R")
