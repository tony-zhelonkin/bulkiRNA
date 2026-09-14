# Cut tests/fixtures/coresh-chunk-micro.rds from the real CoReSh compendium.
#
# Four real dataset objects, 300 genes each, chosen so that between them they
# carry every structural property that has broken this package: repeated Entrez
# ids, missing Entrez ids, and a matrix whose columns are principal components
# rather than samples. Hand-built fixtures contain only what their author
# thought to include; the NA ids present in 0.7% of real datasets passed the
# entire suite and were found by a compendium sweep instead.
#
# Requires the refcache mounted read-only and `qs2`, which is not in the image.
#
# Run:  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache -e R_LIBS=/rlib \
#         -v <msigdb-cache>:/cache -v <scratch-lib>:/rlib:ro \
#         -v /data2/users/shared/refcache/coresh/current:/coresh:ro \
#         -v <repo>:/pkg -w /pkg scdock-r-dev:v0.5.13 \
#         Rscript tests/fixtures/make_coresh_micro.R
#
# `qs2` into the scratch library first:
#   install.packages("qs2", lib = "/rlib", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages(library(qs2))

CHUNK_DIR <- "/coresh/preprocessed_chunks/mmu"   # snapshot syn66227307_20260721
OUT       <- "tests/fixtures/coresh-chunk-micro.rds"
N_GENES   <- 300L
N_SPECIAL <- 60L   # rows kept because they carry the property, not for size

# One object per property. `plain` is the control.
classify <- function(o) {
  if (anyNA(o$rownames)) {
    "na_ids"
  } else if (anyDuplicated(o$rownames) > 0) {
    "duplicate_ids"
  } else if (ncol(o$E1024) < length(o$samples)) {
    "pca_reduced"
  } else {
    "plain"
  }
}

# Keep the rows that carry the property first, then fill with ordinary rows.
# Taking a plain prefix of the matrix loses the duplicates and the NAs, which is
# how the first attempt at this fixture silently lost three of its four cases.
pick_rows <- function(o, key, n = N_GENES) {
  ids <- o$rownames
  special <- switch(
    key,
    na_ids        = which(is.na(ids)),
    duplicate_ids = which(ids %in% ids[duplicated(ids) & !is.na(ids)]),
    integer(0)
  )
  special <- utils::head(special, N_SPECIAL)
  rest <- setdiff(which(!is.na(ids)), special)
  sort(unique(c(special, utils::head(rest, n - length(special)))))
}

paths <- sort(list.files(CHUNK_DIR, pattern = "full_objects\\.qs2$",
                         full.names = TRUE))
if (!length(paths)) {
  stop("No chunk files under ", CHUNK_DIR, ". Mount the refcache read-only.",
       call. = FALSE)
}

want <- list()
for (path in paths) {
  chunk <- qs_read(path)
  for (o in chunk) {
    key <- classify(o)
    if (!is.null(want[[key]])) next
    keep <- pick_rows(o, key)
    # Columns are deliberately NOT subset. They are at most 20 wide, and cutting
    # them would destroy the `ncol(E1024) < nsamples` property that marks a
    # PCA-reduced matrix. `samples` and `nsamples` stay as the source wrote them.
    E <- o$E1024[keep, , drop = FALSE]
    want[[key]] <- list(
      gseId    = o$gseId,
      gplId    = o$gplId,
      E1024    = E,                       # integer, so the /1024 path is real
      rownames = o$rownames[keep],
      samples  = o$samples,
      nsamples = o$nsamples,
      # Recomputed for the subset: the stored value describes the full matrix.
      totalVar = sum((E / 1024)^2)
    )
  }
  if (length(want) == 4L) break
}

missing <- setdiff(c("plain", "duplicate_ids", "na_ids", "pca_reduced"),
                   names(want))
if (length(missing)) {
  stop("No dataset found carrying: ", paste(missing, collapse = ", "),
       ". Widen the search or pick the property differently.", call. = FALSE)
}

want <- want[c("plain", "duplicate_ids", "na_ids", "pca_reduced")]
for (key in names(want)) {
  o <- want[[key]]
  message(sprintf(
    "%-14s %-11s %3dx%-3d nsamples=%-4d reduced=%-5s NA=%-4d dup=%-5s totalVar=%.3f",
    key, o$gseId, nrow(o$E1024), ncol(o$E1024), o$nsamples,
    ncol(o$E1024) < o$nsamples, sum(is.na(o$rownames)),
    anyDuplicated(o$rownames) > 0, o$totalVar))
}

# version = 2 so the fixture stays readable on older R than the image ships.
saveRDS(want, OUT, version = 2)
message("wrote ", OUT, " (", file.size(OUT), " bytes)")
