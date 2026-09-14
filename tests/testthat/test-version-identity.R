version_identity_repo_root <- function() {
  candidates <- c(
    "../..",
    ".",
    testthat::test_path("..", "..")
  )
  for (path in unique(candidates)) {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, ".git"))) {
      return(path)
    }
  }
  NULL
}

version_identity_git <- function(root, args) {
  git <- unname(Sys.which("git"))
  if (!nzchar(git)) return(NULL)

  output <- suppressWarnings(system2(
    git,
    c("-C", root, args),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = unname(output))
}

version_identity_context <- function() {
  root <- version_identity_repo_root()
  if (is.null(root)) {
    skip(paste(
      "Source-tree-only version identity audit: DESCRIPTION and .git are",
      "unavailable in a built tarball or installed package."
    ))
  }

  probe <- version_identity_git(root, c("rev-parse", "--show-toplevel"))
  if (is.null(probe)) {
    skip(paste(
      "Source-tree-only version identity audit: the git executable is",
      "unavailable, so repository identity cannot be checked."
    ))
  }
  if (probe$status != 0L) {
    stop(
      "Git metadata is present but unreadable: ",
      paste(probe$output, collapse = "\n"),
      call. = FALSE
    )
  }
  root
}

version_identity_git_output <- function(root, args) {
  result <- version_identity_git(root, args)
  if (is.null(result)) {
    stop("The git executable became unavailable during the audit.",
         call. = FALSE)
  }
  if (result$status != 0L) {
    stop(
      "Git command failed: ", paste(args, collapse = " "), "\n",
      paste(result$output, collapse = "\n"),
      call. = FALSE
    )
  }
  result$output
}

version_identity_description <- function(root) {
  unname(read.dcf(
    file.path(root, "DESCRIPTION"),
    fields = "Version"
  )[1L, 1L])
}

version_identity_release_tags <- function(root) {
  tags <- version_identity_git_output(root, c("tag", "--list"))
  grep("^v[0-9]+\\.[0-9]+\\.[0-9]+$", tags, value = TRUE)
}

version_identity_tag_commit <- function(root, tag) {
  version_identity_git_output(
    root,
    c("rev-parse", paste0("refs/tags/", tag, "^{commit}"))
  )[[1L]]
}

version_identity_tag_description <- function(root, tag) {
  result <- version_identity_git(
    root,
    c("show", paste0(tag, "^{commit}:DESCRIPTION"))
  )
  if (is.null(result)) {
    stop("The git executable became unavailable during the audit.",
         call. = FALSE)
  }
  if (result$status != 0L) return(NA_character_)

  version_line <- grep(
    "^Version:[[:space:]]*",
    result$output,
    value = TRUE
  )
  if (length(version_line) != 1L) return(NA_character_)
  sub("^Version:[[:space:]]*", "", version_line)
}

version_identity_expect_allowlist <- function(
    observed, allowlist, info = NULL) {
  expected <- names(allowlist)
  if (is.null(expected)) expected <- character(0L)
  expect_identical(sort(as.character(observed)), sort(expected), info = info)
}

# These are assertions only when the source checkout and its git metadata are
# visible. A skip means "could not check", not "checked and correct", so every
# unavailable-source path above reports an explicit reason.

test_that("a release version identifies its release commit and a clean tree", {
  root <- version_identity_context()
  version <- version_identity_description(root)
  is_release <- grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)
  release_tags <- version_identity_release_tags(root)
  expected_tag <- paste0("v", version)

  tag_exists <- !is_release || expected_tag %in% release_tags
  expect_true(
    tag_exists,
    info = paste0("Release version ", version, " needs tag ", expected_tag, ".")
  )

  if (is_release && tag_exists) {
    head <- version_identity_git_output(root, c("rev-parse", "HEAD"))[[1L]]
    tag_commit <- version_identity_tag_commit(root, expected_tag)
    expect_identical(
      head,
      tag_commit,
      info = paste0("Release version ", version, " must be checked out at ",
                    expected_tag, ".")
    )

    status <- version_identity_git_output(
      root,
      c("status", "--porcelain", "--untracked-files=all")
    )
    expect_identical(
      status,
      character(0L),
      info = "A release version requires a clean working tree."
    )
  }
})

test_that("post-release commits identify the next development version", {
  root <- version_identity_context()
  version <- version_identity_description(root)
  release_tags <- version_identity_release_tags(root)
  release_commits <- vapply(
    release_tags,
    function(tag) version_identity_tag_commit(root, tag),
    character(1L)
  )
  head <- version_identity_git_output(root, c("rev-parse", "HEAD"))[[1L]]
  head_is_release <- head %in% release_commits

  # This is deliberately stricter than the common x.y.z.9000-after-releasing-
  # x.y.z convention. Naming the next release makes that intent explicit and
  # lets this audit detect a development tree carrying a stale version.
  is_development <- grepl(
    "^[0-9]+\\.[0-9]+\\.[0-9]+\\.9000$",
    version
  )
  expect_true(
    head_is_release || is_development,
    info = paste0(
      "HEAD is not a release commit, so Version must be A.B.C.9000; got ",
      version, "."
    )
  )

  if (!head_is_release && is_development) {
    # expect_gt() takes no `info`, so the reason rides on expect_true().
    expect_true(
      length(release_tags) > 0L,
      info = "The development-version comparison requires release tags."
    )
    development_base <- numeric_version(sub("\\.9000$", "", version))
    release_versions <- numeric_version(sub("^v", "", release_tags))
    stale_tags <- release_tags[!(development_base > release_versions)]
    expect_identical(
      stale_tags,
      character(0L),
      info = paste(
        "The development base must be strictly greater than every release",
        "tag; it is not greater than:", paste(stale_tags, collapse = ", ")
      )
    )
  }
})

test_that("release tags contain their exact DESCRIPTION version", {
  root <- version_identity_context()
  release_tags <- version_identity_release_tags(root)
  tag_versions <- stats::setNames(vapply(
    release_tags,
    function(tag) version_identity_tag_description(root, tag),
    character(1L)
  ), release_tags)
  expected_versions <- sub("^v", "", names(tag_versions))
  violations <- names(tag_versions)[
    is.na(tag_versions) | unname(tag_versions) != expected_versions
  ]

  historical_exceptions <- c(
    v0.2.0 = paste(
      "This tag predates the repository's conversion into an installed R",
      "package; its commit has no DESCRIPTION file."
    )
  )
  expect_true(
    all(nzchar(historical_exceptions)) &&
      !any(grepl("[\r\n]", historical_exceptions)),
    info = "Every historical tag exception needs its own one-line reason."
  )

  violation_details <- vapply(violations, function(tag) {
    observed <- tag_versions[[tag]]
    if (is.na(observed)) observed <- "no readable Version"
    paste0("`", tag, "` has `", observed, "`.")
  }, character(1L))
  version_identity_expect_allowlist(
    violations,
    historical_exceptions,
    info = paste(
      paste(violation_details, collapse = " "),
      "Each release tag must contain its exact version or have a reasoned",
      "per-tag historical exception."
    )
  )
  version_identity_expect_allowlist(
    character(0L),
    NULL,
    info = "The historical exception check must be shape-stable when empty."
  )
})
