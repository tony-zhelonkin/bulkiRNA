#' Extract gene loadings for one CoReSh hit
#'
#' Projects every gene in a CoReSh dataset onto the unit-normalized profile of
#' the query genes and returns the genes with the largest absolute loadings.
#' Query IDs and duplicated reference IDs are matched with [base::match()], so
#' each Entrez ID contributes at most its first reference row, matching the
#' CoReSh search calculation. Missing reference Entrez IDs remain eligible for
#' the loading calculation but are omitted later when loadings are mapped to
#' symbols by [coresh_sets()].
#'
#' @param chunk_path Path to one `*_full_objects.qs2` CoReSh chunk.
#' @param gse_id A single GEO series accession present in the chunk.
#' @param query A non-empty integer vector of Entrez IDs. Duplicate IDs are
#'   removed with a message before loading extraction.
#' @param top_n Positive whole number of loadings to retain.
#' @param gpl Optional GEO platform accession. When omitted and the GSE occurs
#'   on more than one platform in the chunk, the first platform in
#'   locale-independent radix order is used with a warning.
#' @return A tibble with `gse`, `gpl`, `entrez`, signed `loading`, and
#'   absolute-loading `rank`, ordered from largest to smallest absolute
#'   loading.
#' @examples
#' \dontrun{
#' coresh_loadings(
#'   "001_full_objects.qs2", "GSE123", c(1L, 2L, 3L), gpl = "GPL570"
#' )
#' }
#' @export
coresh_loadings <- function(chunk_path, gse_id, query, top_n = 50L,
                            gpl = NULL) {
  if (!is.character(chunk_path) || length(chunk_path) != 1L ||
      is.na(chunk_path) || !nzchar(chunk_path)) {
    stop("`chunk_path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(chunk_path)) {
    stop("`chunk_path` does not exist: ", chunk_path, ".", call. = FALSE)
  }
  if (!is.character(gse_id) || length(gse_id) != 1L || is.na(gse_id) ||
      !nzchar(gse_id)) {
    stop("`gse_id` must be one non-empty GEO series accession.",
         call. = FALSE)
  }
  if (!is.null(gpl) &&
      (!is.character(gpl) || length(gpl) != 1L || is.na(gpl) ||
       !nzchar(gpl))) {
    stop("`gpl` must be one non-empty GEO platform accession or `NULL`.",
         call. = FALSE)
  }
  if (!is.integer(query) || !length(query) || anyNA(query)) {
    stop("`query` must be a non-empty integer vector of Entrez IDs.",
         call. = FALSE)
  }
  deduplicated <- unique(query)
  n_duplicates <- length(query) - length(deduplicated)
  if (n_duplicates) {
    message("Dropped ", n_duplicates, " duplicate Entrez ID",
            if (n_duplicates == 1L) "" else "s", " from `query`.")
  }
  query <- deduplicated
  top_n <- .coresh_positive_integer(top_n, "top_n")

  chunk <- .coresh_read_chunk(chunk_path)
  if (!is.list(chunk) || !length(chunk)) {
    stop("CoReSh `chunk_path` did not contain a non-empty object list: ",
         chunk_path, ".", call. = FALSE)
  }
  hits <- which(vapply(chunk, function(obj) {
    is.list(obj) && length(obj$gseId) == 1L &&
      identical(as.character(obj$gseId), gse_id)
  }, logical(1L)))
  if (!length(hits)) {
    stop("`gse_id` ", gse_id, " was not found in `chunk_path` ",
         chunk_path, ". Choose a GSE present in that chunk.", call. = FALSE)
  }
  for (hit in hits) {
    .coresh_validate_object(chunk[[hit]], paste0("GSE ", gse_id))
  }
  available <- sort(unique(vapply(
    chunk[hits], function(obj) as.character(obj$gplId), character(1L)
  )), method = "radix")
  selected_gpl <- .coresh_select_gpl(gse_id, available, gpl)
  selected_hits <- hits[vapply(chunk[hits], function(obj) {
    identical(as.character(obj$gplId), selected_gpl)
  }, logical(1L))]
  if (length(selected_hits) > 1L) {
    stop("GSE ", gse_id, " on platform ", selected_gpl, " occurs ",
         length(selected_hits), " times in ", chunk_path,
         "; the `(gse, gpl)` key is not unique.", call. = FALSE)
  }
  obj <- chunk[[selected_hits[[1L]]]]

  query_idx <- match(query, obj$rownames)
  query_idx <- query_idx[!is.na(query_idx)]
  if (length(query_idx) < 3L) {
    stop("`query` for `gse_id` ", gse_id, " has only ", length(query_idx),
         "/", length(query), " genes present; supply at least 3 matching ",
         "Entrez IDs for loading extraction.", call. = FALSE)
  }

  E <- obj$E1024 / 1024
  profile <- colSums(E[query_idx, , drop = FALSE])
  norm <- sqrt(sum(profile^2))
  if (!is.finite(norm) || norm <= 0) {
    stop("`query` for `gse_id` ", gse_id,
         " has a zero or non-finite profile; supply genes with non-zero ",
         "loadings in that dataset.", call. = FALSE)
  }
  loadings <- as.numeric(E %*% (profile / norm))
  ordering <- order(-abs(loadings), method = "radix", na.last = TRUE)
  keep <- utils::head(ordering, min(top_n, length(ordering)))

  tibble::tibble(
    gse = gse_id,
    gpl = selected_gpl,
    entrez = obj$rownames[keep],
    loading = loadings[keep],
    rank = seq_along(keep)
  )
}

#' Select one platform for a CoReSh accession
#'
#' @param gse A GEO series accession.
#' @param available Character vector of available platform accessions.
#' @param gpl Optional requested platform accession.
#' @return The requested platform, or the radix-first available platform.
#' @keywords internal
.coresh_select_gpl <- function(gse, available, gpl = NULL) {
  available <- sort(unique(as.character(available)), method = "radix")
  if (!is.null(gpl)) {
    if (!gpl %in% available) {
      stop("`gpl` ", gpl, " is unavailable for `gse` ", gse,
           "; choose one of: ", paste(available, collapse = ", "), ".",
           call. = FALSE)
    }
    return(gpl)
  }
  selected <- available[[1L]]
  if (length(available) > 1L) {
    warning(
      "GSE ", gse, " is available on multiple platforms: ",
      paste(available, collapse = ", "), ". Using ", selected,
      " because it is first in locale-independent radix order; specify `gpl` ",
      "to select explicitly.",
      call. = FALSE
    )
  }
  selected
}

#' Validate inputs for CoReSh set construction
#'
#' @param top_hits A data frame of ranked CoReSh hits.
#' @param queries Named list of Entrez queries.
#' @return `NULL`, invisibly.
#' @keywords internal
.coresh_validate_set_inputs <- function(top_hits, queries) {
  required <- c("query_name", "gse", "rank")
  if (!is.data.frame(top_hits)) {
    stop("`top_hits` must be a data frame returned by `coresh_search()`.",
         call. = FALSE)
  }
  missing <- setdiff(required, names(top_hits))
  if (length(missing)) {
    stop("`top_hits` is missing column(s): ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(top_hits$query_name) || any(!nzchar(top_hits$query_name)) ||
      anyNA(top_hits$gse) || any(!nzchar(top_hits$gse))) {
    stop("`top_hits$query_name` and `top_hits$gse` must contain non-missing, ",
         "non-empty strings.", call. = FALSE)
  }
  if ("gpl" %in% names(top_hits) &&
      (!is.character(top_hits$gpl) || anyNA(top_hits$gpl) ||
       any(!nzchar(top_hits$gpl)))) {
    stop("`top_hits$gpl` must contain non-missing, non-empty strings when ",
         "present.", call. = FALSE)
  }
  if (!is.numeric(top_hits$rank) || anyNA(top_hits$rank) ||
      any(!is.finite(top_hits$rank)) || any(top_hits$rank < 1) ||
      any(top_hits$rank != floor(top_hits$rank))) {
    stop("`top_hits$rank` must contain positive whole-number ranks.",
         call. = FALSE)
  }
  if (!is.list(queries) || !length(queries) || is.null(names(queries)) ||
      anyNA(names(queries)) || any(!nzchar(names(queries))) ||
      anyDuplicated(names(queries))) {
    stop("`queries` must be a non-empty named list with unique, non-empty ",
         "names.", call. = FALSE)
  }
  unknown <- setdiff(unique(top_hits$query_name), names(queries))
  if (length(unknown)) {
    stop("`top_hits$query_name` is absent from `queries`: ",
         paste0("\"", unknown, "\"", collapse = ", "), ".",
         call. = FALSE)
  }
  for (query_name in unique(top_hits$query_name)) {
    query <- queries[[query_name]]
    if (!is.integer(query) || anyNA(query) || length(unique(query)) < 3L) {
      stop("`queries` entry ", encodeString(query_name, quote = "\""),
           " must contain at least 3 unique, non-missing integer Entrez IDs.",
           call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Find the indexed chunk for one CoReSh hit
#'
#' @param row One row from `top_hits`.
#' @param index A tibble returned by [coresh_chunks()].
#' @return A list containing the unique `chunk` path and selected `gpl`.
#' @keywords internal
.coresh_hit_chunk <- function(row, index) {
  gse <- row$gse[[1L]]
  candidates <- index[index$gse == gse, , drop = FALSE]
  if (!nrow(candidates)) {
    stop("`row$gse` ", gse, " was not found in `index`; use a hit returned ",
         "by `coresh_search()` for the same chunk snapshot.", call. = FALSE)
  }
  requested_gpl <- if ("gpl" %in% names(row)) row$gpl[[1L]] else NULL
  selected_gpl <- .coresh_select_gpl(gse, candidates$gpl, requested_gpl)
  selected <- candidates[candidates$gpl == selected_gpl, , drop = FALSE]
  paths <- sort(unique(as.character(selected$chunk)), method = "radix")
  if (length(paths) != 1L || nrow(selected) != 1L) {
    stop("GSE ", gse, " on platform ", selected_gpl, " occurs ",
         nrow(selected), " times in the CoReSh chunk index; the `(gse, gpl)` ",
         "key is not unique.", call. = FALSE)
  }
  list(chunk = paths[[1L]], gpl = selected_gpl)
}

#' Remove redundant CoReSh sets by Jaccard overlap
#'
#' Input is already in explicit priority order. A later set is removed when
#' its Jaccard overlap with any retained earlier set is strictly greater than
#' `threshold`, so the lower CoReSh rank wins independently of caller order.
#'
#' @param sets Named list of character gene vectors.
#' @param threshold Numeric Jaccard threshold in `[0, 1]`.
#' @return Logical vector selecting retained sets.
#' @keywords internal
.coresh_dedupe_sets <- function(sets, threshold) {
  keep <- rep(TRUE, length(sets))
  for (i in seq_along(sets)) {
    if (!keep[[i]]) next
    for (j in seq_len(i - 1L)) {
      if (!keep[[j]]) next
      union_size <- length(union(sets[[i]], sets[[j]]))
      overlap <- if (union_size) {
        length(intersect(sets[[i]], sets[[j]])) / union_size
      } else {
        0
      }
      if (overlap > threshold) {
        keep[[i]] <- FALSE
        break
      }
    }
  }
  keep
}

#' Summarise CoReSh queries for database provenance
#'
#' Database provenance identifies caller-owned queries by name and unique
#' Entrez-ID count without copying their full definitions into every database.
#'
#' @param queries A validated named list of Entrez queries.
#' @return A character scalar in the form `name(size), name(size)`.
#' @keywords internal
.coresh_query_summary <- function(queries) {
  sizes <- vapply(queries, function(query) length(unique(query)), integer(1L))
  paste0(names(queries), "(", sizes, ")", collapse = ", ")
}

#' Build gene sets from ranked CoReSh hits
#'
#' For every hit, [coresh_loadings()] extracts the strongest gene loadings and
#' [entrez_to_gene()] maps them to symbols. Sets outside the size bounds are
#' dropped. Remaining sets are ordered explicitly by CoReSh rank, query name,
#' GSE, and platform before Jaccard de-duplication; therefore a more highly
#' ranked hit wins, rather than whichever row the caller happened to supply
#' first.
#'
#' Chunk lookup failures, extraction or mapping failures, set-name collisions,
#' and size-filtered drops are counted and reported separately. If every
#' attempted extraction fails, the function stops. A run that extracts hits
#' successfully but produces no in-range sets returns an empty database with
#' an always-on message.
#'
#' @param top_hits A data frame returned by [coresh_search()], usually filtered
#'   to the desired number of hits per query. Required columns are
#'   `query_name`, `gse`, and `rank`; `gpl` selects the exact dataset when
#'   present. For a historical table without `gpl`, an ambiguous accession
#'   warns and uses the radix-first available platform.
#' @param queries A non-empty named list of integer Entrez vectors.
#' @param chunk_dir Optional explicit path to an `hsa` or `mmu` chunk directory.
#' @param species One of `"human"`, `"hsa"`, `"mouse"`, or `"mmu"`.
#' @param top_n Positive whole number of absolute loadings to retain per hit.
#' @param min_size Positive whole-number minimum set size.
#' @param max_size Positive whole-number maximum set size.
#' @param jaccard_threshold Numeric threshold in `[0, 1]`. A later, lower
#'   priority set is removed when overlap is strictly greater than this value.
#' @param verbose Logical. Report the reason for each skipped hit in addition
#'   to the always-reported category counts.
#' @return A [gs_db()] with database-level `provenance` and a set-keyed
#'   `set_provenance` tibble. Database provenance records query names and their
#'   unique Entrez-ID counts, the snapshot, and the set-building parameters.
#'   It does not store the exact Entrez IDs; reproducing the database requires
#'   the caller's own query definitions. Set provenance contains `set_name`,
#'   `query_name`, `gse`, `gpl`, `chunk_path`, `loading_cutoff`, and
#'   `rank_in_coresh` for each retained set.
#' @examples
#' \dontrun{
#' db <- coresh_sets(
#'   hits[hits$rank <= 5L, ],
#'   list(iron = c(7037L, 4891L, 55240L)),
#'   species = "human"
#' )
#' }
#' @export
coresh_sets <- function(top_hits, queries, chunk_dir = NULL,
                        species = "human", top_n = 50L,
                        min_size = 15L, max_size = 500L,
                        jaccard_threshold = 0.8, verbose = FALSE) {
  .coresh_validate_set_inputs(top_hits, queries)
  query_summary <- .coresh_query_summary(queries)
  top_n <- .coresh_positive_integer(top_n, "top_n")
  min_size <- .coresh_positive_integer(min_size, "min_size")
  max_size <- .coresh_positive_integer(max_size, "max_size")
  if (min_size > max_size) {
    stop("`min_size` must not exceed `max_size`.", call. = FALSE)
  }
  if (min_size > top_n) {
    stop("`min_size` (", min_size, ") must not exceed `top_n` (", top_n,
         "); symbol mapping can only retain or reduce the `top_n` extracted ",
         "genes.", call. = FALSE)
  }
  if (!is.numeric(jaccard_threshold) || length(jaccard_threshold) != 1L ||
      is.na(jaccard_threshold) || !is.finite(jaccard_threshold) ||
      jaccard_threshold < 0 || jaccard_threshold > 1) {
    stop("`jaccard_threshold` must be one finite number in [0, 1].",
         call. = FALSE)
  }
  .coresh_logical_scalar(verbose, "verbose")
  species_info <- .species(species)

  index <- coresh_chunks(chunk_dir = chunk_dir, species = species)
  index_provenance <- attr(index, "provenance")
  top_hits$.input_order <- seq_len(nrow(top_hits))
  # Recycled explicitly: `order()` rejects a length-1 key beside length-n keys.
  gpl_order <- if ("gpl" %in% names(top_hits)) {
    top_hits$gpl
  } else {
    rep("", nrow(top_hits))
  }
  priority <- order(
    top_hits$rank,
    top_hits$query_name,
    top_hits$gse,
    gpl_order,
    top_hits$.input_order,
    na.last = TRUE,
    method = "radix"
  )
  top_hits <- top_hits[priority, , drop = FALSE]
  top_hits$.name_with_gpl <- FALSE
  if ("gpl" %in% names(top_hits) && nrow(top_hits)) {
    pair_key <- paste0(
      nchar(top_hits$query_name), ":", top_hits$query_name, ":",
      nchar(top_hits$gse), ":", top_hits$gse
    )
    platform_rows <- !duplicated(top_hits[c("query_name", "gse", "gpl")])
    platform_key <- pair_key[platform_rows]
    multiple <- duplicated(platform_key) |
      duplicated(platform_key, fromLast = TRUE)
    top_hits$.name_with_gpl <- pair_key %in% platform_key[multiple]
  }

  sets <- list()
  provenance_rows <- list()
  chunk_lookup_failures <- character()
  extraction_failures <- character()
  name_collisions <- character()
  size_drops <- character()
  successful_extractions <- 0L
  for (i in seq_len(nrow(top_hits))) {
    row <- top_hits[i, , drop = FALSE]
    label <- paste0(row$query_name[[1L]], "/", row$gse[[1L]])
    location <- tryCatch(.coresh_hit_chunk(row, index), error = identity)
    if (inherits(location, "error")) {
      chunk_lookup_failures <- c(
        chunk_lookup_failures,
        paste0(label, ": ", conditionMessage(location))
      )
      next
    }

    built <- tryCatch({
      loadings <- coresh_loadings(
        location$chunk,
        row$gse[[1L]],
        queries[[row$query_name[[1L]]]],
        top_n = top_n,
        gpl = location$gpl
      )
      used_gpl <- unique(loadings$gpl)
      if (length(used_gpl) != 1L || is.na(used_gpl) || !nzchar(used_gpl)) {
        stop("`coresh_loadings()` did not identify one platform.",
             call. = FALSE)
      }
      mapped_ids <- loadings$entrez[!is.na(loadings$entrez)]
      symbols <- entrez_to_gene(mapped_ids, species = species)
      genes <- unique(unname(symbols))
      list(loadings = loadings, genes = genes, gpl = used_gpl)
    }, error = identity)
    if (inherits(built, "error")) {
      extraction_failures <- c(
        extraction_failures,
        paste0(label, ": ", conditionMessage(built))
      )
      next
    }
    successful_extractions <- successful_extractions + 1L
    if (length(built$genes) < min_size || length(built$genes) > max_size) {
      size_drops <- c(
        size_drops,
        paste0(label, ": produced ", length(built$genes),
               " genes, outside [", min_size, ", ", max_size, "]")
      )
      next
    }

    set_name <- paste0(
      "CORESH_", row$query_name[[1L]], "_", row$gse[[1L]],
      if (row$.name_with_gpl[[1L]]) paste0("_", built$gpl) else ""
    )
    if (set_name %in% names(sets)) {
      name_collisions <- c(
        name_collisions,
        paste0(label, ": set name duplicates an earlier platform hit")
      )
      next
    }
    sets[[set_name]] <- built$genes
    provenance_rows[[set_name]] <- tibble::tibble(
      set_name = set_name,
      query_name = as.character(row$query_name[[1L]]),
      gse = as.character(row$gse[[1L]]),
      gpl = as.character(built$gpl),
      chunk_path = unname(as.character(location$chunk)),
      loading_cutoff = min(abs(built$loadings$loading)),
      rank_in_coresh = as.integer(row$rank[[1L]]),
      # Carried so coresh_labels() can build a web-UI-style label without
      # rejoining to the search ranking. `query_size` is the query size the
      # search scored; `n_genes` is the derived set size. They are not the same
      # number and a figure that conflates them misreports both.
      pct_var = as.numeric(row$pct_var[[1L]] %||% NA_real_),
      p_value = as.numeric(row$p_value[[1L]] %||% NA_real_),
      query_size = as.integer(row$size[[1L]] %||% NA_integer_),
      n_genes = length(built$genes)
    )
  }

  if (nrow(top_hits) && successful_extractions == 0L) {
    failures <- c(chunk_lookup_failures, extraction_failures)
    stop("`coresh_sets()` failed for all ", nrow(top_hits), " attempted hits: ",
         paste(failures, collapse = "; "), ".", call. = FALSE)
  }
  if (length(chunk_lookup_failures)) {
    message("coresh_sets(): chunk lookup failed for ",
            length(chunk_lookup_failures), " of ", nrow(top_hits), " hits.")
    if (verbose) {
      message(paste0("  - ", chunk_lookup_failures, collapse = "\n"))
    }
  }
  if (length(extraction_failures)) {
    message("coresh_sets(): loading extraction or symbol mapping failed for ",
            length(extraction_failures), " of ", nrow(top_hits), " hits.")
    if (verbose) {
      message(paste0("  - ", extraction_failures, collapse = "\n"))
    }
  }
  if (length(name_collisions)) {
    message("coresh_sets(): set-name collision skipped ",
            length(name_collisions), " of ", nrow(top_hits), " hits.")
    if (verbose) {
      message(paste0("  - ", name_collisions, collapse = "\n"))
    }
  }
  if (length(size_drops)) {
    message("coresh_sets(): size filter dropped ", length(size_drops),
            " of ", nrow(top_hits), " hits outside [", min_size, ", ",
            max_size, "] genes.")
    if (verbose) {
      message(paste0("  - ", size_drops, collapse = "\n"))
    }
  }

  if (length(sets)) {
    keep <- .coresh_dedupe_sets(sets, jaccard_threshold)
    sets <- sets[keep]
    provenance_rows <- provenance_rows[keep]
  }
  if (!length(sets)) {
    message("coresh_sets(): ", nrow(top_hits),
            " hits attempted, 0 sets produced.")
  }
  set_provenance <- if (length(provenance_rows)) {
    dplyr::bind_rows(provenance_rows)
  } else {
    tibble::tibble(
      set_name = character(),
      query_name = character(),
      gse = character(),
      gpl = character(),
      chunk_path = character(),
      loading_cutoff = numeric(),
      rank_in_coresh = integer()
    )
  }

  provenance_path <- index_provenance$path
  if (is.null(provenance_path)) {
    provenance_path <- if (is.null(chunk_dir)) {
      NA_character_
    } else {
      unname(as.character(chunk_dir))
    }
  }
  database_provenance <- list(
    source = index_provenance$source %||% "coresh",
    snapshot = index_provenance$snapshot %||% NA_character_,
    chunk_dir = provenance_path,
    species = species_info$scientific,
    n_chunks = index_provenance$n_chunks %||%
      as.integer(length(unique(index$chunk))),
    queries = query_summary,
    top_hits = if (nrow(top_hits)) {
      as.integer(max(top_hits$rank))
    } else {
      0L
    },
    top_n = top_n,
    min_size = min_size,
    max_size = max_size,
    jaccard_threshold = jaccard_threshold
  )

  # Display names come from coresh_labels(), not from the set id. The id
  # encodes the query name, which for a DE-seeded query is the caller's own
  # contrast label -- shown against an unrelated GEO accession, that reads as
  # if the bar were about that contrast. Titles are not available offline, so
  # the label is accession/platform/size/percent-of-variation; a caller with a
  # title lookup can re-derive richer labels from `set_provenance`.
  gs_db(
    sets,
    database = "coresh",
    species = species_info$scientific,
    pathway_names = if (nrow(set_provenance)) {
      coresh_labels(set_provenance)
    } else {
      NULL
    },
    database_label = "CoReSh-derived gene sets",
    set_provenance = set_provenance,
    provenance = database_provenance
  )
}
