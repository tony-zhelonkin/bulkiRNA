# Colour palettes vendored from the `ltc` package.
#
#   https://github.com/loukesio/ltc-color-palettes
#   Copyright (c) 2021 ltc authors. Released under the MIT licence; the full
#   notice ships as inst/LICENSE.note.
#
# Vendored rather than depended on: these are 32 vectors of hex codes, and a
# package dependency for data is a build-and-install cost with nothing to
# maintain on the other side. Kept as plain R source (not `data/`, which
# .Rbuildignore excludes, and not sysdata.rda) so the values stay greppable and
# a colour change shows up as a readable diff.

#' Colour palettes available to `bulkiRNA` renderers
#'
#' The palettes are vendored from the `ltc` package (MIT, see
#' `inst/LICENSE.note`).
#'
#' `hat` (10 colours) is the **qualitative default** -- per-pathway running-sum
#' curves, and any per-family or per-class grouping. Subsets are spread across
#' it rather than taken as a prefix, because it is ordered as a hue wheel; see
#' [gs_palette()].
#'
#' The **diverging default is not from this table**. It stays the house
#' blue-white-orange ramp, which measurably survives colour-vision deficiency
#' better than the red-ended alternatives here. `heatmap2` (ColorBrewer RdBu)
#' and `heatmap3` (RdYlBu) are available for callers who want them: both have a
#' near-white midpoint, so `bulki_palettes("heatmap2")[c(5, 3, 1)]` is a valid
#' length-3 `palette` for [gs_plot_bar()], [gs_plot_dot()] and
#' [gs_plot_heatmap()]. `heatmap0` (9 colours) suits a sequential scale.
#'
#' @param name Optional palette name. `NULL` (default) returns the whole named
#'   list, which is what you want when auditing or comparing palettes.
#' @return A named list of character vectors when `name` is `NULL`, otherwise
#'   one unnamed character vector of hex colours.
#' @examples
#' names(bulki_palettes())
#' bulki_palettes("hat")
#' @export
bulki_palettes <- function(name = NULL) {
  pals <- .BULKI_PALETTES
  if (is.null(name)) return(pals)
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("`name` must be a single palette name, or NULL.", call. = FALSE)
  }
  if (!name %in% names(pals)) {
    stop("`name` is not a known palette: ", encodeString(name, quote = "\""),
         ". Available: ", paste(names(pals), collapse = ", "), ".",
         call. = FALSE)
  }
  pals[[name]]
}

#' The vendored palette table
#'
#' @return A named list of character vectors of hex colours.
#' @keywords internal
.BULKI_PALETTES <- list(
  paloma = c("#83AF9B", "#C8C8A9", "#f8da8a", "#f7bf95", "#fe8ca1"),
  maya = c("#3d5a80", "#98c1d9", "#e0fbfc", "#ee6c4d", "#293241"),
  dora = c("#52777A", "#542437", "#C02942", "#D95B43", "#ECD078"),
  ploen = c("#3F5671", "#83A1C3", "#CEB5C8", "#FAC898", "#B17776"),
  olga = c("#c9e3c2", "#8bc8cb", "#eccd80", "#f5ab70", "#9c87a1"),
  mterese = c("#f7ddaa", "#fac3ad", "#f897a1", "#9298BA", "#9cbeed"),
  gaby = c("#fceaab", "#f1a890", "#a8c4cc", "#82A0C2", "#85496F"),
  franscoise = c("#5980B1", "#b96a8d", "#A55062", "#E05256", "#E9A986"),
  fernande = c("#ff7676", "#F9D662", "#7cab7d", "#75B7D1"),
  sylvie = c("#E8B961", "#E88170", "#C6BDE8", "#5DB7C4", "#FD95BC"),
  expevo = c("#FC4E07", "#E7B800", "#00AFBB", "#8B4769", "#1d457f",
             "#808080"),
  minou = c("#00798c", "#d1495b", "#edae49", "#66a182", "#2e4057",
            "#8d96a3"),
  kiss = c("#FF7C7E", "#FEC300", "#9E3F71", "#31BCBA", "#E20035"),
  hat = c("#efb306", "#eb990c", "#e8351e", "#cd023d", "#852f88",
          "#4e54ac", "#0f8096", "#7db954", "#17a769", "#000000"),
  reading = c("#EFBC68", "#919F89", "#EDBDAE", "#57717C", "#5F97A4",
              "#CAEAC8", "#95A1AE", "#C8CFD6"),
  alger = c("#000000", "#1A5B5B", "#ACC8BE", "#F4AB5C", "#D1422F"),
  trio1 = c("#0E7175", "#FD7901", "#C35BCA"),
  trio2 = c("#89973D", "#E8B92F", "#A45E41"),
  trio3 = c("#E69F00", "#56B4E9", "#009E73"),
  trio4 = c("#94475E", "#364C54", "#E5A11F"),
  heatmap0 = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
               "#EE9B00", "#CA6702", "#AE2012", "#9B2226"),
  pantone23 = c("#7A92A5", "#1F2C43", "#FFB000", "#842c48", "#46483d"),
  remains = c("#69326E", "#EEEDC0", "#FF6D1F", "#EED455"),
  midnight = c("#16232A", "#FF5B04", "#075056", "#E4EEF0"),
  lincoln = c("#EEE9DF", "#C9C1B1", "#2C3B4D", "#FFB162", "#A35139",
              "#1B2632"),
  luminaries = c("#FF5B04", "#075056", "#233038", "#FDF6E3", "#F4D47C",
                 "#D3DBDD"),
  seafarer = c("#013D5A", "#FCF3E3", "#BDD3CE", "#708C69", "#E4A25B"),
  shuggie = c("#5B5F8D", "#9BB29E", "#DA6B51", "#F1DCBA", "#484149"),
  heatmap1 = c("#4d7799", "#7fa4c4", "#c5c8d4", "#d48e95", "#b5515b"),
  heatmap2 = c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0"),
  heatmap3 = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
  casa_natal = c("#245E55", "#ED773C", "#808BC5", "#C63F3E", "#EAC119",
                 "#EAA7C7", "#9ED6DF", "#1D1D1B", "#EAE4DA")
)

#' The default qualitative palette name
#'
#' @return A character scalar.
#' @keywords internal
.BULKI_QUALITATIVE <- "hat"

#' Expand a palette to an arbitrary length
#'
#' Spreads the picks across the palette when it is long enough, and
#' interpolates when it is not.
#'
#' Spreading, rather than taking the first `n`: a qualitative palette is
#' usually ordered as a hue wheel, so a contiguous prefix is a set of
#' *neighbouring* hues. `hat`'s first five are amber, orange, red, crimson and
#' purple -- all warm, and the first two are nearly indistinguishable on a
#' plotted line. Taking every other entry instead gives five hues that are
#' actually telling apart. The whole point of a qualitative scale is
#' separation, so separation is what the subset optimises.
#'
#' @param pal Character vector of hex colours.
#' @param n Number of colours required.
#' @return A character vector of `n` colours.
#' @keywords internal
.bulki_expand_palette <- function(pal, n) {
  if (!length(pal)) {
    stop("`pal` must contain at least one colour.", call. = FALSE)
  }
  if (n <= 0) return(character(0L))
  if (n == 1L) return(pal[[1L]])
  if (n <= length(pal)) {
    idx <- unique(round(seq(1, length(pal), length.out = n)))
    # Rounding can collide on short palettes; fall back to a prefix rather than
    # return fewer colours than were asked for.
    if (length(idx) < n) idx <- seq_len(n)
    return(pal[idx])
  }
  grDevices::colorRampPalette(pal)(n)
}

#' A colour per identifier, stable for a given set of identifiers
#'
#' Returns colours **named by identifier**, which is what every `bulkiRNA`
#' renderer maps its colour aesthetic to. Assignment follows `sort(ids)`, so
#' the same identifiers always receive the same colours whatever order they
#' were requested in.
#'
#' # The stability contract, stated rather than implied
#'
#' This is stable for a *fixed set* of identifiers, not across different sets.
#' Ask for the top 5 pathways of one contrast and the top 5 of another and the
#' two figures will disagree, because the sets differ. A renderer cannot fix
#' this: it only ever sees the identifiers it was asked to draw.
#'
#' To get one colour per pathway across a whole project, call `gs_palette()`
#' **once** over the full universe of identifiers you will ever plot and pass
#' the result to every call:
#'
#' ```
#' pal <- gs_palette(all_pathway_ids)
#' gs_plot_running(res_a, ranks_a, db = db, palette = pal)
#' gs_plot_running(res_b, ranks_b, db = db, palette = pal)
#' ```
#'
#' A palette with more names than a figure plots is fine -- renderers subset it
#' by identifier.
#'
#' @param ids Character vector of identifiers: pathway ids, TE families, or any
#'   other categorical key. Duplicates are collapsed.
#' @param palette A palette name from [bulki_palettes()], or an explicit
#'   character vector of colours.
#' @return A character vector of hex colours, named by identifier, in the order
#'   `ids` was given.
#' @examples
#' gs_palette(c("SET_B", "SET_A"))
#'
#' # Order-invariant: the same ids get the same colours either way round.
#' identical(
#'   gs_palette(c("A", "B"))[["A"]],
#'   gs_palette(c("B", "A"))[["A"]]
#' )
#' @export
gs_palette <- function(ids, palette = .BULKI_QUALITATIVE) {
  if (!is.character(ids) && !is.factor(ids)) {
    stop("`ids` must be a character or factor vector.", call. = FALSE)
  }
  ids <- unique(as.character(ids))
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop("`ids` must not contain missing or empty identifiers.", call. = FALSE)
  }
  if (!length(ids)) return(stats::setNames(character(0L), character(0L)))

  pal <- if (is.character(palette) && length(palette) == 1L &&
               palette %in% names(.BULKI_PALETTES)) {
    bulki_palettes(palette)
  } else if (is.character(palette) && length(palette) > 0L) {
    palette
  } else {
    stop("`palette` must be a palette name from `bulki_palettes()` or a ",
         "character vector of colours.", call. = FALSE)
  }

  # Sorted assignment is what makes this order-invariant. Radix ordering keeps
  # it independent of LC_COLLATE, so a figure does not change colour with the
  # locale the container happens to run under.
  keys <- ids[order(ids, method = "radix")]
  cols <- .bulki_expand_palette(pal, length(keys))
  stats::setNames(cols, keys)[ids]
}
