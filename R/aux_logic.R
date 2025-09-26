#' Auxiliary logic functions
#' @param gr A GRanges object with exon annotations
#' @param left_exon A GRanges object with left exon(s) to match
#' @param right_exon A GRanges object with right exon(s) to match
#' @param type The type of overlap to consider when identifying matches.
#' @return A tibble with columns from gr and two additional columns:
#' 'match_left' and 'match_right' indicating the number of overlaps with
#' the left and right exons, respectively.
#' @importFrom dplyr mutate as_tibble
#' @importFrom plyranges count_overlaps
#' @importFrom magrittr %>%
get_matcher <- function(gr, left_exon, right_exon, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  matcher <- switch(
    type,
    "over" = {
        gr %>%
        plyranges::mutate(
            match_left  = plyranges::count_overlaps(., left_exon) >= 1,
            match_right = plyranges::count_overlaps(., right_exon) >= 1
        ) |>
        as_tibble()
    },
    "in" = {
      gr %>%
        plyranges::mutate(
          match_left  = . == left_exon,
          match_right = . == right_exon
        ) |>
        as_tibble()
    },
    "boundary" = {
      gr %>%
        plyranges::mutate(
          match_left  = GenomicRanges::end(.)   %in% GenomicRanges::end(left_exon),
          match_right = GenomicRanges::start(.) %in% GenomicRanges::start(right_exon)
        ) |>
        as_tibble()
    }
  )
return(matcher)
}

# Filter candidates that do not overlap any pos_exons, and are internal
# Then get the left and right exons for each candidate
# Return a named list with three GRanges objects: candidates, left_exons, right_exons
#' @param neg_exons A GRanges object with candidate negative exons (neg_exons)
#' @param pos_exons A GRanges object with positive exons (pos_exons)
#' @param gr The original GRanges object with all exons (for looking up left/right)
#' @return A named list with three GRanges objects: candidates, left_exons, right_exons
#' @importFrom plyranges filter_by_non_overlaps_directed slice
#' @importFrom GenomicRanges GRanges  
candidates_by_non_overlap <- function(neg_exons, pos_exons, gr, type) {
  # filter candidates that do not overlap any pos_exons
  candidates <- neg_exons |>
    plyranges::filter_by_non_overlaps_directed(pos_exons)

  # early return: no candidates -> empty list with GRanges objects
  if (length(candidates) == 0L) {
    return(list(
      candidates = candidates,
      left_exons  = GenomicRanges::GRanges(), 
      right_exons = GenomicRanges::GRanges()
    ))
  }

  candidates <- candidates |>
    plyranges::filter(internal)

  # keys of exons to the left and right of candidates
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1L)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1L)

  # get the actual exons for the candidates (preserves order of keys)
  left_exons  <- gr |> plyranges::slice(match(left_keys, key))
  right_exons <- gr |> plyranges::slice(match(right_keys, key))

  # return all three objects in a named list
  list(
    candidates = candidates,
    left_exons  = left_exons,
    right_exons = right_exons
  )
}

#' For a given candidate exon, find matching left and right exons in pos_exons
#' @param pos_exons A GRanges object with positive exons (pos_exons)
#' @param left_exon A GRanges object with the left exon to match
#' @param right_exon A GRanges object with the right exon to match
#' @param type The type of overlap to consider when identifying matches.
#' @return A list with two tibbles: left_tbl and right_tbl
#' @importFrom dplyr inner_join filter distinct rename
#' @importFrom plyranges filter_by_overlaps
#' @importFrom magrittr %>%
match_left_right <- function(pos_exons, left_exon, right_exon, type) {

    # split by tx_id (as tibbles) and check each transcript separately
    temp_pos <- get_matcher(
      pos_exons,
      left_exon  = left_exon,
      right_exon = right_exon,
      type = type
    )
    # left/right exon ranks per tx
    left_tbl <- temp_pos |>
      dplyr::filter(match_left == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(l = exon_rank)

    right_tbl <- temp_pos |>
      dplyr::filter(match_right == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(r = exon_rank)

    list(
      left_tbl = left_tbl, 
      right_tbl = right_tbl)
}