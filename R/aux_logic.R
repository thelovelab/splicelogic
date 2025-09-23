#' Auxiliary logic functions
#' @param gr A GRanges object with exon annotations
#' @param left_exon A GRanges object with left exon(s) to match
#' @param right_exon A GRanges object with right exon(s) to match
#' @param type The type of overlap to consider when identifying matches.
#' Currently only "boundary" is implemented.
#' @return A tibble with columns from gr and two additional columns:
#' 'match_left' and 'match_right' indicating the number of overlaps with
#' the left and right exons, respectively.
#' @importFrom dplyr mutate as_tibble
#' @importFrom plyranges count_overlaps
#' @importFrom magrittr %>%
#' @export  
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

