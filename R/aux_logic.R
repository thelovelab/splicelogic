#' Takes a GRanges object and left/right exons to compute matches.
#' i.e takes the pos_exons (GRanges) and left and right exons and returns
#' a tibble that is the same as pos_exons but with two additional columns
#' 'match_left' and 'match_right' indicating whether each exon matches the left
#'  and right exons respectively. The matching is done based on the 'type'
#' parameter which can be "in", "over", or "boundary".
#' @param gr A GRanges object with exon annotations
#' @param left_exon A GRanges object with left exon(s) to match
#' @param right_exon A GRanges object with right exon(s) to match
#' @param type The type of overlap to consider when identifying matches.
#' @return A tibble created from gr with two additional columns:
#' 'match_left' and 'match_right' indicating the number of overlaps with
#' the left and right exons, respectively.
#' @importFrom dplyr mutate as_tibble
#' @importFrom plyranges count_overlaps
#' @importFrom magrittr %>%
compute_matches <- function(gr, left_exon, right_exon, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  gr_matched <- switch(
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
return(gr_matched)
}

#' Filter candidates that do not overlap any pos_exons, and are internal
#' Then get the left and right exons for each candidate
#' Return a named list with three GRanges objects: candidates, left_exons, right_exons
#' @param neg_exons A GRanges object with candidate negative exons (neg_exons)
#' @param pos_exons A GRanges object with positive exons (pos_exons)
#' @param gr The original GRanges object with all exons (for looking up left/right)
#' @return A named list with three GRanges objects: candidates, left_exons, right_exons 
#' candidates, left_exons, right_exons are all from neg_exons set
#' @importFrom plyranges filter_by_non_overlaps_directed slice
#' @importFrom GenomicRanges GRanges  
candidates_by_non_overlap_directed <- function(neg_exons, pos_exons, gr, type) {
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
  # filter candidates that are internal (have both left and right exons)
  candidates <- candidates |>
    plyranges::filter(internal) #internal col is added in preprocess_input

  # keys of exons to the left and right of candidates
  # TO DO : check if this works for all cases if the strand is - instead of +. ie in a - would the left exon be exon_rank +1?
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1L)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1L)

  # get the actual exons for the candidates (preserves order of keys)
  # exon to the left of the candidates from the neg_exons set
  left_exons  <- gr |> plyranges::slice(match(left_keys, key))
  # exon to the right of the candidates from the neg_exons set
  right_exons <- gr |> plyranges::slice(match(right_keys, key))

  # returns a list of GRanges of same length:
  # ‘candidates’ - neg exons that do not overlap any pos exons, and are internal
  # ‘left_exons’ - left exons of candidates
  # ‘right_exons’ - right exons of candidates
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

    # check for matches to left and right exons in the pos_exons set
    pos_exons_matched <- compute_matches(
      pos_exons,
      left_exon  = left_exon,
      right_exon = right_exon,
      type = type
    )
    # filter and get left/right matched exons (tx_id, exon_rank) as tibbles
    left_tbl <- pos_exons_matched |>
      dplyr::filter(match_left == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(l = exon_rank)

    right_tbl <- pos_exons_matched |>
      dplyr::filter(match_right == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(r = exon_rank)
    # returns a list with two tibbles: left_tbl and right_tbl
    # left_tbl: tx_id, l (exon_rank)
    # right_tbl: tx_id, r (exon_rank)
    # there are as many rows as there are matches between 1 candidate
    # and all the isoforms in pos_exons
    list(
      left_tbl = left_tbl, 
      right_tbl = right_tbl)
}


#' function to find introns given a GRanges object of exons
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coef'.
#' @return A GRanges object with introns as ranges and metadata (tx_id, gene_idß).
find_introns <- function(gr) {
  gr <- gr |> plyranges::arrange(tx_id, start)
   # introns are between exons - use the start of the next exon and the end of the current exon to define their start/end
  gr <- gr |>
    plyranges::group_by(gene_id, tx_id) |>
    plyranges::mutate(
      intron_start = end + 1L,
      intron_end   = dplyr::lead(start) - 1L
    ) |>
    plyranges::filter(!is.na(intron_start) & !is.na(intron_end) & intron_end >= intron_start) |>
    plyranges::ungroup() 
  # create a GRanges object for the introns with the same metadata as the tx_id and gene_id they belong to
  GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr),
      ranges   = IRanges::IRanges(start = gr$intron_start, end = gr$intron_end),
      strand   = GenomicRanges::strand(gr),
      gene_id = gr$gene_id,
      tx_id   = gr$tx_id,
      coefs  = gr$coefs,
      intron = TRUE
    )
  # TO DO : include case where no introns are found eg create_mock_data(1,1,1) and check that the output is an empty GRanges object with the correct metadata columns
}