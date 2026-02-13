#' Takes a GRanges object and left/right exons to compute matches.
#' i.e takes the pos_exons (GRanges) and left and right exons and returns
#' a tibble that is the same as pos_exons but with two additional columns
#' 'match_left' and 'match_right' indicating whether each exon matches the left
#' and right exons respectively. The matching is done based on the 'type'
#' parameter which can be "in", "over", or "boundary".
#' @param gr A GRanges object with exon annotations
#' @param left_exon A GRanges object with left exon(s) to match
#' @param right_exon A GRanges object with right exon(s) to match
#' @param type The type of overlap to consider when identifying matches.
#' @return A tibble created from gr with two additional columns:
#' 'match_left' and 'match_right' indicating the number of overlaps with
#' the left and right exons, respectively.
#' @importFrom magrittr %>%
#' @keywords internal
compute_matches <- function(gr, left_exon, right_exon, type = c("over","in","boundary")) {
  type <- match.arg(type)
  gr_matched <- switch(
    type,
    "over" = {
        gr %>%
        dplyr::mutate(
            match_left  = plyranges::count_overlaps(., left_exon) >= 1,
            match_right = plyranges::count_overlaps(., right_exon) >= 1
        ) |>
        tibble::as_tibble()
    },
    "in" = {
      gr %>%
        dplyr::mutate(
          match_left  = . == left_exon,
          match_right = . == right_exon
        ) |>
        tibble::as_tibble()
    },
    "boundary" = {
      gr %>%
        dplyr::mutate(
          match_left  = GenomicRanges::end(.)   %in% GenomicRanges::end(left_exon),
          match_right = GenomicRanges::start(.) %in% GenomicRanges::start(right_exon)
        ) |>
        tibble::as_tibble()
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
#' @param type The type of overlap to consider when identifying matches.
#' @return A named list with three GRanges objects: candidates, left_exons, right_exons 
#' candidates, left_exons, right_exons are all from neg_exons set
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
    dplyr::filter(internal) #internal col is added in preprocess_input

  # keys of exons to the left and right of candidates
  # TO DO : check if this works for all cases if the strand is - instead of +. ie in a - would the left exon be exon_rank +1?
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1L)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1L)

  # get the actual exons for the candidates (preserves order of keys)
  # exon to the left of the candidates from the neg_exons set
  left_exons  <- gr |> dplyr::slice(match(left_keys, key))
  # exon to the right of the candidates from the neg_exons set
  right_exons <- gr |> dplyr::slice(match(right_keys, key))

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
#' @importFrom magrittr %>%
#' @keywords internal
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
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coefs'.
#' @return A GRanges object with introns as ranges and metadata (tx_id, gene_id)
#' @keywords internal
find_introns <- function(gr) {
  gr <- gr |> dplyr::arrange(tx_id, start)
  gr <- gr |> dplyr::arrange(tx_id, start)
   # introns are between exons - use the start of the next exon and the end of the current exon to define their start/end
  gr <- gr |>
    dplyr::group_by(gene_id, tx_id) |>
    dplyr::mutate(
      intron_start = end + 1L,
      intron_end   = dplyr::lead(start) - 1L
    ) |>
    dplyr::filter(!is.na(intron_start) & !is.na(intron_end) & intron_end >= intron_start) |>
    dplyr::ungroup() 
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


#' Filter candidates based on their presence in transcripts
#' Then get the left and right exons for each candidate
#' Return a named list with three GRanges objects: candidates, left_exons, right
#' @param gr A GRanges object with all exons
#' @param coef_col The name of the coefficient column in gr
#' @return A named list with three GRanges objects: candidates, left_exons, right_exons 
#' candidates, left_exons, right_exons are all from neg_exons set
#' @keywords internal
candidates_by_presence <- function(gr, coef_col) {
  coef <- rlang::sym(coef_col)

  # filter internal exons and put n_txp per gene
  gr_internal <- gr |> 
    dplyr::group_by(gene_id) |>
    dplyr::mutate(n_txp = dplyr::n_distinct(tx_id)) |>
    dplyr::ungroup()

  # full data
  keys <- data.frame(seqnames = GenomicRanges::seqnames(gr_internal), 
                      start = IRanges::start(gr_internal),
                      end = IRanges::end(gr_internal),
                      gene_id = gr_internal$gene_id,
                      n_txp = gr_internal$n_txp)
                      
  # count occurrences of each exon 
  key_count <- vctrs::vec_count(keys) |>
    dplyr::filter((count / key$n_txp) < 1)

  # match it back to GRanges
   candidates <- gr_internal[vctrs::vec_in(keys, key_count$key)] |>
    dplyr::filter(internal & (!!coef < 0 )) 
    
  # early return: no candidates -> empty list with GRanges objects
  if (length(candidates) == 0L) {
    return(list(
      candidates = candidates,
      left_exons  = GenomicRanges::GRanges(), 
      right_exons = GenomicRanges::GRanges()
    ))
  }
  # keys of exons to the left and right of candidates
  # TO DO : check if this works for all cases if the strand is - instead of +. ie in a - would the left exon be exon_rank +1?
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1L)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1L)

  # get the actual exons for the candidates (preserves order of keys)
  # exon to the left of the candidates from the neg_exons set
  left_exons  <- gr |> dplyr::slice(match(left_keys, key))
  # exon to the right of the candidates from the neg_exons set
  right_exons <- gr |> dplyr::slice(match(right_keys, key))

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

candidates_by_presence_v2 <- function(neg_exons, pos_exons) {
  pos_exons <- pos_exons |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(n_txp_pos = dplyr::n_distinct(tx_id)) |>
    dplyr::ungroup()
  count <- neg_exons |> plyranges::count_overlaps(pos_exons) 

  neg_exons <- neg_exons |>
    dplyr::mutate(overlap_count = count,
                n_txp_pos = pos_exons$n_txp_pos[match(gene_id, pos_exons$gene_id)]) 
  candidates <- neg_exons |>
    dplyr::filter(internal & (overlap_count / n_txp_pos) < 1)

 # early return: no candidates -> empty list with GRanges objects
  if (length(candidates) == 0L) {
    return(list(
      candidates = candidates,
      left_exons  = GenomicRanges::GRanges(), 
      right_exons = GenomicRanges::GRanges()
    ))
  }
  # keys of exons to the left and right of candidates
  # TO DO : check if this works for all cases if the strand is - instead of +. ie in a - would the left exon be exon_rank +1?
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1L)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1L)

  # get the actual exons for the candidates (preserves order of keys)
  # exon to the left of the candidates from the neg_exons set
  left_exons  <- gr |> dplyr::slice(match(left_keys, key))
  # exon to the right of the candidates from the neg_exons set
  right_exons <- gr |> dplyr::slice(match(right_keys, key))

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
