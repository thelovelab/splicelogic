#' Calculate skipped exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess_input().
#' @param type The type of overlap to consider when identifying skipped exons.
#' @return A GRanges object with an additional 'event' metadata column indicating skipped exons.
#' @export
calc_skipped_exons <- function(gr, type = c("over","in", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  check_preprocessed(gr)

  # separate positive and negative exons
  pos_exons <- gr |> dplyr::filter(sign(estimates) == 1)
  neg_exons <- gr |> dplyr::filter(sign(estimates) == -1)

  # candidates_by_presence_v2 returns GRanges
  filter_results <- candidates_by_presence_v2(gr, neg_exons, pos_exons)

  candidates <- filter_results$candidates
  left_exons <- filter_results$left_exons
  right_exons <- filter_results$right_exons

  if (length(candidates) == 0L) {
    return(GRanges())
  }

  # batch findOverlaps: match pos_exons against all left/right exons at once
  ##############################################################
  left_hits  <- find_matches_batch(pos_exons, left_exons, type) #returns Hits object with queryHits = pos_exons index, subjectHits = left_exons index
  right_hits <- find_matches_batch(pos_exons, right_exons, type)
  pos_tbl <- tibble::as_tibble(pos_exons)
  cand_tbl <- tibble::as_tibble(candidates)
  cand_tbl$cand_idx <- seq_len(nrow(cand_tbl))

  
  # build tibbles from overlap hits: pos_exon index -> (tx_id, exon_rank, cand_idx)
  left_match_tbl <- tibble::tibble(
    pos_idx  = S4Vectors::queryHits(left_hits),
    cand_idx = S4Vectors::subjectHits(left_hits)
  ) |>
    dplyr::mutate(
      tx_id = pos_tbl$tx_id[pos_idx],
      l     = pos_tbl$exon_rank[pos_idx],
      gene_id_pos = pos_tbl$gene_id[pos_idx]
    )

  right_match_tbl <- tibble::tibble(
    pos_idx  = S4Vectors::queryHits(right_hits),
    cand_idx = S4Vectors::subjectHits(right_hits)
  ) |>
    dplyr::mutate(
      tx_id = pos_tbl$tx_id[pos_idx],
      r     = pos_tbl$exon_rank[pos_idx],
      gene_id_pos = pos_tbl$gene_id[pos_idx]
    )

  # restrict matches to same gene as candidate
  left_match_tbl <- left_match_tbl |>
    dplyr::filter(gene_id_pos == cand_tbl$gene_id[cand_idx])
  right_match_tbl <- right_match_tbl |>
    dplyr::filter(gene_id_pos == cand_tbl$gene_id[cand_idx])

  # join left and right by (cand_idx, tx_id), filter adjacent exons
  pairs <- dplyr::inner_join(
    left_match_tbl |> dplyr::select(cand_idx, tx_id, l),
    right_match_tbl |> dplyr::select(cand_idx, tx_id, r),
    by = c("cand_idx", "tx_id")
  ) |>
    dplyr::filter(abs(l - r) == 1) |>
    dplyr::distinct(cand_idx, tx_id)

  if (nrow(pairs) == 0L) {
    return(GRanges())
  }

  # build result tibble: one row per (candidate, tx_event) pair
  hits_tbl <- cand_tbl[pairs$cand_idx, ] |>
    dplyr::mutate(
      event    = "skipped_exon",
      tx_event = pairs$tx_id
    )

  # convert back to GRanges for return 
  GenomicRanges::GRanges(
    seqnames = hits_tbl$seqnames,
    ranges   = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand   = hits_tbl$strand,
    hits_tbl |> dplyr::select(-seqnames, -start, -end, -width, -strand, -cand_idx)
  )
}

#' Calculate mutually exclusive exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess_input().
#' @param type The type of overlap to consider when identifying mutually exclusive exons.
#' @return A GRanges object with an additional 'event' metadata column indicating mutually exclusive exons.
#' @export
calc_mutually_exclusive <- function(gr, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  check_preprocessed(gr)

  # separate positive and negative exons
  pos_exons <- gr |> dplyr::filter(sign(estimates) == 1)
  neg_exons <- gr |> dplyr::filter(sign(estimates) == -1)
  
  # candidates_by_presence_v2 returns GRanges (pairwise-compatible)
  filter_results <- candidates_by_presence_v2(gr, neg_exons, pos_exons)
  candidates <- filter_results$candidates
  left_exons <- filter_results$left_exons
  right_exons <- filter_results$right_exons

    if (length(candidates) == 0L) {
        return(GRanges())
    }

    # batch findOverlaps: match pos_exons against all left/right exons at once
    ##############################################################
    left_hits  <- find_matches_batch(pos_exons, left_exons, type)
    right_hits <- find_matches_batch(pos_exons, right_exons, type)
    pos_tbl <- tibble::as_tibble(pos_exons)
    cand_tbl <- tibble::as_tibble(candidates)
    cand_tbl$cand_idx <- seq_len(nrow(cand_tbl))

    # build tibbles from overlap hits: pos_exon index -> (tx_id, exon_rank, cand_idx)
    left_match_tbl <- tibble::tibble(
        pos_idx  = S4Vectors::queryHits(left_hits),
        cand_idx = S4Vectors::subjectHits(left_hits)
    ) |>
        dplyr::mutate(
            tx_id = pos_tbl$tx_id[pos_idx],
            l     = pos_tbl$exon_rank[pos_idx],
            gene_id_pos = pos_tbl$gene_id[pos_idx]
        )

    right_match_tbl <- tibble::tibble(
        pos_idx  = S4Vectors::queryHits(right_hits),
        cand_idx = S4Vectors::subjectHits(right_hits)
    ) |>
        dplyr::mutate(
            tx_id = pos_tbl$tx_id[pos_idx],
            r     = pos_tbl$exon_rank[pos_idx],
            gene_id_pos = pos_tbl$gene_id[pos_idx]
        )

    # restrict matches to same gene as candidate
    left_match_tbl <- left_match_tbl |>
        dplyr::filter(gene_id_pos == cand_tbl$gene_id[cand_idx])
    right_match_tbl <- right_match_tbl |>
        dplyr::filter(gene_id_pos == cand_tbl$gene_id[cand_idx])

    # join left and right by (cand_idx, tx_id), filter for mx exons (l-r==2)
    pairs <- dplyr::inner_join(
        left_match_tbl |> dplyr::select(cand_idx, tx_id, l),
        right_match_tbl |> dplyr::select(cand_idx, tx_id, r),
        by = c("cand_idx", "tx_id")
    ) |>
        dplyr::filter(abs(l - r) == 2) |>
        dplyr::mutate(middle_rank = pmin(l, r) + 1L)

    if (nrow(pairs) == 0L) {
        return(GRanges())
    }

    # look up the middle pos exon for each pair
    pos_lookup <- pos_tbl |>
        dplyr::mutate(pos_row = dplyr::row_number()) |>
        dplyr::select(tx_id, exon_rank, pos_row)

    pairs <- pairs |>
        dplyr::inner_join(
            pos_lookup,
            by = c("tx_id", "middle_rank" = "exon_rank")
        )

    if (nrow(pairs) == 0L) {
        return(GRanges())
    }

    # build candidate hit rows
    cand_hits <- cand_tbl[pairs$cand_idx, ] |>
        dplyr::mutate(
            event    = "mutually_exclusive",
            tx_event = pairs$tx_id
        )

    # build middle pos exon hit rows
    pos_hits <- pos_tbl[pairs$pos_row, ] |>
        dplyr::mutate(
            event    = "mutually_exclusive",
            tx_event = cand_tbl$tx_id[pairs$cand_idx]
        )

    # interleave candidate and pos hits (cand1, pos1, cand2, pos2, ...)
    hits_tbl <- dplyr::bind_rows(cand_hits, pos_hits) |>
        dplyr::mutate(.pair_order = rep(seq_len(nrow(pairs)), 2)) |>
        dplyr::arrange(.pair_order) |>
        dplyr::select(-.pair_order)

    # convert back to GRanges for return
    GenomicRanges::GRanges(
        seqnames = hits_tbl$seqnames,
        ranges   = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
        strand   = hits_tbl$strand,
        hits_tbl |> dplyr::select(
            -seqnames, -start, -end, -width, -strand,
            -dplyr::any_of("cand_idx")
        )
    )
}

#' Function to calculate retained introns given a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coef'.
#' @return A GRanges object with an additional 'event' metadata column indicating retained introns.
#' @export      
calc_retained_introns <- function(gr){
  # TO DO : fix coef_col argument and add to preprocess_input and then rename coef_col to coefs

  # if preprocessing didn't happen
  check_preprocessed(gr)

  # separate positive and negative exons
  # find introns in the negative coef transcripts
  introns <- gr |> dplyr::filter(estimates < 0) |> find_introns()

  # filter by overlap on the positive coef transcripts
  candidates <- gr |> dplyr::filter(estimates > 0) |> 
                      plyranges::filter_by_overlaps(introns)

  hits <- GRanges()
  if (length(candidates) == 0L) {
      return(hits) #return if no candidates
  }
  for (i in seq_along(candidates)) {
  cand <- candidates[i]  # a length-1 GRanges

  # restric to the same gene
  cand_introns <- introns |> 
                dplyr::filter(gene_id == cand$gene_id)

  # find which transcript pairs the event is happening with 
  txp_events <- cand_introns |> 
              plyranges::filter_by_overlaps(cand)|> 
              tibble::as_tibble()|> 
              dplyr::select(gene_id, tx_id)

  if (nrow(txp_events) > 0) {
    txs <- unique(txp_events$tx_id)
    cand_rep <- rep(cand, length(txs)) |>
      dplyr::mutate(
        event    = "retained_intron",
        tx_event = txs
      )
    hits <- c(hits, cand_rep)
    }
  }
  hits
  }


#' Function to calcualte 5' and 3' alternative splice sites given a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coef'.
#' @return A GRanges object with an additional 'event' metadata column indicating retained introns.
#' @export      
calc_a3ss_a5ss <- function(gr ){
  # if preprocessing didn't happen
  check_preprocessed(gr)

  # separate positive and negative exons
  pos_exons <- gr |> dplyr::filter(sign(estimates) == 1)
  neg_exons <- gr |> dplyr::filter(sign(estimates) == -1)

  # candidates are pos exons that do not exacly match neg_exons (%in%)
  # filter candidates pos_exons that are exactly the same as any neg_exons
  # and keep only those that overlap any neg_exons
  # TO DO: this is not pairwise comparison - need to restrict to the same gene and then compare to all neg exons in that gene pair by pair
  # rn if any neg exon is the same as the candidate, it will be filtered out, even if there is another neg exon that overlaps but is not the same and could be a valid a3ss/a5ss event.
  candidates <- pos_exons |> 
                dplyr::filter(!(pos_exons %in% neg_exons)) |>
                plyranges::filter_by_overlaps_directed(neg_exons)

  hits <- GRanges()
  if (length(candidates) == 0L) {
    return(hits) #return unchanged if no candidates
  }

 for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges
    # restric to check on the same gene in neg_exons
    cand_neg_exons <- neg_exons |> 
                  dplyr::filter(gene_id == cand$gene_id) 
                  # dplyr::filter(!(. %in%c cand))
                  # dplyr::filter_by_overlaps(cand)


      matches <- cand_neg_exons %>%
        dplyr::mutate(
          match_start  = GenomicRanges::start(.) %in% GenomicRanges::start(cand),
          match_end = GenomicRanges::end(.) %in% GenomicRanges::end(cand)
        ) |>
         # only keep those that match on one end but not the other
        dplyr::filter(match_start != match_end)
        
      if (length(matches) > 0L) {
        event_type <- ifelse(matches$match_start, "a5ss", "a3ss")
        txs <- matches$tx_id
        cand_rep <- rep(cand, length(txs)) |>
          dplyr::mutate(
            event    = event_type,
            tx_event = txs
          )
        hits <- c(hits, cand_rep)
        }
  }
  hits
}
