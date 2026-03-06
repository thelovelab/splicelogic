#' Calculate skipped exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess_input().
#' @param type The type of overlap to consider when identifying skipped exons.
#' @param inverse If TRUE, identifies included exons instead of skipped exons.
#' @return A GRanges object with an additional 'event' metadata column indicating skipped exons.
#' @export
calc_skipped_exons <- function(gr, type = c("boundary","over","in"), inverse = FALSE) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  check_preprocessed(gr)

  if (inverse) {
    gr <- gr |> dplyr::mutate(estimates = -estimates)
    event_name <- "included_exon"
  } else {
    event_name <- "skipped_exon"
  }

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
      event    = event_name,
      tx_event = pairs$tx_id
    )

  # convert back to GRanges for return
  res <- GenomicRanges::GRanges(
    seqnames = hits_tbl$seqnames,
    ranges   = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand   = hits_tbl$strand,
    hits_tbl |> dplyr::select(gene_id, tx_id, exon_rank, estimates, event, tx_event)
  )
  # preserve seqinfo from input
  GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
  res
}

#' Calculate included exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess_input().
#' @param type The type of overlap to consider when identifying included exons.
#' @return A GRanges object with an additional 'event' metadata column indicating included exons.
#' @export  
calc_included_exons <- function(gr, type = c("boundary","over","in")) {
  calc_skipped_exons(gr, type, inverse = TRUE)
}

#' Calculate mutually exclusive exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess_input().
#' @param type The type of overlap to consider when identifying mutually exclusive exons.
#' @return A GRanges object with an additional 'event' metadata column indicating mutually exclusive exons.
#' @export
calc_mx_exons <- function(gr, type = c("boundary","in", "over")) {
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
        ) |>
        # exclude pairs where the middle pos exon has the same coordinates
        # as the candidate neg exon (not a true MX — just an overlap/SE)
        # FP scenario: candidate neg exon (e.g. 21-25) is absent from some pos
        # txps but present in another (tx3). tx3 has flanking matches with gap=2,
        # and the middle pos exon IS the candidate itself — not a true MX pair.
        dplyr::filter(
            pos_tbl$start[pos_row] != cand_tbl$start[cand_idx] |
            pos_tbl$end[pos_row] != cand_tbl$end[cand_idx]
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
    res <- GenomicRanges::GRanges(
        seqnames = hits_tbl$seqnames,
        ranges   = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
        strand   = hits_tbl$strand,
        hits_tbl |> dplyr::select(gene_id, tx_id, exon_rank, estimates, event, tx_event)
    )
    # preserve seqinfo from input
    GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
    res
}

#' Function to calculate retained introns given a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coef'.
#' @return A GRanges object with an additional 'event' metadata column indicating retained introns.
#' @export      
calc_retained_introns <- function(gr){
  # if preprocessing didn't happen
  check_preprocessed(gr)

    # find introns in the negative coef transcripts
    neg_exons <- gr |> dplyr::filter(sign(estimates) == -1)
    pos_exons <- gr |> dplyr::filter(sign(estimates) == 1)
    introns <- neg_exons |> find_introns()

    if (length(introns) == 0L || length(pos_exons) == 0L) {
        return(GRanges())
    }

    # batch findOverlaps: intron must be fully within the pos exon
    hits <- GenomicRanges::findOverlaps(introns, pos_exons,
                                        type = "within")

    if (length(hits) == 0L) {
        return(GRanges())
    }

    # build match tibble from hits
    intron_tbl <- tibble::as_tibble(introns)
    pos_tbl <- tibble::as_tibble(pos_exons)

    match_tbl <- tibble::tibble(
        intron_idx = S4Vectors::queryHits(hits),
        pos_idx    = S4Vectors::subjectHits(hits)
    ) |>
        dplyr::mutate(
            gene_id_intron = intron_tbl$gene_id[intron_idx],
            gene_id_pos    = pos_tbl$gene_id[pos_idx],
            tx_id_intron   = intron_tbl$tx_id[intron_idx]
        ) |>
        # restrict to same gene
        dplyr::filter(gene_id_intron == gene_id_pos) |>
        dplyr::distinct(pos_idx, tx_id_intron)

    if (nrow(match_tbl) == 0L) {
        return(GRanges())
    }

    # build result: one row per (pos exon, intron transcript) pair
    hits_tbl <- pos_tbl[match_tbl$pos_idx, ] |>
        dplyr::mutate(
            event    = "retained_intron",
            tx_event = match_tbl$tx_id_intron
        )

    # convert back to GRanges for return
    res <- GenomicRanges::GRanges(
        seqnames = hits_tbl$seqnames,
        ranges   = IRanges::IRanges(start = hits_tbl$start,
                                     end = hits_tbl$end),
        strand   = hits_tbl$strand,
        hits_tbl |> dplyr::select(gene_id, tx_id, exon_rank, estimates, event, tx_event)
    )
    # preserve seqinfo from input
    GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
    res
}


#' Function to calculate alternative splice sites
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id',
#' 'tx_id', and 'coef'. Must be preprocessed with preprocess_input().
#' @param by_start If TRUE, detects a5ss (same exon start, different end).
#' If FALSE, detects a3ss (same end, different start).
#' @return A GRanges object with annotated events for alternative splice sites.
#' @keywords internal
calc_alt_ss <- function(gr, by_start = TRUE) {
    # if preprocessing didn't happen
    check_preprocessed(gr)

    event_name <- if (by_start) "a5ss" else "a3ss"

    # separate positive and negative exons
    pos_exons <- gr |> dplyr::filter(sign(estimates) == 1)
    neg_exons <- gr |> dplyr::filter(sign(estimates) == -1)

    # candidates: pos exons not exactly identical to any neg exon,
    # but overlapping at least one neg exon (directed)
    candidates <- pos_exons |>
        dplyr::filter(!(pos_exons %in% neg_exons)) |>
        plyranges::filter_by_overlaps_directed(neg_exons)

    if (length(candidates) == 0L) {
        return(GRanges())
    }

    # batch findOverlaps: candidates against neg_exons
    hits <- GenomicRanges::findOverlaps(candidates, neg_exons)

    if (length(hits) == 0L) {
        return(GRanges())
    }

    cand_tbl <- tibble::as_tibble(candidates)
    neg_tbl <- tibble::as_tibble(neg_exons)

    # build match tibble with boundary checks
    match_tbl <- tibble::tibble(
        cand_idx = S4Vectors::queryHits(hits),
        neg_idx  = S4Vectors::subjectHits(hits)
    ) |>
        dplyr::mutate(
            gene_id_cand = cand_tbl$gene_id[cand_idx],
            gene_id_neg  = neg_tbl$gene_id[neg_idx],
            match_start  = cand_tbl$start[cand_idx] ==
                           neg_tbl$start[neg_idx],
            match_end    = cand_tbl$end[cand_idx] ==
                           neg_tbl$end[neg_idx],
            tx_id_neg    = neg_tbl$tx_id[neg_idx]
        ) |>
        # restrict to same gene
        dplyr::filter(gene_id_cand == gene_id_neg)

    # exclude candidates that overlap 2+ neg exons in the same neg
    # transcript — those are retained introns, not alternative splice sites
    multi_overlap <- match_tbl |>
        dplyr::count(cand_idx, tx_id_neg) |>
        dplyr::filter(n >= 2)

    match_tbl <- match_tbl |>
        dplyr::anti_join(multi_overlap, by = c("cand_idx", "tx_id_neg")) |>
        # one boundary matches but not the other (XOR)
        dplyr::filter(match_start != match_end) |>
        dplyr::filter(match_start == by_start) |>
        dplyr::mutate(
            event = event_name
        )

    if (nrow(match_tbl) == 0L) {
        return(GRanges())
    }

    # build result: one row per (candidate, neg exon match)
    hits_tbl <- cand_tbl[match_tbl$cand_idx, ] |>
        dplyr::mutate(
            event    = match_tbl$event,
            tx_event = match_tbl$tx_id_neg
        )

    # convert back to GRanges for return
    res <- GenomicRanges::GRanges(
        seqnames = hits_tbl$seqnames,
        ranges   = IRanges::IRanges(start = hits_tbl$start,
                                     end = hits_tbl$end),
        strand   = hits_tbl$strand,
        hits_tbl |> dplyr::select(gene_id, tx_id, exon_rank, estimates, event, tx_event)
    )
    # preserve seqinfo from input
    GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
    res

}


#' Calculate alternative 5' splice sites from a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id',
#' 'tx_id', and 'coef'. Must be preprocessed with preprocess_input().
#' @return A GRanges object with an additional 'event' metadata column
#' indicating a5ss events.
#' @export
calc_a5ss <- function(gr) {
    calc_alt_ss(gr, by_start = TRUE)
}

#' Calculate alternative 3' splice sites from a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id',
#' 'tx_id', and 'coef'. Must be preprocessed with preprocess_input().
#' @return A GRanges object with an additional 'event' metadata column
#' indicating a3ss events.
#' @export
calc_a3ss <- function(gr) {
    calc_alt_ss(gr, by_start = FALSE)
}

#' Calculate all splicing events from a GRanges object
#' @description Runs all event detection functions and returns a single
#' concatenated GRanges with results from each.
#' @param gr A GRanges object with exon annotations, preprocessed with
#' preprocess_input().
#' @param type The type of overlap to consider for skipped exons, included
#' exons, and mutually exclusive exons.
#' @return A GRanges object combining all detected events, with an 'event'
#' metadata column indicating the event type.
#' @export
calc_all_events <- function(gr, type = c("boundary", "over", "in")) {
    type <- match.arg(type)
    check_preprocessed(gr)

    results <- list()

    message("Calculating skipped exon events...")
    results$se <- calc_skipped_exons(gr, type)

    message("Calculating included exon events...")
    results$ie <- calc_included_exons(gr, type)

    message("Calculating mutually exclusive exon events...")
    results$mx <- calc_mx_exons(gr, type)

    message("Calculating retained intron events...")
    results$ri <- calc_retained_introns(gr)

    message("Calculating alternative 5' splice site events...")
    results$a5ss <- calc_a5ss(gr)

    message("Calculating alternative 3' splice site events...")
    results$a3ss <- calc_a3ss(gr)

    # keep only non-empty results
    results <- Filter(function(x) length(x) > 0L, results)

    if (length(results) == 0L) {
    message("No events detected.")
    return(GenomicRanges::GRanges())
    }

    combined <- do.call(plyranges::bind_ranges, results)
    message("Done! ", length(combined), " events detected.")
    combined
}