#' @rdname find_events
#' @name find_events
#'
#' @title Find splice events from annotated exons
#'
#' @description Functions to find different types of alternative splicing
#' events from preprocessed GRanges exon data. Events include skipped exon (se),
#' included exon (ie), mutatualy exclusive exons (mxe), retained intron (ri),
#' and alternative 5' and 3' splice sites (a5ss / a3ss).
#'
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess().

NULL

#' @rdname find_events
#' @param type The type of overlap to consider when
#' identifying events.
#' @param inverse If TRUE, identifies included exons
#' instead of skipped exons.
#' @return A GRanges object with an additional column `event` indicating:
#' 
#' `find_se()`: skipped exons
#' @export
#' @examples
#'
#' # make some mock data and run the function
#' gr <- create_mock_data(n_genes = 2, n_tx = 4, n_exons = 4) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_skipped_exons(1)
#'
#' # this should find the skipped exon events we generated
#' find_se(gr, type = "boundary")
#'
#' find_ie(gr, type = "boundary")
#'
find_se <- function(
  gr,
  type = c("boundary", "over", "in"),
  inverse = FALSE
) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  check_preprocessed(gr)
  # keep user's original metadata columns (minus preprocessing intermediates)
  keep_cols <- keep_cols(gr)
  if (inverse) {
    factor <- -1
    event_name <- "ie" # included exon
  } else {
    factor <- 1
    event_name <- "se" # skipped exon
  }

  matches <- find_candidates_and_flanks(gr, type, factor)
  #check if early return is needed
  if (length(matches) == 0L) {
    return(GenomicRanges::GRanges())
  }
  #unlist results
  pos_tbl <- matches$pos_tbl
  cand_tbl <- matches$cand_tbl
  left_match_tbl <- matches$left_match_tbl
  right_match_tbl <- matches$right_match_tbl

  # join left and right by (cand_idx, tx_id), filter adjacent exons
  pairs <- dplyr::inner_join(
    left_match_tbl |> dplyr::select(cand_idx, tx_id, l),
    right_match_tbl |> dplyr::select(cand_idx, tx_id, r),
    by = c("cand_idx", "tx_id")
  ) |>
    # for SE/included exon, flanking exons must be
    # adjacent (rank difference of 1)
    dplyr::filter(abs(l - r) == 1) |>
    dplyr::distinct(cand_idx, tx_id)

  if (nrow(pairs) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build result tibble: one row per (candidate, tx_event) pair
  hits_tbl <- cand_tbl[pairs$cand_idx, ] |>
    dplyr::mutate(
      event = event_name,
      tx_event = pairs$tx_id
    )
  # convert back to GRanges for return
  res <- GenomicRanges::GRanges(
    seqnames = hits_tbl$seqnames,
    ranges = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand = hits_tbl$strand,
    hits_tbl |> dplyr::select(dplyr::all_of(keep_cols), event, tx_event)
  )
  # preserve seqinfo from input
  GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
  res
}

#' @rdname find_events
#' @return `find_ie()`: included exons
#' @export
find_ie <- function(gr, type = c("boundary", "over", "in")) {
  find_se(gr, type, inverse = TRUE)
}

#' @rdname find_events
#' @return `find_mxe()`: mutually exclusive exons
#' @examples
#'
#' # detect mutually exclusive exons
#' gr_mx <- create_mock_data(
#'   n_genes = 2, n_tx = 4, n_exons = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_mx(1)
#' 
#' find_mxe(gr_mx, type = "boundary")
#'
#' @export
find_mxe <- function(gr, type = c("boundary", "in", "over")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  check_preprocessed(gr)
  keep_cols <- keep_cols(gr)

  matches <- find_candidates_and_flanks(gr, type, 1)

  #check if early return is needed
  if (length(matches) == 0L) {
    return(GenomicRanges::GRanges())
  }
  #unlist results
  pos_tbl <- matches$pos_tbl
  cand_tbl <- matches$cand_tbl
  left_match_tbl <- matches$left_match_tbl
  right_match_tbl <- matches$right_match_tbl

  # join left and right by (cand_idx, tx_id), filter for mx exons (l-r==2)
  pairs <- dplyr::inner_join(
    left_match_tbl |> dplyr::select(cand_idx, tx_id, l),
    right_match_tbl |> dplyr::select(cand_idx, tx_id, r),
    by = c("cand_idx", "tx_id")
  ) |>
    # for MX, flanking exons must have a gap of
    # exactly 2 (one exon in between)
    dplyr::filter(abs(l - r) == 2) |>
    dplyr::mutate(middle_rank = pmin(l, r) + 1L)

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
    return(GenomicRanges::GRanges())
  }

  # build candidate hit rows
  cand_hits <- cand_tbl[pairs$cand_idx, ] |>
    dplyr::mutate(
      event = "mxe", # mutually exclusive exon
      tx_event = pairs$tx_id
    )

  # build middle pos exon hit rows
  pos_hits <- pos_tbl[pairs$pos_row, ] |>
    dplyr::mutate(
      event = "mxe", # mutually exclusive exon
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
    ranges = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand = hits_tbl$strand,
    hits_tbl |> dplyr::select(dplyr::all_of(keep_cols), event, tx_event)
  )
  # preserve seqinfo from input
  GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
  res
}

#' @rdname find_events
#' @return `find_ri()`: retained introns
#' @examples
#'
#' # detect retained introns
#' gr_ri <- create_mock_data(
#'   n_genes = 2, n_tx = 4, n_exons = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_retained_introns(1)
#' 
#' find_ri(gr_ri)
#'
#' @export
find_ri <- function(gr) {
  # if preprocessing didn't happen
  check_preprocessed(gr)
  keep_cols <- keep_cols(gr)

  # find introns in the negative coef transcripts
  neg_exons <- gr |> dplyr::filter(sign(estimate) == -1)
  pos_exons <- gr |> dplyr::filter(sign(estimate) == 1)
  introns <- neg_exons |> find_introns()

  if (length(introns) == 0L || length(pos_exons) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # batch findOverlaps: intron must be fully within the pos exon
  hits <- GenomicRanges::findOverlaps(introns, pos_exons, type = "within")

  if (length(hits) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build match tibble from hits
  intron_tbl <- tibble::as_tibble(introns)
  pos_tbl <- tibble::as_tibble(pos_exons)

  match_tbl <- tibble::tibble(
    intron_idx = S4Vectors::queryHits(hits),
    pos_idx = S4Vectors::subjectHits(hits)
  ) |>
    dplyr::mutate(
      gene_id_intron = intron_tbl$gene_id[intron_idx],
      gene_id_pos = pos_tbl$gene_id[pos_idx],
      tx_id_intron = intron_tbl$tx_id[intron_idx]
    ) |>
    # restrict to same gene
    dplyr::filter(gene_id_intron == gene_id_pos) |>
    dplyr::distinct(pos_idx, tx_id_intron)

  if (nrow(match_tbl) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build result: one row per (pos exon, intron transcript) pair
  hits_tbl <- pos_tbl[match_tbl$pos_idx, ] |>
    dplyr::mutate(
      event = "ri", # retained intron
      tx_event = match_tbl$tx_id_intron
    )

  # convert back to GRanges for return
  res <- GenomicRanges::GRanges(
    seqnames = hits_tbl$seqnames,
    ranges = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand = hits_tbl$strand,
    hits_tbl |> dplyr::select(dplyr::all_of(keep_cols), event, tx_event)
  )
  # preserve seqinfo from input
  GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
  res
}


#' @rdname find_events
#' @param by_start If TRUE, detects a5ss (same exon start, different end).
#' If FALSE, detects a3ss (same end, different start).
#' @return `find_alt_ss()`: alternative splice sites
#' @noRd 
find_alt_ss <- function(gr, by_start = TRUE) {
  # if preprocessing didn't happen
  check_preprocessed(gr)
  keep_cols <- keep_cols(gr)

  event_name <- if (by_start) "a5ss" else "a3ss"

  # separate positive and negative exons
  pos_exons <- gr |> dplyr::filter(sign(estimate) == 1)
  neg_exons <- gr |> dplyr::filter(sign(estimate) == -1)

  # candidates: pos exons not exactly identical to any neg exon,
  # but overlapping at least one neg exon (directed)
  candidates <- pos_exons |>
    dplyr::filter(!(pos_exons %in% neg_exons)) |>
    plyranges::filter_by_overlaps_directed(neg_exons)

  if (length(candidates) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # batch findOverlaps: candidates against neg_exons
  hits <- GenomicRanges::findOverlaps(candidates, neg_exons)

  if (length(hits) == 0L) {
    return(GenomicRanges::GRanges())
  }

  cand_tbl <- tibble::as_tibble(candidates)
  neg_tbl <- tibble::as_tibble(neg_exons)

  # build match tibble with boundary checks
  match_tbl <- tibble::tibble(
    cand_idx = S4Vectors::queryHits(hits),
    neg_idx = S4Vectors::subjectHits(hits)
  ) |>
    dplyr::mutate(
      gene_id_cand = cand_tbl$gene_id[cand_idx],
      gene_id_neg = neg_tbl$gene_id[neg_idx],
      match_start = cand_tbl$start[cand_idx] == neg_tbl$start[neg_idx],
      match_end = cand_tbl$end[cand_idx] == neg_tbl$end[neg_idx],
      tx_id_neg = neg_tbl$tx_id[neg_idx]
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
    return(GenomicRanges::GRanges())
  }

  # build result: one row per (candidate, neg exon match)
  hits_tbl <- cand_tbl[match_tbl$cand_idx, ] |>
    dplyr::mutate(
      event = match_tbl$event,
      tx_event = match_tbl$tx_id_neg
    )

  # convert back to GRanges for return
  res <- GenomicRanges::GRanges(
    seqnames = hits_tbl$seqnames,
    ranges = IRanges::IRanges(start = hits_tbl$start, end = hits_tbl$end),
    strand = hits_tbl$strand,
    hits_tbl |> dplyr::select(dplyr::all_of(keep_cols), event, tx_event)
  )
  # preserve seqinfo from input
  GenomicRanges::seqinfo(res) <- GenomicRanges::seqinfo(gr)
  res
}


#' @rdname find_events
#' @return `find_a5ss()`: alternative 5' splice sites
#' @examples
#'
#' # detect alternative 5' splice sites
#' gr_a5 <- create_mock_data(
#'   n_genes = 2, n_tx = 4, n_exons = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_a5ss(1)
#' 
#' find_a5ss(gr_a5)
#'
#' @export
find_a5ss <- function(gr) {
  find_alt_ss(gr, by_start = TRUE)
}

#' @rdname find_events
#' @return `find_a3ss()`: : alternative 3' splice sites
#' @examples
#'
#' # detect alternative 3' splice sites
#' gr_a3 <- create_mock_data(
#'   n_genes = 2, n_tx = 4, n_exons = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_a3ss(1)
#' find_a3ss(gr_a3)
#'
#' @export
find_a3ss <- function(gr) {
  find_alt_ss(gr, by_start = FALSE)
}

#' @rdname find_events
#' @param verbose If TRUE, prints progress messages. Default TRUE.
#' @return `find_all_events()`: all detected events
#' 
#' @examples
#'
#' # detect all event types at once
#' gr_all <- create_mock_data(
#'   n_genes = 2, n_tx = 4, n_exons = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_skipped_exons(1)
#' 
#' find_all_events(gr_all, type = "boundary", verbose = FALSE)
#'
#' @export
find_all_events <- function(
  gr,
  type = c("boundary", "over", "in"),
  verbose = TRUE
) {
  type <- match.arg(type)
  check_preprocessed(gr)
  msg <- if (verbose) message else function(...) invisible (NULL)

  results <- list()

  msg("Calculating skipped exon events...")
  results$se <- find_se(gr, type)

  msg("Calculating included exon events...")
  results$ie <- find_ie(gr, type)

  msg("Calculating mutually exclusive exon events...")
  results$mx <- find_mxe(gr, type)

  msg("Calculating retained intron events...")
  results$ri <- find_ri(gr)

  msg("Calculating alternative 5' splice site events...")
  results$a5ss <- find_a5ss(gr)

  msg("Calculating alternative 3' splice site events...")
  results$a3ss <- find_a3ss(gr)

  # keep only non-empty results
  results <- Filter(function(x) length(x) > 0L, results)

  if (length(results) == 0L) {
    msg("No events detected.")
    return(GenomicRanges::GRanges())
  }

  combined <- do.call(plyranges::bind_ranges, results)
  msg("Done! ", length(combined), " events detected.")
  combined
}
