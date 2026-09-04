#' @rdname find_events
#' @name find_events
#'
#' @title Find splice events from annotated exons
#'
#' @description Functions to find different types of alternative splicing
#' events from preprocessed GRanges exon data. Events include skipped exon (se),
#' included exon (ie), mutatualy exclusive exons (mxe), retained intron (ri),
#' alternative 5' and 3' splice sites (a5ss / a3ss), and alternative
#' transcription start and end sites (aTSS / aTES).
#'
#' `find_a5ss()` / `find_a3ss()` do not require either boundary of an exon
#' to be shared with its partner, so an exon that moved both boundaries is
#' reported as an a5ss *and* an a3ss. A first exon has no acceptor and a
#' last exon has no donor, so those boundaries never produce an a3ss /
#' a5ss; `find_aTSS()` / `find_aTES()` report them instead by comparing
#' where each transcript begins and ends. `exon_rank` runs in
#' transcription order, so this holds on both strands.
#'
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns and preprocessed with preprocess().
#'
#' @return A GRanges object with the detected exon ranges and the following
#' additional metadata columns:
#' \describe{
#'   \item{\code{event_type}}{The type of splicing event detected (e.g.
#'     \code{"se"}, \code{"ie"}, \code{"mxe"}, \code{"ri"}, \code{"a5ss"},
#'     \code{"a3ss"}, \code{"aTSS"}, \code{"aTES"}).}
#'   \item{\code{event_tx_id}}{Transcript ID of the paired transcript
#'     involved in the event.}
#'   \item{\code{event_estimate}}{DTU coefficient of the paired transcript.}
#'   \item{\code{event_<col>}}{One column per name in
#'     \code{metadata(gr)$additional_columns}, prefixed with \code{event_},
#'     carrying the corresponding value from the paired transcript.}
#' }

NULL

# for find_se
utils::globalVariables(c("cand_idx", "l", "r", "event_type", "event_tx_id"))

#' @rdname find_events
#' @param type The type of overlap to consider when
#' identifying events.
#' @param inverse If TRUE, identifies included exons
#' instead of skipped exons.
#' @return `find_se()`: skipped exons
#' @export
#' @examples
#'
#' # make some mock data and run the function
#' gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_se(n_events = 1)
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
    dplyr::mutate(event_type = event_name, event_tx_id = pairs$tx_id)
  # convert back to GRanges for return
  tbl_to_granges(hits_tbl, keep_cols, gr)
}

#' @rdname find_events
#' @export
find_skipped_exons <- find_se

#' @rdname find_events
#' @return `find_ie()`: included exons
#' @export
find_ie <- function(gr, type = c("boundary", "over", "in")) {
  find_se(gr, type, inverse = TRUE)
}


#' @rdname find_events
#' @export
find_included_exons <- find_ie

# for find_mxe
utils::globalVariables(c(
  "cand_idx", "l", "r", "pos_row", ".pair_order", "event_type", "event_tx_id"
))

#' @rdname find_events
#' @return `find_mxe()`: mutually exclusive exons
#' @examples
#'
#' # detect mutually exclusive exons
#' gr_mx <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_mxe(n_events = 1)
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
    # exclude pairs where the middle pos exon overlaps the candidate neg
    # exon at all — mutually exclusive exons are disjoint by definition, so
    # the middle exon must lie entirely before or after the candidate.
    # false positive 1: candidate neg exon (e.g. 21-25) is absent from some
    # pos txps but present in another (tx3). tx3 has flanking matches with
    # gap=2, and the middle pos exon IS the candidate — not a true MX pair.
    # false positive 2: the middle pos exon is the candidate with one boundary
    # shifted (e.g. 21-25 vs 23-25). Truncating an exon does not change how
    # many exons sit between the flanks, so the gap stays 2 and an
    # alternative 5'/3' splice site would otherwise be reported as MX too.
    # similar to filter by non overlapping middle pos exon
    dplyr::filter(
      pos_tbl$end[pos_row] < cand_tbl$start[cand_idx] |
        pos_tbl$start[pos_row] > cand_tbl$end[cand_idx]
    )

  if (nrow(pairs) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build candidate hit rows
  cand_hits <- cand_tbl[pairs$cand_idx, ] |>
    dplyr::mutate(
      event_type = "mxe", # mutually exclusive exon
      event_tx_id = pairs$tx_id
    )

  # build middle pos exon hit rows
  pos_hits <- pos_tbl[pairs$pos_row, ] |>
    dplyr::mutate(
      event_type = "mxe", # mutually exclusive exon
      event_tx_id = cand_tbl$tx_id[pairs$cand_idx]
    )
  # interleave candidate and pos hits (cand1, pos1, cand2, pos2, ...)
  hits_tbl <- dplyr::bind_rows(cand_hits, pos_hits) |>
    dplyr::mutate(.pair_order = rep(seq_len(nrow(pairs)), 2)) |>
    dplyr::arrange(.pair_order) |>
    dplyr::select(-.pair_order)
  # convert back to GRanges for return
  tbl_to_granges(hits_tbl, keep_cols, gr)
}

#' @rdname find_events
#' @export
find_mutually_exclusive_exons <- find_mxe

# for find_ri
utils::globalVariables(c(
  "estimate", "intron_idx", "pos_idx", "gene_id_intron",
  "gene_id_pos", "tx_id_intron", "event_type", "event_tx_id"
))

#' @rdname find_events
#' @return `find_ri()`: retained introns
#' @examples
#'
#' # detect retained introns
#' gr_ri <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_ri(n_events = 1)
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
      event_type = "ri", # retained intron
      event_tx_id = match_tbl$tx_id_intron
    )

  # convert back to GRanges for return
  tbl_to_granges(hits_tbl, keep_cols, gr)
}

#' @rdname find_events
#' @export
find_retained_introns <- find_ri

# for find_alt_ss
utils::globalVariables(c(
  "estimate", "cand_idx", "neg_idx", "gene_id_cand", "gene_id_neg",
  "tx_id_neg", "n", "match_start", "match_end", "event_type", "event_tx_id",
  "on_minus", "is_first", "is_last", "site_is_start", "site_matches"
))

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

  # flag each transcript's first and last exon. a first exon has no
  # acceptor and a last exon has no donor, so a boundary moving there is
  # a transcription start / end change, not a splice site change --
  # find_aTSS() / find_aTES() report those instead. exon_rank follows
  # transcription, so rank 1 is the first exon on either strand and this
  # rule needs no strand handling of its own. min/max rather than
  # 1/nexons so that CDS-style offset ranks still work.
  gr_flagged <- gr |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      is_first = exon_rank == min(exon_rank),
      is_last = exon_rank == max(exon_rank)
    ) |>
    dplyr::ungroup()

  # separate positive and negative exons
  pos_exons <- gr_flagged |> dplyr::filter(sign(estimate) == 1)
  neg_exons <- gr_flagged |> dplyr::filter(sign(estimate) == -1)

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
      # which boundary carries the donor depends on the strand: on "+" the
      # exon end is the donor (5' splice site) and the start the acceptor,
      # on "-" transcription runs the other way and the two swap. candidates
      # are found with directed overlaps, so the candidate's strand is also
      # the partner's. "*" falls through as "+".
      on_minus = as.character(cand_tbl$strand[cand_idx]) == "-",
      is_first = cand_tbl$is_first[cand_idx],
      is_last = cand_tbl$is_last[cand_idx],
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
    dplyr::mutate(
      # this call only cares about one boundary: the donor for a5ss, the
      # acceptor for a3ss. the donor is the exon end on "+" and the exon
      # start on "-"; the acceptor is the other way round. nothing here
      # looks at the *other* boundary, so an exon that moved both is
      # reported by find_a5ss() and find_a3ss() alike.
      site_is_start = xor(!by_start, on_minus),
      site_matches = dplyr::if_else(site_is_start, match_start, match_end)
    ) |>
    dplyr::filter(
      # the site moved
      !site_matches,
      # and it exists: the last exon has no donor, the first no acceptor
      if (by_start) !is_last else !is_first
    )

  if (nrow(match_tbl) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build result: one row per (candidate, moved boundary)
  hits_tbl <- cand_tbl[match_tbl$cand_idx, ] |>
    dplyr::mutate(event_type = event_name, event_tx_id = match_tbl$tx_id_neg)
  # convert back to GRanges for return
  tbl_to_granges(hits_tbl, keep_cols, gr)
}


#' @rdname find_events
#' @return `find_a5ss()`: alternative 5' splice sites
#' @examples
#'
#' # detect alternative 5' splice sites
#' gr_a5 <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_a5ss(n_events = 1)
#' 
#' find_a5ss(gr_a5)
#'
#' @export
find_a5ss <- function(gr) {
  find_alt_ss(gr, by_start = TRUE)
}

#' @rdname find_events
#' @export
find_alternative_5_prime_splice_sites <- find_a5ss

#' @rdname find_events
#' @return `find_a3ss()`: alternative 3' splice sites
#' @examples
#'
#' # detect alternative 3' splice sites
#' gr_a3 <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_a3ss(n_events = 1)
#' find_a3ss(gr_a3)
#'
#' @export
find_a3ss <- function(gr) {
  find_alt_ss(gr, by_start = FALSE)
}

#' @rdname find_events
#' @export
find_alternative_3_prime_splice_sites <- find_a3ss

#' Find alternative transcription start / end sites
#'
#' Unlike the splice-site finders this does not go through overlaps at
#' all: it takes the first (or last) exon of every transcript and simply
#' compares where transcription begins (or ends).
#'
#' @param gr A preprocessed GRanges object
#' @param at_start If TRUE, compares the first exon of each transcript
#'   (aTSS); if FALSE, the last exon (aTES).
#' @return A GRanges of detected events
#' @noRd
find_alt_ts <- function(gr, at_start = TRUE) {
  # if preprocessing didn't happen
  check_preprocessed(gr)
  keep_cols <- keep_cols(gr)

  event_name <- if (at_start) "aTSS" else "aTES"

  # flag each transcript's first and last exon, exactly as find_alt_ss()
  # does -- but there the flags only *veto* a label, while here they
  # select which exons are considered at all. min/max rather than
  # 1/nexons so that CDS-style offset ranks still work.
  gr_flagged <- gr |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      is_first = exon_rank == min(exon_rank),
      is_last = exon_rank == max(exon_rank)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(if (at_start) is_first else is_last)

  # separate positive and negative exons
  pos_exons <- gr_flagged |> dplyr::filter(sign(estimate) == 1)
  neg_exons <- gr_flagged |> dplyr::filter(sign(estimate) == -1)

  # candidates: every terminal pos exon. unlike find_alt_ss() there is no
  # overlap requirement and no identity test -- an alternative start / end site
  # can sit anywhere in the gene, so a first exon that does not touch its
  # partner at all still counts.
  if (length(pos_exons) == 0L || length(neg_exons) == 0L) {
    return(GenomicRanges::GRanges())
  }

  cand_tbl <- tibble::as_tibble(pos_exons)
  neg_tbl <- tibble::as_tibble(neg_exons)

  # index pairs: every candidate against every negative terminal exon of
  # the same gene. this is what findOverlaps() supplies in find_alt_ss();
  # with no overlap requirement the pairing is a join on gene_id, which
  # also applies the same-gene restriction that find_alt_ss() has to
  # filter for afterwards.
  hits <- dplyr::inner_join(
    cand_tbl |>
      dplyr::mutate(cand_idx = dplyr::row_number()) |>
      dplyr::select(cand_idx, gene_id),
    neg_tbl |>
      dplyr::mutate(neg_idx = dplyr::row_number()) |>
      dplyr::select(neg_idx, gene_id),
    by = "gene_id",
    relationship = "many-to-many"
  )

  if (nrow(hits) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build match tibble with boundary checks
  match_tbl <- tibble::tibble(
    cand_idx = hits$cand_idx,
    neg_idx = hits$neg_idx
  ) |>
    dplyr::mutate(
      match_start = cand_tbl$start[cand_idx] == neg_tbl$start[neg_idx],
      match_end = cand_tbl$end[cand_idx] == neg_tbl$end[neg_idx],
      # which boundary carries the transcription start depends on the
      # strand: on "+" transcription begins at the exon start and ends at
      # the exon end, on "-" the two swap. "*" falls through as "+".
      on_minus = as.character(cand_tbl$strand[cand_idx]) == "-",
      tx_id_neg = neg_tbl$tx_id[neg_idx]
    )

  match_tbl <- match_tbl |>
    dplyr::mutate(
      # this call only cares about one boundary: the  start
      # for aTSS, the end for aTES. the same line in
      # find_alt_ss() reads xor(!by_start, on_minus), because there the
      # a5ss donor is the exon *end* on "+" rather than the start.
      site_is_start = xor(at_start, on_minus),
      site_matches = dplyr::if_else(site_is_start, match_start, match_end)
    ) |>
    dplyr::filter(
      # the site moved
      !site_matches
    )

  if (nrow(match_tbl) == 0L) {
    return(GenomicRanges::GRanges())
  }

  # build result: one row per (candidate, negative terminal exon)
  hits_tbl <- cand_tbl[match_tbl$cand_idx, ] |>
    dplyr::mutate(event_type = event_name, event_tx_id = match_tbl$tx_id_neg)
  # convert back to GRanges for return
  tbl_to_granges(hits_tbl, keep_cols, gr)
}

#' @rdname find_events
#' @return `find_aTSS()`: alternative transcription start sites — first
#'   exons that begin transcription at a different coordinate
#' @examples
#'
#' gr_tss <- data.frame(
#'   seqnames = "chr1",
#'   start = c(1, 21, 41, 5, 21, 41),
#'   end = c(10, 30, 50, 10, 30, 50),
#'   strand = "+",
#'   gene_id = "g1",
#'   tx_id = rep(c("down", "up"), each = 3),
#'   exon_rank = rep(seq_len(3), 2),
#'   estimate = rep(c(-1, 1), each = 3)
#' ) |>
#'   plyranges::as_granges() |>
#'   preprocess(coef_col = "estimate")
#'
#' find_aTSS(gr_tss)
#'
#' @export
find_aTSS <- function(gr) {
  find_alt_ts(gr, at_start = TRUE)
}

#' @rdname find_events
#' @export
find_alternative_start_sites <- find_aTSS

#' @rdname find_events
#' @return `find_aTES()`: alternative transcription end sites — last
#'   exons that end transcription at a different coordinate
#' @export
find_aTES <- function(gr) {
  find_alt_ts(gr, at_start = FALSE)
}

#' @rdname find_events
#' @export
find_alternative_end_sites <- find_aTES

#' @rdname find_events
#' @param verbose If TRUE, prints progress messages. Default TRUE.
#' @return `find_all_events()`: all detected events
#' 
#' @examples
#'
#' # detect all event types at once
#' gr_all <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' ) |>
#'   preprocess(coef_col = "estimate") |>
#'   generate_se(n_events = 1)
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

  msg("Calculating alternative transcription start site events...")
  results$aTSS <- find_aTSS(gr)

  msg("Calculating alternative transcription end site events...")
  results$aTES <- find_aTES(gr)

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
