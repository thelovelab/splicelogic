#' Calculate skipped exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns.
#' @param coef_col The name of the metadata column indicating upregulated (+1) and downregulated (-1) exons.
#' @param type The type of overlap to consider when identifying skipped exons.
#' @return A GRanges object with an additional 'event' metadata column indicating skipped exons.
#' @export
calc_skipped_exons <- function(gr, coef_col, type = c("over","in", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(GenomicRanges::mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  # separate positive and negative exons
  coef <- rlang::sym(coef_col)
  pos_exons <- gr |> dplyr::filter(sign(!!coef) == 1)
  neg_exons <- gr |> dplyr::filter(sign(!!coef) == -1)

  # filter_results <- candidates_by_non_overlap_directed(neg_exons, pos_exons, gr)
  filter_results <- candidates_by_presence(gr, coef_col) # TO DO : check by overlap not in 
  # filter_results <- candidates_by_presence_v2(neg_exons, pos_exons)

  candidates <- filter_results$candidates
  left_exons <- filter_results$left_exons
  right_exons <- filter_results$right_exons

  hits <- GRanges()
  if (length(candidates) == 0L) {
    return(hits) #return unchanged if no candidates
  }

  for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges
    # restric to the same gene
    cand_pos_exons <- pos_exons |> 
                  dplyr::filter(gene_id == cand$gene_id)

    matches <- match_left_right(
      cand_pos_exons,
      left_exon  = left_exons[i],
      right_exon = right_exons[i],
      type = type
    )
    left_tbl <- matches$left_tbl
    right_tbl <- matches$right_tbl

    # join left and right by tx_id, filter for adjacent exons (l-r==1)
    pairs <- dplyr::inner_join(left_tbl, right_tbl, by = "tx_id") |>
        dplyr::filter(abs(l - r) == 1)

    if (nrow(pairs) > 0) {
      txs <- unique(pairs$tx_id)
      cand_rep <- rep(cand, length(txs)) |>
        dplyr::mutate(
          event    = "skipped_exon",
          tx_event = txs
        )
      hits <- c(hits, cand_rep)
      }
    }
  hits
  }

#' Calculate mutually exclusive exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns.
#' @param coef_col The name of the metadata column indicating upregulated (+1) and downregulated (-1) exons.
#' @param type The type of overlap to consider when identifying mutually exclusive exons.
#' @return A GRanges object with an additional 'event' metadata column indicating mutually exclusive exons.
#' @export
calc_mutually_exclusive <- function(gr, coef_col, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(GenomicRanges::mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  # separate positive and negative exons
  coef <- rlang::sym(coef_col)
  pos_exons <- gr |> dplyr::filter(sign(!!coef) == 1)
  neg_exons <- gr |> dplyr::filter(sign(!!coef) == -1)
  
  # returns a list of GRanges of same length:
  # ‘candidates’ - neg exons that do not overlap any pos exons, and are internal
  # ‘left_exons’ - left exons of candidates
  # ‘right_exons’ - right exons of candidates
  filter_results <- candidates_by_non_overlap_directed(neg_exons, pos_exons, gr)
  candidates <- filter_results$candidates
  left_exons <- filter_results$left_exons
  right_exons <- filter_results$right_exons

  hits <- GRanges()
  if (length(candidates) == 0L) {
    return(hits) #return unchanged if no candidates
  }

  for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges
    
    # restric to the same gene
    cand_pos_exons <- pos_exons |> 
                  dplyr::filter(gene_id == cand$gene_id)
    matches <- match_left_right(
      cand_pos_exons,
      left_exon  = left_exons[i],
      right_exon = right_exons[i],
      type = type
    )
    left_tbl <- matches$left_tbl
    right_tbl <- matches$right_tbl

    # join left and right by tx_id, filter for mx exons (l-r==2)
    pairs <- dplyr::inner_join(left_tbl, right_tbl, by = "tx_id") |>
        dplyr::filter(abs(l - r) == 2)

    for (i in seq_along(pairs)) {
      tx_event <- pairs$tx_id[i]
      exon_rank_event <- pairs$r[i] - 1  # middle exon rank
      mx_pos_exon <- cand_pos_exons |>
        dplyr::filter(tx_id == tx_event & exon_rank == exon_rank_event)
      if (length(mx_pos_exon) == 1L) {
        cand_hit <- cand |>
          dplyr::mutate(
            event    = "mutually_exclusive",
            tx_event = tx_event
          )
        pos_hit <- mx_pos_exon |>
          dplyr::mutate(
            event    = "mutually_exclusive",
            tx_event = cand_hit$tx_id
          )
        hits <- c(hits, cand_hit, pos_hit)
      }
    }
  }
  hits
}

#' Function to calculate retained introns given a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', and 'coef'.
#' @return A GRanges object with an additional 'event' metadata column indicating retained introns.
#' @export      
calc_retained_introns <- function(gr, coef_col){
  # TO DO : fix coef_col argument and add to preprocess_input and then rename coef_col to coefs
  if (!all(c("key", "nexons", "internal", "event") %in% names(GenomicRanges::mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
    }
  # separate positive and negative exons
  coef <- rlang::sym(coef_col)
  # find introns in the negative coef transcripts
  introns <- gr |> dplyr::filter(!!coef < 0) |> find_introns()

  # filter by overlap on the positive coef transcripts
  candidates <- gr |> dplyr::filter(!!coef > 0) |> 
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
calc_a3ss_a5ss <- function(gr, coef_col ){
  # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(GenomicRanges::mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  
  # separate positive and negative exons
  coef <- rlang::sym(coef_col)
  pos_exons <- gr |> dplyr::filter(sign(!!coef) == 1)
  neg_exons <- gr |> dplyr::filter(sign(!!coef) == -1)

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