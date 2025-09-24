#' Calculate skipped exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'tx_id', 'exon',
#' and 'coef_col' metadata columns.
#' @param coef_col The name of the metadata column indicating upregulated (+1) and downregulated (-1) exons.
#' @param type The type of overlap to consider when identifying skipped exons.
#' @return A GRanges object with an additional 'event' metadata column indicating skipped exons.
#' @import GenomicRanges
#' @importFrom dplyr group_by mutate ungroup filter
#' @importFrom plyranges filter_by_non_overlaps_directed
#' @export
calc_skipped_exons <- function(gr, coef_col, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  # separate positive and negative exons
  var <- rlang::sym(coef_col)
  pos_exons <- gr |> plyranges::filter(sign(!!var) == 1)
  neg_exons <- gr |> plyranges::filter(sign(!!var) == -1)

  res <- candidates_by_non_overlap(neg_exons, pos_exons, gr)
  candidates <- res$candidates
  left_exons <- res$left_exons
  right_exons <- res$right_exons

  hits <- GRanges()
  if (length(candidates) == 0L) {
    return(hits) #return unchanged if no candidates
  }

  for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges

    matches <- match_left_right(
      pos_exons,
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
        plyranges::mutate(
          event    = "skipped_exon",
          tx_event = txs
        )
      hits <- c(hits, cand_rep)
      }
    }
  hits
  }


calc_mutually_exclusive <- function(gr, coef_col, type = c("in", "over", "boundary")) {
  type <- match.arg(type)
  # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  # separate positive and negative exons
  var <- rlang::sym(coef_col)
  pos_exons <- gr |> plyranges::filter(sign(!!var) == 1)
  neg_exons <- gr |> plyranges::filter(sign(!!var) == -1)

  res <- candidates_by_non_overlap(neg_exons, pos_exons, gr)
  candidates <- res$candidates
  left_exons <- res$left_exons
  right_exons <- res$right_exons


  hits <- GRanges()
  if (length(candidates) == 0L) {
    return(hits) #return unchanged if no candidates
  }

  for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges

    matches <- match_left_right(
      pos_exons,
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
      mx_pos_exon <- pos_exons |>
        plyranges::filter(tx_id == tx_event & exon_rank == exon_rank_event)
      if (length(mx_pos_exon) == 1L) {
        cand_hit <- cand |>
          plyranges::mutate(
            event    = "mutually_exclusive",
            tx_event = tx_event
          )
        pos_hit <- mx_pos_exon |>
          plyranges::mutate(
            event    = "mutually_exclusive",
            tx_event = cand_hit$tx_id
          )
        hits <- c(hits, cand_hit, pos_hit)
      }
    }
  }
  hits
}
