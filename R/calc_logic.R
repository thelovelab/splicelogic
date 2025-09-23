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

  candidates <- neg_exons |>
    plyranges::filter_by_non_overlaps_directed(pos_exons)
  
  if (length(candidates) == 0L) {
    return(gr) #return unchanged if no candidates
  }

  candidates <- candidates |>
    plyranges::filter(internal)

  # keys of exons to the left and right of candidates
  left_keys  <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank - 1)
  right_keys <- paste0(candidates$tx_id, "-", 
                        candidates$exon_rank + 1)

  # get the actual exons for the candidates
  left_exons <- gr |> plyranges::slice(match(left_keys, key))
  right_exons <- gr |> plyranges::slice(match(right_keys, key))

  hits <- GRanges()
  for (i in seq_along(candidates)) {
    cand <- candidates[i]  # a length-1 GRanges

    # split by tx_id (as tibbles) and check each transcript separately
    temp_pos <- get_matcher(
      pos_exons,
      left_exon  = left_exons[i],
      right_exon = right_exons[i],
      type = type
    )
    # dplyr::glimpse(temp_pos)  

    # left/right exon ranks per tx
    left_tbl <- temp_pos |>
      dplyr::filter(match_left == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(l = exon_rank)

    right_tbl <- temp_pos |>
      dplyr::filter(match_right == TRUE) |>
      dplyr::distinct(tx_id, exon_rank) |>
      dplyr::rename(r = exon_rank)

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