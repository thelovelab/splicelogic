#' Calculate skipped exons from a GRanges object
#' @param gr A GRanges object with exon annotations, including 'txp', 'rank',
#' and 'coef' metadata columns.
#' The 'coef' column should contain 1 for upregulated exons and -1 for downregulated exons.
#' @return A GRanges object with an additional 'event' metadata column indicating skipped exons.
#' @import GenomicRanges
#' @importFrom dplyr group_by mutate ungroup filter
#' @importFrom plyranges filter_by_non_overlaps_directed
#' @export          
calc_skipped_exons <- function(gr) {

  pos_exons <- gr |> filter(coef == 1)
  neg_exons <- gr |> filter(coef == -1)

  tmp <- neg_exons |> filter_by_non_overlaps_directed(pos_exons)
  candidates <- if (length(tmp) == 0L) {
    tmp # return the empty GRanges unchanged to avoid error in next step
  } else {
    tmp |>
      plyranges::mutate(SE = rep_len(TRUE, length(tmp))) |>
      plyranges::filter(internal)
  }

  left_keys  <- paste0(candidates$txp, "-", 
                        candidates$rank - 1)
  right_keys <- paste0(candidates$txp, "-", 
                        candidates$rank + 1)

  left_exons  <- gr |> filter(key %in% left_keys)
  right_exons <- gr |> filter(key %in% right_keys)

  candidates <- candidates |>
    mutate(
      left_and_right =
        left_exons %in% pos_exons &
        right_exons %in% pos_exons
    )

    # # create skipped exons and bind to original gr OLD WAY
    # skipped_exons <- candidates |> 
    # plyranges::filter(left_and_right == TRUE) |> 
    # plyranges::mutate(event = ifelse(left_and_right %in% TRUE,
    # "skipped_exon", NA_character_))
    # gr <- plyranges::bind_ranges(gr, skipped_exons)

  # find keys of candidates that are true skipped exons
  skipped_exons_keys <- candidates$key[candidates$left_and_right %in% TRUE]
  # tag them in the original gr
  gr$event[gr$key %in% skipped_exons_keys] <- "skipped_exon"

  return(gr)
}