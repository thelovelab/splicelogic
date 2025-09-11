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

  # if preprocessing didn't happen
  if (!all(c("key","nexons","internal","event") %in% names(mcols(gr)))) {
    gr <- gr | preprocess_input(gr)
  }
  
  pos_exons <- gr |> filter(coef == 1)
  neg_exons <- gr |> filter(coef == -1)

  candidates <- neg_exons |>
    filter_by_non_overlaps_directed(pos_exons)
  
  if (length(candidates) == 0L) {
    return(
      candidates |>
        plyranges::mutate(skipped_exon = NA)
    )
  }
  
  candidates <- candidates |>
    plyranges::filter(internal)

  # keys of exons to the left and right of candidates
  left_keys  <- paste0(candidates$txp, "-", 
                        candidates$rank - 1)
  right_keys <- paste0(candidates$txp, "-", 
                        candidates$rank + 1)

  # get the actual exons for the candidates
  left_exons  <- gr |> filter(key %in% left_keys)
  right_exons <- gr |> filter(key %in% right_keys)

  # determine if the left and right exons of candidates
  # were somewhere in the positive exon set
  candidates <- candidates |>
    filter(
        left_exons %in% pos_exons &
        right_exons %in% pos_exons
    )

  gr |>
    mutate(
      event = ifelse(key %in% candidates$key, "skipped_exon", event)
    )

}
