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
calc_skipped_exons <- function(gr, coef_col, type = c("in","over","boundary")) {
  type <- match.arg(type) 
  # if preprocessing didn't happen
  if (!all(c("key","nexons","internal","event") %in% names(mcols(gr)))) {
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

  candidates <- switch(
    type,
    "in" = {
      # determine if the left and right exons of candidates are in the positive set
      candidates <- candidates |>
        plyranges::filter(
            left_exons %in% pos_exons &
            right_exons %in% pos_exons
        )
    },
    "over" = {
      candidates <- candidates |>
        plyranges::filter(
            left_exons %over% pos_exons &
            right_exons %over% pos_exons
        )
    },
    "boundary" = {
      has_end_in   <- function(x, y) end(x)   %in% end(y)
      has_start_in <- function(x, y) start(x) %in% start(y)
      candidates <- candidates |>
        plyranges::filter(
          (left_exons  |> has_end_in(pos_exons)) &
          (right_exons |> has_start_in(pos_exons))
        )
    }
  )
  gr |>
    plyranges::mutate(
      event = ifelse(key %in% candidates$key, "skipped_exon", event)
    )

}
