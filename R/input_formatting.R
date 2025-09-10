
#' Check that input is a valid GRanges object with required metadata columns
#' @param gr A GRanges object
#' @return TRUE if input is valid, otherwise throws an error
#' @importFrom GenomicRanges mcols
#' @export
check_input <- function(gr) {
  if (!is(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  
  required_cols <- c("rank", "gene", "txp", "coef")
  missing_cols <- setdiff(required_cols, names(mcols(gr)))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
  }
  
  TRUE
}

#' Preprocess input GRanges object for splicing event calculation
#'
#' This function checks that the input is a valid GRanges object with required metadata columns,
#' then adds a unique key, the number of exons per transcript, and an 'internal' flag for each exon.
#' It also initializes an 'event' column for downstream splicing event annotation.
#'
#' @param gr A GRanges object with metadata columns: 'rank', 'gene', 'txp', 'coef'.
#' @return A GRanges object with added 'key', 'nexons', 'internal', and 'event' columns.
#' @importFrom plyranges group_by mutate ungroup
preprocess_input <- function(gr) {
  check_input(gr) # check metadata columns are present

  # include key nexons and internal columns
  gr <- gr |>
    group_by(txp) |> 
    mutate(
      key     = paste0(txp, "-", rank),
      nexons  = length(rank),
      internal = rank > 1 & rank < nexons
    ) |> 
    ungroup()

  # initialize event column
  gr$event <- NA_character_
  return(gr)
}