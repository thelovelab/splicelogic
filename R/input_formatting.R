
#' Check that input is a valid GRanges object with required metadata columns
#' @param gr A GRanges object
#' @return TRUE if input is valid, otherwise throws an error
#' @importFrom GenomicRanges mcols
#' @importFrom plyranges bind_ranges
#' @export
check_input <- function(gr) {
  if (!is(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  
  required_cols <- c("exon_rank", "gene_id", "tx_id", "coef")
  missing_cols <- setdiff(required_cols, names(mcols(gr)))

  #check if coef column is present and valid
  if ("coef" %in% names(mcols(gr))) {
    if (any(!mcols(gr)$coef %in% c(1, -1))) {
      stop("The 'coef' metadata column must contain only 1 (contrast) or -1 (reference) values.")
    }
  }

  if (length(missing_cols) > 0) {
    stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
  }
  
  TRUE
}

combine_gr_input <- function(gr1, gr2) {
  if (!is(gr1, "GRanges") || !is(gr2, "GRanges")) {
    stop("Both inputs must be GRanges objects.")
  }
  #check if they have coef metadata column, add it if missing 
  if (!"coef" %in% names(mcols(gr1))) {
    mcols(gr1)$coef <- +1 #gr1 is the contrast
  }
  if (!"coef" %in% names(mcols(gr2))) {
    mcols(gr2)$coef <- -1 #gr2 is the reference
  }

  gr <- plyranges::bind_ranges(gr1,gr2)
  return(gr)
}

#' Preprocess input GRanges object for splicing event calculation
#'
#' This function checks that the input is a valid GRanges object with required metadata columns,
#' then adds a unique key, the number of exons per transcript, and an 'internal' flag for each exon.
#' It also initializes an 'event' column for downstream splicing event annotation.
#'
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id', 'tx_id', 'coef'.
#' @return A GRanges object with added 'key', 'nexons', 'internal', and 'event' columns.
#' @importFrom plyranges group_by mutate ungroup
preprocess_input <- function(gr) {
  check_input(gr) # check metadata columns are present

  # include key nexons and internal columns
  gr <- gr |>
    group_by(tx_id) |> 
    mutate(
      key     = paste0(tx_id, "-", exon_rank),
      nexons  = length(exon_rank),
      internal = exon_rank > 1 & exon_rank < nexons
    ) |> 
    ungroup()

  # initialize event column if needed 
  if (!"event" %in% colnames(mcols(gr))) {
    gr$event <- NA_character_
  }
  return(gr)
}