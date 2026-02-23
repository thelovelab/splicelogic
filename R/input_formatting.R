#' Check that input is a valid GRanges object with required metadata columns
#' exon_rank", "gene_id", "tx_id", coef_col
#' @param gr A GRanges object
#' @param coef_col Name of the coefficient metadata column (string)
#' @return TRUE if input is valid, otherwise throws an error
#' @keywords internal
check_input <- function(gr, coef_col) {
  if (!is(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  #check for required metadata columns
  required_cols <- c("exon_rank", "gene_id", "tx_id", coef_col)
  missing_cols <- setdiff(required_cols, names(GenomicRanges::mcols(gr)))
  if (length(missing_cols) > 0) {
      stop(paste("Missing required metadata columns:", paste(missing_cols, collapse = ", ")))
    }
  #check if coef_col is present and valid
  if (coef_col %in% names(GenomicRanges::mcols(gr))) {
    vals <- GenomicRanges::mcols(gr)[[coef_col]]
    if (any(vals < -1 | vals > 1)) {
      stop(sprintf("The '%s' metadata column must contain values only in the range [-1, 1].", coef_col))
    }
  }

  TRUE
}

#' Combine two GRanges objects into one for splicing event analysis
#' This function takes two GRanges objects, typically representing
#' positive and negative sets of exons, and combines them into a single GRanges object.
#' @param gr1 A GRanges object (e.g., positive set)
#' @param gr2 A GRanges object (e.g., negative set)
#' @param coef_col Name of the coefficient metadata column (string)
#' @return A combined GRanges object with appropriate coef metadata
#' @keywords internal
combine_gr_input <- function(gr1, gr2, coef_col) {
  if (!is(gr1, "GRanges") || !is(gr2, "GRanges")) {
    stop("Both inputs must be GRanges objects.")
  }
  coef <- rlang::sym(coef_col)
  #check if they have coef metadata column, add it if missing 
  if (!coef_col %in% names(GenomicRanges::mcols(gr1))) {
    GenomicRanges::mcols(gr1)[[coef_col]] <- +1 #gr1 is the contrast
  }
  if (!coef_col %in% names(GenomicRanges::mcols(gr2))) {
    GenomicRanges::mcols(gr2)[[coef_col]] <- -1 #gr2 is the reference
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
#' @param coef_col blah
#' @param method_string blah
#' 
#' @return A GRanges object with added 'key', 'nexons', 'internal', and 'event' columns.
#' @export
preprocess_input <- function(gr, coef_col, method_string=NULL) {
  check_input(gr, coef_col) # check metadata columns are present

  # include key nexons and internal columns
  gr <- gr |>
    dplyr::group_by(tx_id) |> 
    dplyr::mutate(
      key     = paste0(tx_id, "-", exon_rank),
      nexons  = length(exon_rank),
      internal = exon_rank > 1 & exon_rank < nexons
    ) |> 
    dplyr::ungroup()

  metadata(gr)$splicelogic_preprocessed <- TRUE
  if (!is.null(method_string)) {
    metadata(gr)$method_string <- method_string
  }

  return(gr)
}