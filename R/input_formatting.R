#' Check that input is a valid GRanges object with required
#' metadata columns provided by the user
#' exon_rank", "gene_id", "tx_id", coef_col
#' @param gr A GRanges object
#' @param coef_col Name of the coefficient metadata column (string)
#' @return TRUE if input is valid, otherwise throws an error
#' @importFrom methods is
check_input <- function(gr, coef_col) {
  if (!is(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  #check for required metadata columns
  required_cols <- c("exon_rank", "gene_id", "tx_id", coef_col)
  missing_cols <- setdiff(required_cols, names(GenomicRanges::mcols(gr)))
  if (length(missing_cols) > 0) {
    stop(paste(
      "Missing required metadata columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }
  #check if coef_col is present and valid
  if (coef_col %in% names(GenomicRanges::mcols(gr))) {
    vals <- GenomicRanges::mcols(gr)[[coef_col]]
    if (!is.numeric(vals)) {
      stop(sprintf(
        "The '%s' metadata column must contain numeric values.",
        coef_col
      ))
    }
  }

  TRUE
}

#' Combine two GRanges objects into one for splicing event analysis
#' This function takes two GRanges objects, typically representing
#' positive and negative sets of exons, and combines them into
#' a single GRanges object.
#' @param gr1 A GRanges object (e.g., positive set)
#' @param gr2 A GRanges object (e.g., negative set)
#' @param coef_col Name of the coefficient metadata column (string)
#' @return A combined GRanges object with appropriate coef metadata
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

  gr <- plyranges::bind_ranges(gr1, gr2)
  return(gr)
}

#' Preprocess input GRanges object for splicing event calculation
#'
#' This function checks that the input is a valid GRanges
#' object with required metadata columns, then adds a unique
#' key, the number of exons per transcript, and an 'internal'
#' flag for each exon. It also initializes an 'event' column
#' for downstream splicing event annotation.
#'
#' @param gr A GRanges object with metadata columns:
#' 'exon_rank', 'gene_id', 'tx_id', 'coef'.
#' @param coef_col The name of the metadata column indicating
#' upregulated (+1) and downregulated (-1) exons.
#' @param method_string The Differential Transcript Usage (DTU)
#' method used to obtain the coef_col, for annotation purposes
#' (optional).
#'
#' @return A GRanges object with added 'key', 'nexons',
#' 'internal', and 'event' columns.
#' @export
preprocess_input <- function(gr, coef_col, method_string = NULL) {
  check_input(gr, coef_col) # check metadata columns are present

  # include key nexons and internal columns
  gr <- gr |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      key = paste0(tx_id, "-", exon_rank),
      nexons = length(exon_rank),
      internal = exon_rank > 1 & exon_rank < nexons,
      # always use "estimates" as the column name
      # for coef values in downstream functions
      estimates = !!rlang::sym(coef_col)
    ) |>
    dplyr::ungroup()

  S4Vectors::metadata(gr)$splicelogic_preprocessed <- TRUE
  if (!is.null(method_string)) {
    S4Vectors::metadata(gr)$method_string <- method_string
  }
  S4Vectors::metadata(gr)$splicelogic_preprocessed <- TRUE
  return(gr)
}


check_preprocessed <- function(gr) {
  if (!isTRUE(S4Vectors::metadata(gr)$splicelogic_preprocessed)) {
    stop(
      "Input has not been preprocessed with preprocess_input().\n",
      "  Please run preprocess_input() on your GRanges object before\n",
      "  calculating splicing events."
    )
  }
}

#' Prepare exon ranges from a TxDb and DTU results table
#'
#' Extracts exon ranges from a TxDb object, merges them with
#' differential transcript usage (DTU) results, and returns a flat
#' GRanges ready for \code{\link{preprocess_input}}.
#'
#' @param txdb A \code{TxDb} object (from GenomicFeatures).
#' @param dtu_table A data.frame or tibble with DTU results. Must
#'   contain columns for transcript ID, gene ID, and a coefficient.
#' @param coef_col Column name in \code{dtu_table} with the
#'   coefficient / effect size values.
#' @param tx_id_col Column name in \code{dtu_table} with transcript
#'   IDs matching the TxDb transcript names. Default \code{"tx_id"}.
#' @param gene_id_col Column name in \code{dtu_table} with gene IDs.
#'   Default \code{"gene_id"}.
#' @param verbose Whether to print progress messages. Default \code{TRUE}.
#' @return A GRanges object with metadata columns: \code{gene_id},
#'   \code{tx_id}, \code{exon_rank}, the coefficient column, and any
#'   additional columns from \code{dtu_table}.
#' @export
prepare_exons <- function(
  txdb,
  dtu_table,
  coef_col,
  tx_id_col = "tx_id",
  gene_id_col = "gene_id",
  verbose = TRUE
) {
  msg <- if (verbose) message else function(...) invisible(NULL)

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "Package 'GenomicFeatures' is required. ",
      "Install with: BiocManager::install('GenomicFeatures')"
    )
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Package 'AnnotationDbi' is required. ",
      "Install with: BiocManager::install('AnnotationDbi')"
    )
  }

  dtu_table <- tibble::as_tibble(dtu_table)
  required <- c(tx_id_col, gene_id_col, coef_col)
  missing_cols <- setdiff(required, colnames(dtu_table))
  if (length(missing_cols) > 0) {
    stop("Missing columns in dtu_table: ", paste(missing_cols, collapse = ", "))
  }

  # extract exons grouped by transcript
  msg("Extracting exons from TxDb...")
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx")

  # map TxDb internal TXID to TXNAME
  msg("Mapping transcript IDs...")
  tx_map <- suppressMessages(AnnotationDbi::select(
    txdb,
    keys = AnnotationDbi::keys(txdb, "TXID"),
    columns = "TXNAME",
    keytype = "TXID"
  )) |>
    tibble::as_tibble()

  # map TXID to TXNAME
  idx <- match(names(ebt), tx_map$TXID)
  names(ebt) <- tx_map$TXNAME[idx]

  # check that dtu_table tx_ids match TxDb transcript names
  keep <- names(ebt) %in% dtu_table[[tx_id_col]]
  if (!any(keep)) {
    stop(
      "No matching transcript IDs between TxDb and dtu_table$",
      tx_id_col,
      ". Check that tx_id_col contains TXNAME values ",
      "(e.g. ENST...), not internal TXID integers."
    )
  }
  n_missing <- sum(!dtu_table[[tx_id_col]] %in% names(ebt))
  if (n_missing > 0) {
    msg(
      n_missing,
      " transcript(s) in dtu_table not found in TxDb ",
      "and will be excluded."
    )
  }
  ebt <- ebt[keep]

  # flatten GRangesList to GRanges
  exons <- unlist(ebt)
  exons$tx_id <- names(exons)
  names(exons) <- exons$exon_name

  # merge dtu_table columns onto exons by tx_id
  msg("Merging DTU results onto exons...")
  txp_idx <- match(exons$tx_id, dtu_table[[tx_id_col]])
  add_cols <- dtu_table[txp_idx, ] |>
    dplyr::select(-dplyr::any_of(tx_id_col))
  merged_DF <- cbind(
    GenomicRanges::mcols(exons),
    S4Vectors::DataFrame(add_cols)
  )
  GenomicRanges::mcols(exons) <- merged_DF

  # rename to standard column names expected by preprocess_input
  if (gene_id_col != "gene_id") {
    col_names <- names(GenomicRanges::mcols(exons))
    col_names[col_names == gene_id_col] <- "gene_id"
    names(GenomicRanges::mcols(exons)) <- col_names
  }

  msg(
    "Done. Returned ",
    length(exons),
    " exon ranges from ",
    length(unique(exons$tx_id)),
    "unique transcripts."
  )
  exons
}
