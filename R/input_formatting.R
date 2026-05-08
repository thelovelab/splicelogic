#' Check that input is a valid GRanges object with required
#' metadata columns provided by the user
#' exon_rank", "gene_id", "tx_id", coef_col
#' @param gr A GRanges object
#' @param coef_col Name of the coefficient metadata column (string)
#' @return TRUE if input is valid, otherwise throws an error
#' @importFrom methods is
#' @noRd
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


# for preprocess
utils::globalVariables(c("nexons"))

#' Preprocess input GRanges object for splicing event calculation
#'
#' This function checks that the input is a valid GRanges
#' object with required metadata columns, then adds a unique
#' key, the number of exons per transcript, and an 'internal'
#' flag for each exon.
#'
#' @param gr A GRanges object with metadata columns:
#' 'exon_rank', 'gene_id', 'tx_id', 'coef'.
#' @param coef_col The name of the metadata column indicating
#' upregulated (+1) and downregulated (-1) exons.
#' @param method_string The Differential Transcript Usage (DTU)
#' method used to obtain the coef_col, for annotation purposes
#' (optional).
#' @param additional_columns A character vector of metadata column
#' names to record for downstream use. Stored in
#' \code{metadata(result)$additional_columns} (optional).
#'
#' @return A GRanges object with added 'key', 'nexons',
#' and 'internal' columns.
#' @export
#' @examples
#'
#' # create mock data and run preprocessing
#' gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4) |>
#'   preprocess(coef_col = "estimate", method_string = "mock_method")
#'
preprocess <- function(gr, coef_col, method_string = NULL,
                       additional_columns = NULL) {
  check_input(gr, coef_col) # check metadata columns are present

  gr_seqinfo <- GenomicRanges::seqinfo(gr)
  # use tibble for faster group_by/mutate to add key, nexons, internal
  tbl <- gr |> tibble::as_tibble()  |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      key = paste0(tx_id, "-", exon_rank),
      nexons = length(exon_rank),
      internal = exon_rank > 1 & exon_rank < nexons,
      # always use "estimate" as the column name
      # for coef values in downstream functions
      estimate = !!rlang::sym(coef_col)
    ) |>
    dplyr::ungroup()

  keep_cols <- setdiff(names(tbl), c("seqnames", "start", "end",
                                      "width", "strand"))
  gr <- GenomicRanges::GRanges(
    seqnames = tbl$seqnames,
    ranges = IRanges::IRanges(start = tbl$start, end = tbl$end),
    strand = tbl$strand,
    tbl[, keep_cols]
  )
  GenomicRanges::seqinfo(gr) <- gr_seqinfo

  S4Vectors::metadata(gr)$splicelogic_preprocessed <- TRUE
  if (!is.null(method_string)) {
    S4Vectors::metadata(gr)$method_string <- method_string
  }
  if (!is.null(additional_columns)) {
    # check that additional_columns is a character vector
    if (!is.character(additional_columns)) {
      stop("'additional_columns' must be a character vector of column names.")
    }
    # check that additional_columns are present in the GRanges metadata
    missing_cols <- setdiff(additional_columns, names(GenomicRanges::mcols(gr)))
    if (length(missing_cols) > 0) {
      stop(paste(
        "The following additional columns are missing in GRanges metadata:",
        paste(missing_cols, collapse = ", ")
      ))
    }
    S4Vectors::metadata(gr)$additional_columns <- additional_columns
  }
  return(gr)
}


#' Prepare exons from two transcript partitions
#'
#' Combines two transcript partitions (up- and down-regulated) and assigns
#' an \code{estimate} coefficient: \code{+1} to \code{up} and \code{-1} to
#' \code{down}. Accepts either GRanges objects or character vectors of
#' transcript IDs (in which case \code{txdb} is required to look up exon
#' coordinates). The result is ready to pass to \code{\link{preprocess}}
#' with \code{coef_col = "estimate"}.
#'
#' When \code{up} and \code{down} are GRanges, both must have
#' \code{exon_rank}, \code{gene_id}, and \code{tx_id} metadata columns.
#' Extra columns are kept; if one object lacks a column present in the
#' other, those entries receive \code{NA}.
#'
#' @param up A GRanges object or character vector of transcript IDs for
#'   the upregulated partition (assigned \code{estimate = +1}).
#' @param down A GRanges object or character vector of transcript IDs for
#'   the downregulated partition (assigned \code{estimate = -1}).
#' @param txdb A \code{TxDb} object (from GenomicFeatures). Required when
#'   \code{up} and \code{down} are character vectors of transcript IDs.
#' @param verbose Whether to print progress messages. Default \code{TRUE}.
#' @return A combined GRanges object with an \code{estimate} column
#'   (\code{+1} for \code{up}, \code{-1} for \code{down}),
#'   ready for \code{\link{preprocess}}.
#' @export
#' @examples
#'
#' # GRanges input
#' gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
#' gr <- generate_se(gr, n_events = 1)
#' gr_down <- gr[gr$estimate < 0]
#' gr_up <- gr[gr$estimate > 0]
#' prepare_exons_by_partition(gr_up, gr_down) |>
#'   preprocess(coef_col = "estimate")
#'
prepare_exons_by_partition <- function(up, down, txdb = NULL, verbose = TRUE) {
  both_gr   <- methods::is(up, "GRanges") && methods::is(down, "GRanges")
  both_char <- is.character(up) && is.character(down)
  if (!both_gr && !both_char) {
    stop(
      "'up' and 'down' must both be GRanges objects or ",
      "both be character vectors of tx IDs."
    )
  }
  # input are GRanges
  if (methods::is(up, "GRanges")) {
    # create estimate column with +1 for up and -1 for down
    up$estimate <- 1L
    check_input(up, "estimate")
    down$estimate <- -1L
    check_input(down, "estimate")
    plyranges::bind_ranges(up, down)
    } 
    # input are character vectors of tx IDs
    else if (is.character(up)) {
      if (is.null(txdb)) {
        stop("'txdb' is required when 'up' and 'down' are tx ID vectors.")
      }
      check_bioc_packages()
      # extract gene IDs for the provided transcript IDs
      tx_gene <- AnnotationDbi::select(
        txdb,
        keys    = c(up, down),
        columns = "GENEID",
        keytype = "TXNAME"
      ) |> tibble::as_tibble()
      # create a combined table of transcript IDs, gene IDs, and estimate values
      # for up (+1) and down (-1) partitions, and join with tx_gene to get gene IDs
      dtu_table <- tibble::tibble(
        tx_id    = c(up, down),
        estimate = c(rep(1L, length(up)), rep(-1L, length(down)))
      ) |>
        dplyr::left_join(
          dplyr::rename(tx_gene, tx_id = TXNAME, gene_id = GENEID),
          by = "tx_id"
        )
      # call prepare_exons to extract exon ranges given the txdb and dtu_table
      prepare_exons(
        txdb, dtu_table,
        coef_col = "estimate", verbose = verbose
      )
    } else {
      stop(
        "'up' and 'down' must be GRanges objects or ",
        "character vectors of tx IDs."
      )
      }
}


check_preprocessed <- function(gr) {
  if (!isTRUE(S4Vectors::metadata(gr)$splicelogic_preprocessed)) {
    stop(
      "Input has not been preprocessed with preprocess().\n",
      "  Please run preprocess() on your GRanges object before\n",
      "  calculating splicing events."
    )
  }
}

#' @noRd
check_bioc_packages <- function() {
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
}

#' @noRd
check_dtu_cols <- function(dtu_table, tx_id_col, gene_id_col, coef_col) {
  required <- c(tx_id_col, gene_id_col, coef_col)
  missing_cols <- setdiff(required, colnames(dtu_table))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in dtu_table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

#' @noRd
extract_named_ebt <- function(txdb) {
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx")
  tx_map <- AnnotationDbi::select(
    txdb,
    keys = AnnotationDbi::keys(txdb, "TXID"),
    columns = "TXNAME",
    keytype = "TXID"
  ) |>
    tibble::as_tibble()
  idx <- match(names(ebt), tx_map$TXID)
  names(ebt) <- tx_map$TXNAME[idx]
  ebt
}

#' @noRd
filter_ebt <- function(ebt, dtu_table, tx_id_col, msg) {
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
  ebt[keep]
}

#' @noRd
flatten_and_merge <- function(ebt, dtu_table, tx_id_col, gene_id_col) {
  exons <- unlist(ebt)
  exons$tx_id <- names(exons)
  names(exons) <- paste0(exons$tx_id, "-exon", exons$exon_rank)
  txp_idx <- match(exons$tx_id, dtu_table[[tx_id_col]])
  keep_idx <- !is.na(txp_idx)
  exons <- exons[keep_idx]
  txp_idx <- txp_idx[keep_idx]
  add_cols <- dtu_table[txp_idx, ] |>
    dplyr::select(-dplyr::any_of(tx_id_col))
  merged_DF <- cbind(
    GenomicRanges::mcols(exons),
    S4Vectors::DataFrame(add_cols)
  )
  GenomicRanges::mcols(exons) <- merged_DF
  if (gene_id_col != "gene_id") {
    col_names <- names(GenomicRanges::mcols(exons))
    col_names[col_names == gene_id_col] <- "gene_id"
    names(GenomicRanges::mcols(exons)) <- col_names
  }
  exons
}

#' Prepare exon ranges from a TxDb and DTU results table
#'
#' Extracts exon ranges from a TxDb object, merges them with
#' differential transcript usage (DTU) results, and returns a flat
#' GRanges ready for \code{\link{preprocess}}.
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
#' @examples
#'
#' library(AnnotationHub)
#' library(AnnotationDbi)
#' library(GenomicFeatures)
#' library(tibble)
#'
#' ah <- AnnotationHub()
#' txdb <- ah[["AH84134"]] # fly TxDb (Drosophila melanogaster)
#'
#' # build a simulated DTU table from the TxDb transcripts
#' txps <- txdb |>
#'   AnnotationDbi::select(
#'     keys(txdb, "TXID"), c("TXNAME", "GENEID"), "TXID"
#'   ) |>
#'   tibble::as_tibble() |>
#'   dplyr::select(tx_id = TXNAME, gene_id = GENEID)|>
#'   dplyr::filter(!is.na(gene_id))
#'
#' sim_dtu_table <- txps |>
#'   dplyr::mutate(
#'     padj = runif(dplyr::n()),
#'     effect_est = rnorm(dplyr::n())
#'   )
#'
#' fly_exons <- prepare_exons(
#'   txdb, sim_dtu_table, coef_col = "effect_est", verbose = TRUE
#' )
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
  check_bioc_packages()
  dtu_table <- tibble::as_tibble(dtu_table)
  check_dtu_cols(dtu_table, tx_id_col, gene_id_col, coef_col)

  msg("Extracting exons from TxDb...")
  ebt <- extract_named_ebt(txdb)
  msg("Mapping transcript IDs...")
  ebt <- filter_ebt(ebt, dtu_table, tx_id_col, msg)

  msg("Merging DTU results onto exons...")
  exons <- flatten_and_merge(ebt, dtu_table, tx_id_col, gene_id_col)

  msg(
    "Done. Returned ",
    length(exons),
    " exon ranges from ",
    length(unique(exons$tx_id)),
    " unique transcripts."
  )
  exons
}
