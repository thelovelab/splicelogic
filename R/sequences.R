#' Extract sequences for GRanges
#'
#' Given a GRanges, look up the DNA or RNA sequence 
#' of each range in a reference genome. This is a simple
#' helper function that rearranges the `getSeq` arguments
#' and combines it with other helper functions to refer
#' to genomes by name.
#'
#' @param gr A GRanges object.
#' @param genome Either a character string naming an installed BSgenome
#'   package (e.g. \code{"hg38"}, or \code{"BSgenome.Hsapiens.UCSC.hg38"}) 
#'   or a BSgenome object. See \code{?getBSgenome} for details.
#' @param as_rna If \code{TRUE}, return an \code{RNAStringSet}
#'   (\code{T -> U}). Default \code{FALSE} returns a \code{DNAStringSet}.
#' 
#' @return A \code{DNAStringSet} (or \code{RNAStringSet} if
#'   \code{as_rna = TRUE}), one entry per range in \code{gr}, in the
#'   same order. Assign it onto \code{gr} (e.g. \code{gr$seq <- ...})
#'   to keep it as a metadata column.
#' 
#' @examples
#' 
#'   gr <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(start = c(1e6, 1.1e6), width = 50)
#'   )
#' 
#'   suppressPackageStartupMessages(
#'     library(Biostrings)
#'   )
#' 
#'   get_seq(gr, "hg38")
#'   get_seq(gr, "hg38", as_rna = TRUE)
#'
#'   set.seed(123)
#'   gr <- create_mock_data(
#'     n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6
#'   ) |>
#'     generate_se(n_events = 1) |>
#'     GenomicRanges::shift(50e6) # move mock data
#' 
#'   skipped <- find_se(gr)
#'   skipped %>% 
#'     plyranges::flank_upstream(100) %>%
#'     dplyr::mutate(seq = get_seq(., "hg38"))
#' 
#' @export
get_seq <- function(gr, genome, as_rna = FALSE) {
  if (!methods::is(gr, "GRanges")) {
    stop("'gr' must be a GRanges object.")
  }
  if (!is.character(genome) && !methods::is(genome, "BSgenome")) {
    stop(
      "'genome' must be a BSgenome object or a character string ",
      "naming an installed BSgenome package."
    )
  }
  if (is.character(genome) &&
        (length(genome) != 1L || is.na(genome))) {
    stop("'genome' must be a single, non-NA character string.")
  }
  if (!is.logical(as_rna) || length(as_rna) != 1L || is.na(as_rna)) {
    stop("'as_rna' must be a single logical value (TRUE or FALSE).")
  }

  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop(
      "Package 'BSgenome' is required. ",
      "Install with: BiocManager::install('BSgenome')"
    )
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
      "Package 'Biostrings' is required. ",
      "Install with: BiocManager::install('Biostrings')"
    )
  }
  
  if (is.character(genome)) {
    genome <- BSgenome::getBSgenome(genome)
  }

  seqs <- Biostrings::getSeq(genome, gr)
  
  if (as_rna) {
    seqs <- Biostrings::RNAStringSet(seqs)
  }

  seqs

}
