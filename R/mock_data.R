
#' Create a sample GRanges with two transcripts per gene and candidate logic
#' 
#' @param ngenes Number of genes to include (default is 1)
#' @param ntranscripts Number of transcripts per gene (default is 2)
#' 
#' @return A GRanges object with two transcripts per gene and candidate logic
#' @import GenomicRanges
#' @importFrom plyranges as_granges bind_ranges
#' @importFrom dplyr mutate case_when
#' @return A GRanges object with two transcripts per gene
#' @export
get_mock_data <- function() {
  
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41, 51, 61),
    width = 5,
    strand = "+",
    exon_rank = 1:7,
    gene_id = rep(1, 7),
    coefs = rep(runif(1, min = -1, max = 0), 7)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = 1:6,
    gene_id = rep(1, 6),
    coefs = rep(runif(1, min = 0, max = 1), 6)
  )

    df3 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = 1:6,
    gene_id = rep(1, 6),
    coefs = rep(runif(1, min = 0, max = 1), 6)
  )


  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr3 <- plyranges::as_granges(df3)
  gr <- plyranges::bind_ranges(gr1, gr2, gr3) |>
      plyranges::mutate(
      tx_id = c(rep(1, 7), rep(2, 6), rep(3, 6)))

return(gr)
}