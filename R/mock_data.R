
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
    start = c(1, 11, 21, 31),
    width = 5,
    strand = "+",
    exon_rank = 1:4,
    gene_id = rep(1, 4),
    tx_id = rep(1, 4),
    coef = rep(-1, 4)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  31),
    width = 5,
    strand = "+",
    exon_rank = 1:3,
    gene_id = rep(1, 3),
    tx_id = rep(1, 3),
    coef = rep(1, 3)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
      plyranges::mutate(
      tx_id = c(rep(1, 4), rep(2, 3)))

return(gr)
}