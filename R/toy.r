library(GenomicRanges)
library(IRanges)
library(wiggleplotr)
library(plyranges)

gr <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 100, end = 200),
  strand = "+"
)


gr2 <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 201, end = 250),
  strand = "+"
)

  
pltRanges <- function(gr) {

  if (!all(c("event") %in% names(mcols(gr)))) {
    gr$event <- NA
  }

  if (!"tx_id" %in% names(mcols(gr))) {
    mcols(gr)$tx_id <- 1L
  }

  # Split 
  exons <- split(gr, list(gr$tx_id, !is.na(GenomicRanges::mcols(gr)$event)), drop = TRUE)
   
  p <- wiggleplotr::plotTranscripts(
    exons = exons,
    rescale_introns = TRUE
  )
  return(p)
}
  
 