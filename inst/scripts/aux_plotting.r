library(GenomicRanges)
library(IRanges)
library(wiggleplotr)
library(plyranges)
library(patchwork)

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

  #make a plot per gene_id
  plots <- list()
  n_genes <- length(unique(gr$gene_id))

  for (gene in unique(gr$gene_id)) {
    gr_gene <- gr[gr$gene_id == gene]
    if (length(gr_gene) == 0) next
    
    exons <- split(
      gr_gene,
      list(gr_gene$tx_id, !is.na(GenomicRanges::mcols(gr_gene)$event)),
      drop = TRUE
    )
    p <- wiggleplotr::plotTranscripts(
      exons = exons,
      rescale_introns = TRUE
    ) +
      ggplot2::ggtitle(gene)
    
    plots[[gene]] <- p
  }

  # Combine all plots into a grid (adjust ncol for layout)
  combined_plot <- wrap_plots(plots, ncol = round(n_genes/2 + 1))  
  combined_plot

}
  
 