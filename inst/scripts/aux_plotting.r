library(GenomicRanges)
library(IRanges)
library(wiggleplotr)
library(plyranges)
library(patchwork)


# gr <- GRanges(
#   seqnames = "chr1",
#   ranges = IRanges(start = c(100,300), end = c(200,400)),
#   strand = "+"
# )

# gr2 <- GRanges(
#   seqnames = "chr1",
#   ranges = IRanges(start = c(100,350), end = c(200,400)),
#   strand = "+"
# )

# gr |> 
#     plyranges::filter(!(gr %in% gr2)) |>
#     plyranges::filter_by_overlaps_directed(gr2)

   
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
      list(gr_gene$tx_id),
      drop = FALSE
    )

    transcript_annotations <- exons |>
      unlist(use.names = FALSE) |>
      plyranges::mutate(
        transcript_id = rep(names(exons), lengths(exons)),
        color_by   = ifelse(coefs > 0, "blue", "red")
      ) |> as_tibble() |>
      dplyr::select(
        transcript_id, strand, color_by)|> unique()


    p <- myPlotTranscripts(
      exons = exons,
      transcript_annotations = transcript_annotations, 
      rescale_introns = TRUE
    ) +
      ggplot2::ggtitle(gene)
    
    plots[[gene]] <- p
  }

  # Combine all plots into a grid (adjust ncol for layout)
  if (length(plots) == 0L) {
  return(invisible(NULL))
} else if (length(plots) == 1L) {
  print(plots[[1]])
  return(invisible(plots[[1]]))
} else {
  combined_plot <- patchwork::wrap_plots(plots, ncol = max(1, round(n_genes/2 + 1)))
  print(combined_plot)
  return(invisible(combined_plot))
}

}
  
 

pairwise_overlaps
function(gr) {
    idx <- expand.grid(i = seq_along(gr), j = seq_along(gr))
    overlapsAny(gr[idx$i], gr[idx$j])
}