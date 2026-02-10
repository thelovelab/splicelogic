#' Basic function  to visualize easily a gr object with transcript coefficients and events.
#' @param gr A GRanges object with metadata columns: 'gene_id', 'tx_id', 'coefs', and optionally 'event'.
#' @return A ggplot object visualizing the transcripts and events for each gene.
#' @export 
slPltRanges <- function(gr) {

  if (!all(c("event") %in% names(GenomicRanges::mcols(gr)))) {
    gr$event <- NA
  }

  if (!"tx_id" %in% names(GenomicRanges::mcols(gr))) {
    GenomicRanges::mcols(gr)$tx_id <- 1L
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
      dplyr::mutate(
        transcript_id = rep(names(exons), lengths(exons)),
        color_by   = ifelse(coefs > 0, "blue", "red")
      ) |> dplyr::as_tibble() |>
      dplyr::select(
        transcript_id, strand, color_by)|> unique()


    p <- slPlotTranscripts(
      exons = exons,
      transcript_annotations = transcript_annotations, 
      rescale_introns = TRUE
    ) +
      ggplot2::ggtitle(gene)
    
    plots <- append(plots, list(p))
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
