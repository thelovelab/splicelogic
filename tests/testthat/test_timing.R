test_that("timing as the number of transcripts grows large", {

  ngenes <- 10
  # create mock data with 10 genes, 2 transcripts per gene, and 5 exons per transcript 
  gr <- create_mock_data(n_genes = ngenes, n_tx_per_gene=2, n_exons_per_tx=5)
  table(gr$gene_id)

  # generate some skipped exon events in the mock data
  gr <- generate_skipped_exons(gr, n_se = 2)
  # view the plot to check that the events look correct
  # slPltRanges(gr)  
  
  gr <- preprocess_input(gr, coef_col = "coefs")
  res <- calc_skipped_exons(gr)
  
})
