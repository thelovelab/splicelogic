test_that("timing as the number of transcripts grows large", {

  ngenes <- 10
  gr <- create_mock_data(n_genes = ngenes, n_tx_per_gene=2, n_exons_per_tx=5)
  table(gr$gene_id)
  
  gr <- preprocess_input(gr, "coefs")
  
  res <- calc_skipped_exons(gr, "coefs")
  
})
