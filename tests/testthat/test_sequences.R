# Test get_seq on the output of find_se after a plyranges flank
test_that("get_seq adds a seq column to flank_upstream(find_se()) output", {
  skip_if_not_installed("BSgenome")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("plyranges")

  gr <- create_mock_data(
    n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6
  ) |>
    generate_se(n_events = 1)
  skipped <- find_se(gr)

  # mock coords start near 1, so flank_upstream(100) would run off the
  # chromosome; shift into a real region of hg38 before sequence lookup.
  skipped <- GenomicRanges::shift(skipped, 1e6)

  result <- skipped |>
    plyranges::flank_upstream(100) 
  result$seq <- get_seq(result, genome = "hg38", as_rna = TRUE)
  
  # Alternatively, could do in one step if %>% imported
  # result <- skipped |>
  #   plyranges::flank_upstream(100) %>%
  #   dplyr::mutate(seq = get_seq(., genome = "hg38", as_rna = TRUE))

  expect_s4_class(result, "GRanges")
  expect_true("seq" %in% names(GenomicRanges::mcols(result)))
  expect_s4_class(result$seq, "RNAStringSet")
  expect_equal(length(result$seq), length(result))
})