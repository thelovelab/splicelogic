
# test_that 
# tests for se_mock_data from mock_data.R
test_that("se_mock_data returns a GRanges and has required metadata columns", {
  gr <- se_mock_data()
  expect_s4_class(gr, "GRanges")

  required_cols <- c("exon_rank", "gene_id", "tx_id", "coefs")
  expect_true(all(required_cols %in% names(GenomicRanges::mcols(gr))))
})

test_that("coefs metadata is numeric and within [-1, 1]", {
  gr <- se_mock_data()
  coefs <- GenomicRanges::mcols(gr)$coefs
  expect_true(is.numeric(coefs))
  expect_true(all(!is.na(coefs)))
  expect_true(all(coefs >= -1 & coefs <= 1))
})

test_that("check metadata column types", {
  gr <- se_mock_data()
  md <- GenomicRanges::mcols(gr)
  expect_true(is.numeric(md$exon_rank) || is.integer(md$exon_rank))
  expect_true(is.numeric(md$tx_id) || is.integer(md$tx_id))
  expect_true(is.numeric(md$gene_id) || is.integer(md$gene_id))
})