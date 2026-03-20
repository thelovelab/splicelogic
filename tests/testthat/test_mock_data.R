# test_that
# tests for se_mock_data from mock_data.R
test_that("se_mock_data returns a GRanges and has required metadata columns", {
  gr <- se_mock_data()
  expect_s4_class(gr, "GRanges")

  required_cols <- c("exon_rank", "gene_id", "tx_id", "estimate")
  expect_true(all(required_cols %in% names(GenomicRanges::mcols(gr))))
})

test_that("estimate metadata is numeric and within [-1, 1]", {
  gr <- se_mock_data()
  estimate <- GenomicRanges::mcols(gr)$estimate
  expect_true(is.numeric(estimate))
  expect_true(all(!is.na(estimate)))
  expect_true(all(estimate >= -1 & estimate <= 1))
})

test_that("check metadata column types", {
  gr <- se_mock_data()
  md <- GenomicRanges::mcols(gr)
  expect_true(is.numeric(md$exon_rank) || is.integer(md$exon_rank))
  expect_true(is.numeric(md$tx_id) || is.integer(md$tx_id))
  expect_true(is.numeric(md$gene_id) || is.integer(md$gene_id))
})
