# test_that
# --- se_mock_data ---
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

# --- mx_mock_data ---
test_that("mx_mock_data returns GRanges with required columns", {
  gr <- mx_mock_data()
  expect_s4_class(gr, "GRanges")
  required_cols <- c("exon_rank", "gene_id", "tx_id", "estimate")
  expect_true(all(required_cols %in% names(GenomicRanges::mcols(gr))))
  expect_equal(length(unique(gr$tx_id)), 2L)
  expect_true(any(gr$estimate < 0) && any(gr$estimate > 0))
})

# --- no_event_mock_data ---
test_that("no_event_mock_data returns GRanges with matching exon structure", {
  gr <- no_event_mock_data()
  expect_s4_class(gr, "GRanges")
  required_cols <- c("exon_rank", "gene_id", "tx_id", "estimate")
  expect_true(all(required_cols %in% names(GenomicRanges::mcols(gr))))
  expect_equal(length(unique(gr$tx_id)), 2L)
  tx1 <- sort(GenomicRanges::start(gr[gr$tx_id == 1]))
  tx2 <- sort(GenomicRanges::start(gr[gr$tx_id == 2]))
  expect_equal(tx1, tx2)
})

# --- create_mock_data ---
test_that("create_mock_data returns correct dimensions", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 4)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2 * 3 * 4)
  expect_equal(length(unique(gr$gene_id)), 2L)
  expect_equal(length(unique(gr$tx_id)), 6L)
})

test_that("create_mock_data ensures pos and neg estimates per gene", {
  gr <- create_mock_data(n_genes = 3, n_tx_per_gene = 4, n_exons_per_tx = 5)
  for (gid in unique(gr$gene_id)) {
    gene_estimates <- unique(gr$estimate[gr$gene_id == gid])
    expect_true(any(gene_estimates < 0))
    expect_true(any(gene_estimates > 0))
  }
})

test_that("create_mock_data is preprocessed", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 5)
  preprocess_cols <- c("key", "nexons", "internal", "estimate")
  expect_true(all(preprocess_cols %in% names(GenomicRanges::mcols(gr))))
})

# --- generate_se ---
test_that("generate_se removes exons", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_se <- generate_se(gr, n_events = 1)
  expect_s4_class(gr_se, "GRanges")
  expect_true(length(gr_se) < original_len)
})

# --- generate_mxe ---
test_that("generate_mxe removes exons pairwise", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_mx <- generate_mxe(gr, n_events = 1)
  expect_s4_class(gr_mx, "GRanges")
  expect_true(length(gr_mx) < original_len)
})

# --- generate_ri ---
test_that("generate_ri merges exons", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_ri <- generate_ri(gr, n_events = 1)
  expect_s4_class(gr_ri, "GRanges")
  expect_true(length(gr_ri) < original_len)
})

# --- generate_a5ss ---
test_that("generate_a5ss modifies exon ends", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  gr_a5 <- generate_a5ss(gr, n_events = 1)
  expect_s4_class(gr_a5, "GRanges")
  expect_equal(length(gr_a5), length(gr))
  expect_false(all(GenomicRanges::end(gr_a5) == GenomicRanges::end(gr)))
})

# --- generate_a3ss ---
test_that("generate_a3ss modifies exon starts", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  gr_a3 <- generate_a3ss(gr, n_events = 1)
  expect_s4_class(gr_a3, "GRanges")
  expect_equal(length(gr_a3), length(gr))
  expect_false(all(GenomicRanges::start(gr_a3) == GenomicRanges::start(gr)))
})
