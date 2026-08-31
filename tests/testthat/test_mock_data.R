# Test for se_mock_data
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

# Test for mx_mock_data
test_that("mx_mock_data returns GRanges with required columns", {
  gr <- mx_mock_data()
  expect_s4_class(gr, "GRanges")
  required_cols <- c("exon_rank", "gene_id", "tx_id", "estimate")
  expect_true(all(required_cols %in% names(GenomicRanges::mcols(gr))))
  expect_equal(length(unique(gr$tx_id)), 2L)
  expect_true(any(gr$estimate < 0) && any(gr$estimate > 0))
})

# Test for no_event_mock_data
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

# Test for create_mock_data
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

test_that("create_mock_data genes do not overlap when n_exons_per_tx is large", {
  set.seed(42)
  gr <- create_mock_data(n_genes = 5, n_tx_per_gene = 3, n_exons_per_tx = 20)
  gr_se <- generate_se(gr, n_events = 1)
  result <- find_se(gr_se)
  expect_gt(length(result), 0L)
})

# Test for generate_se
test_that("generate_se removes exons", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_se <- generate_se(gr, n_events = 1)
  expect_s4_class(gr_se, "GRanges")
  expect_true(length(gr_se) < original_len)
})

# Test for generate_mxe
test_that("generate_mxe removes exons pairwise", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_mx <- generate_mxe(gr, n_events = 1)
  expect_s4_class(gr_mx, "GRanges")
  expect_true(length(gr_mx) < original_len)
})

# Test for generate_ri
test_that("generate_ri merges exons", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  original_len <- length(gr)
  gr_ri <- generate_ri(gr, n_events = 1)
  expect_s4_class(gr_ri, "GRanges")
  expect_true(length(gr_ri) < original_len)
})

# Test for generate_a5ss
test_that("generate_a5ss modifies exon ends", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  gr_a5 <- generate_a5ss(gr, n_events = 1)
  expect_s4_class(gr_a5, "GRanges")
  expect_equal(length(gr_a5), length(gr))
  expect_false(all(GenomicRanges::end(gr_a5) == GenomicRanges::end(gr)))
})

# Test for generate_a3ss
test_that("generate_a3ss modifies exon starts", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6)
  gr_a3 <- generate_a3ss(gr, n_events = 1)
  expect_s4_class(gr_a3, "GRanges")
  expect_equal(length(gr_a3), length(gr))
  expect_false(all(GenomicRanges::start(gr_a3) == GenomicRanges::start(gr)))
})

# Test the neg/pos pairing guarantee and its one boundary case
test_that("a single transcript per gene warns and is all-negative", {
  # the neg/pos pairing keys off tx_order 1 and 2, so one transcript per
  # gene leaves everything negative and nothing for a generator to edit
  expect_warning(gr <- create_mock_data(2, 1, 5), "n_tx_per_gene < 2")
  expect_true(all(GenomicRanges::mcols(gr)$estimate < 0))

  # two transcripts per gene is fine and does not warn
  expect_silent(gr2 <- create_mock_data(2, 2, 5))
  expect_true(any(GenomicRanges::mcols(gr2)$estimate > 0))
  expect_s4_class(generate_ri(gr2, 1), "GRanges")
})

# Test that create_mock_data honours the strand argument
test_that("create_mock_data places transcripts on the requested strand", {
  gr_plus <- create_mock_data(2, 3, 4, strand = "+")
  gr_minus <- create_mock_data(2, 3, 4, strand = "-")

  expect_true(all(as.character(GenomicRanges::strand(gr_plus)) == "+"))
  expect_true(all(as.character(GenomicRanges::strand(gr_minus)) == "-"))
  expect_equal(length(gr_minus), length(gr_plus))

  # on "-" exon_rank runs the other way: rank 1 is the rightmost exon,
  # so rank and genomic start are anti-correlated within a transcript
  rank_vs_start <- function(gr) {
    tibble::as_tibble(gr) |>
      dplyr::group_by(tx_id) |>
      dplyr::summarise(rho = stats::cor(exon_rank, start), .groups = "drop") |>
      dplyr::pull(rho)
  }
  expect_true(all(rank_vs_start(gr_plus) == 1))
  expect_true(all(rank_vs_start(gr_minus) == -1))

  expect_error(create_mock_data(strand = "?"))
})

# Test that the generators put the event on the strand-correct boundary
test_that("generate_a5ss and generate_a3ss follow the strand", {
  # on "+" the donor is the exon end and the acceptor the exon start;
  # on "-" the two swap. shifting the wrong boundary silently turns an
  # a5ss into an a3ss and vice versa, which is what this guards.
  for (st in c("+", "-")) {
    set.seed(7)
    gr <- create_mock_data(3, 3, 6, strand = st)

    a5 <- generate_a5ss(gr, n_events = 2)
    expect_gt(length(find_a5ss(a5)), 0L)
    expect_equal(length(find_a3ss(a5)), 0L)

    a3 <- generate_a3ss(gr, n_events = 2)
    expect_gt(length(find_a3ss(a3)), 0L)
    expect_equal(length(find_a5ss(a3)), 0L)

    # the moved boundary is the one the strand marks as the splice site
    moved_end <- GenomicRanges::end(a5) != GenomicRanges::end(gr)
    moved_start <- GenomicRanges::start(a5) != GenomicRanges::start(gr)
    if (st == "+") {
      expect_true(any(moved_end) && !any(moved_start))
    } else {
      expect_true(any(moved_start) && !any(moved_end))
    }
  }
})

test_that("generate_se, mxe and ri work on the minus strand", {
  # these three are rank-based rather than coordinate-based, but
  # generate_ri merges an exon with its rank neighbour, which lies to the
  # genomic left on "-" -- extending the end there would invert the range
  for (st in c("+", "-")) {
    set.seed(3)
    gr <- create_mock_data(3, 3, 6, strand = st)

    se <- generate_se(gr, n_events = 2)
    expect_gt(length(find_se(se)), 0L)

    mxe <- generate_mxe(gr, n_events = 2)
    expect_gt(length(find_mxe(mxe)), 0L)

    ri <- generate_ri(gr, n_events = 2)
    expect_true(all(GenomicRanges::width(ri) > 0))
    expect_gt(length(find_ri(ri)), 0L)
  }
})
