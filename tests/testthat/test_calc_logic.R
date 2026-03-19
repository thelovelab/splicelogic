# Test for calc_skipped_exons from calc_logic.R
test_that("calc_skipped_exons works with se_mock_data and preprocess", {
  gr <- se_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr)
  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "skipped_exon"))
})

# # Test for invalid coef column input in calc_skipped_exons
# test_that("calc_skipped_exons errors if coef_col is invalid", {
#   gr <- se_mock_data()
#   expect_error(
#     calc_skipped_exons(gr),
#     regexp = "Missing required metadata columns: foo"
#   )
# })

# Test for no event detected in skipped exon
test_that("calc_skipped_exons returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for skipped exon detection in mx_mock_data
test_that("calc_skipped_exons detects single skipped exon in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), 5L)
  expect_equal(as.character(GenomicRanges::mcols(result)$event), "skipped_exon")
})


# Test for mutually exclusive detection in mx_mock_data
test_that("calc_mx_exons detects mutually exclusive exons in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- calc_mx_exons(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 4L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), c(1, 2, 1, 2))
  expect_equal(
    as.integer(GenomicRanges::mcols(result)$exon_rank),
    c(3, 3, 8, 7)
  )
  expect_equal(
    as.character(GenomicRanges::mcols(result)$event),
    rep("mutually_exclusive", 4L)
  )
})

# Test for no event detected in mx
test_that("calc_mx_exons returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "coefs")
  result <- calc_mx_exons(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for calc retained introns
test_that("calc_retained_introns test", {
  gr <- create_mock_data(3, 6, 3)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_retained_introns(gr, n_events = 3)

  result <- calc_retained_introns(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "retained_intron"))
})

# Test for calc_a3ss
test_that("calc_a3ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_a3ss(gr, n_events = 3)

  result <- calc_a3ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "a3ss"))
})


# Test for calc_a5ss
test_that("calc_a5ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_a5ss(gr, n_events = 3)

  result <- calc_a5ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "a5ss"))
})

# Test for calc_skipped_exons with new mock data
test_that("calc_skipped_exons test with generate_skipped_exons", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_skipped_exons(gr, n_events = 1)

  result <- calc_skipped_exons(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "skipped_exon"))
})

# Test for calc_all_events
test_that("calc_all_events returns combined GRanges with multiple event types", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_skipped_exons(gr, n_events = 1)
  gr <- generate_a3ss(gr, n_events = 1)
  gr <- generate_a5ss(gr, n_events = 1)

  result <- calc_all_events(gr, type = "over", verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_true(length(result) > 0L)
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true("skipped_exon" %in% result$event)
})

test_that("calc_all_events returns empty GRanges when no events exist", {
  gr <- no_event_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")

  result <- calc_all_events(gr, verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})
