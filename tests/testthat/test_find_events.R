# Test for find_se from calc_logic.R
test_that("find_se works with se_mock_data and preprocess", {
  gr <- se_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- find_se(gr)
  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "se"))
})

# # Test for invalid coef column input in find_se
# test_that("find_se errors if coef_col is invalid", {
#   gr <- se_mock_data()
#   expect_error(
#     find_se(gr),
#     regexp = "Missing required metadata columns: foo"
#   )
# })

# Test for no event detected in skipped exon
test_that("find_se returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "coefs")
  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for skipped exon detection in mx_mock_data
test_that("find_se detects single skipped exon in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), 5L)
  expect_equal(as.character(GenomicRanges::mcols(result)$event), "se")
})


# Test for mutually exclusive detection in mx_mock_data
test_that("find_mxe detects mutually exclusive exons in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")
  result <- find_mxe(gr)

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
    rep("mxe", 4L)
  )
})

# Test for no event detected in mx
test_that("find_mxe returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "coefs")
  result <- find_mxe(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for calc retained introns
test_that("find_ri test", {
  gr <- create_mock_data(3, 6, 3)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_retained_introns(gr, n_events = 3)

  result <- find_ri(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "ri"))
})

# Test for find_a3ss
test_that("find_a3ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_a3ss(gr, n_events = 3)

  result <- find_a3ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "a3ss"))
})


# Test for find_a5ss
test_that("find_a5ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_a5ss(gr, n_events = 3)

  result <- find_a5ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "a5ss"))
})

# Test for find_se with new mock data
test_that("find_se test with generate_skipped_exons", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_skipped_exons(gr, n_events = 1)

  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event == "se"))
})

# Test for find_all_events
test_that("find_all_events returns combined GRanges with multiple event types", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "coefs")
  gr <- generate_skipped_exons(gr, n_events = 1)
  gr <- generate_a3ss(gr, n_events = 1)
  gr <- generate_a5ss(gr, n_events = 1)

  result <- find_all_events(gr, type = "over", verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_true(length(result) > 0L)
  expect_true("event" %in% names(GenomicRanges::mcols(result)))
  expect_true("se" %in% result$event)
})

test_that("find_all_events returns empty GRanges when no events exist", {
  gr <- no_event_mock_data()
  gr <- preprocess(gr, coef_col = "coefs")

  result <- find_all_events(gr, verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})
