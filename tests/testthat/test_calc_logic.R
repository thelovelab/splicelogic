# Test for calc_skipped_exons from calc_logic.R
test_that("calc_skipped_exons works with se_mock_data and preprocess_input", {
  gr <- se_mock_data()
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr, coef_col = "coefs")
  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(mcols(result)))
  expect_true(any(result$event == "skipped_exon"))
})

# Test for invalid coef column input in calc_skipped_exons
test_that("calc_skipped_exons errors if coef_col is invalid", {
  gr <- se_mock_data()
  expect_error(
    calc_skipped_exons(gr, coef_col = "foo"),
    regexp = "Missing required metadata columns: foo"
  )
})

# Test for no event detected in skipped exon
test_that("calc_mutually_exclusive returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr, coef_col = "coefs")

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for skipped exon detection in mx_mock_data
test_that("calc_skipped_exons detects single skipped exon in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr, coef_col = "coefs")

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), 5L)
  expect_equal(as.character(GenomicRanges::mcols(result)$event), "skipped_exon")
})


# Test for mutually exclusive detection in mx_mock_data
test_that("calc_mutually_exclusive detects mutually exclusive exons in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_mutually_exclusive(gr, coef_col = "coefs")

  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 4L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), c(1, 2, 1, 2))
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), c(3, 3, 8, 7))
  expect_equal(as.character(GenomicRanges::mcols(result)$event), rep("mutually_exclusive", 4L))
})

# Test for no event detected in mx 
test_that("calc_mutually_exclusive returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_mutually_exclusive(gr, coef_col = "coefs")

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for calc retained introns
test_that("calc_retained_introns test", {
    gr <- create_mock_data(3,6,3)
    gr <- preprocess_input(gr, coef_col = "coefs")
    gr <- generate_retained_introns(gr, n_ri = 3)

    result <- calc_retained_introns(gr, coef_col = "coefs")

    expect_s4_class(result, "GRanges")
    expect_true("event" %in% names(GenomicRanges::mcols(result)))
    expect_true(any(result$event == "retained_intron"))

})

# Test for calc_a3ss_a5ss
test_that("calc_a3ss_a5ss test", {
    gr <- create_mock_data(1,2,4)
    gr <- preprocess_input(gr, coef_col = "coefs")
    gr <- generate_a3ss(gr, n_a3ss = 3)

    result <- calc_a3ss_a5ss(gr, coef_col = "coefs")

    expect_s4_class(result, "GRanges")
    expect_true("event" %in% names(GenomicRanges::mcols(result)))
    expect_true(any(result$event == "a3ss"))
})
