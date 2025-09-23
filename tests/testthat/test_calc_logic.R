# Test for calc_skipped_exons from calc_logic.R
test_that("calc_skipped_exons works with get_mock_data and preprocess_input", {
  gr <- get_mock_data()
  gr <- preprocess_input(gr, coef_col = "coefs")
  result <- calc_skipped_exons(gr, coef_col = "coefs")
  expect_s4_class(result, "GRanges")
  expect_true("event" %in% names(mcols(result)))
  expect_true(any(result$event == "skipped_exon"))
})
# Test for invalid coef column input
test_that("calc_skipped_exons errors if coef_col is invalid", {
  gr <- get_mock_data()
  expect_error(
    calc_skipped_exons(gr, coef_col = "foo"),
    regexp = "Missing required metadata columns: foo"
  )
})
