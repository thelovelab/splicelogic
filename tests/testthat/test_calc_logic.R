
# Test for calc_skipped_exons from calc_logic.R
test_that("calc_skipped_exons works with get_mock_data and preprocess_input", {
  gr <- get_mock_data()
  gr <- preprocess_input(gr)
  result <- calc_skipped_exons(gr)
  expect_true("event" %in% names(mcols(result)))
  expect_true(any(result$event == "skipped_exon"))
})

