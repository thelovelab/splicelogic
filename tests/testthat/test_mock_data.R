
test_that("get_mock_data returns a GRanges object", {
  gr <- get_mock_data()
  expect_s4_class(gr, "GRanges")
})

test_that("get_mock_data returns correct number of ranges", {
  gr <- get_mock_data()
  expect_equal(length(gr), 7) # 4 from df1, 3 from df2
})

test_that("get_mock_data returns correct tx_id values", {
  gr <- get_mock_data()
  expect_equal(gr$tx_id, c(rep(1, 4), rep(2, 3)))
})
