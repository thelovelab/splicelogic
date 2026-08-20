# Test for find_se from calc_logic.R
test_that("find_se works with se_mock_data and preprocess", {
  gr <- se_mock_data()
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_se(gr)
  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "se"))
})

# Test for find_se when exon_rank is offset (CDS-style ranges)
test_that("find_se detects a skipped exon when exon_rank starts at 2", {
  gr <- cds_mock_data()
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_se(gr)

  expect_s4_class(result, "GRanges")

  # tx 2 skips 31-35, which is exon_rank 5 of tx 1. This is only found
  # if 'internal' is computed relative to each transcript's own rank
  # range: the candidate is the second-to-last exon, which the old
  # 'exon_rank < nexons' form marked terminal and discarded.
  expect_equal(length(result), 1L)
  expect_equal(as.integer(GenomicRanges::start(result)), 31L)
  expect_equal(as.integer(GenomicRanges::end(result)), 35L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), 5L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$event_tx_id), 2L)
  expect_equal(as.character(GenomicRanges::mcols(result)$event_type), "se")
})

# Test for no event detected in skipped exon
test_that("find_se returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for skipped exon detection in mx_mock_data
test_that("find_se detects single skipped exon in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$exon_rank), 5L)
  expect_equal(as.character(GenomicRanges::mcols(result)$event_type), "se")
})


# Test for mutually exclusive detection in mx_mock_data
test_that("find_mxe detects mutually exclusive exons in mx_mock_data", {
  gr <- mx_mock_data()
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_mxe(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))

  # expect exactly one detected event on exon_rank = 5 of tx_id 1
  expect_equal(length(result), 4L)
  expect_equal(as.integer(GenomicRanges::mcols(result)$tx_id), c(1, 2, 1, 2))
  expect_equal(
    as.integer(GenomicRanges::mcols(result)$exon_rank),
    c(3, 3, 8, 7)
  )
  expect_equal(
    as.character(GenomicRanges::mcols(result)$event_type),
    rep("mxe", 4L)
  )
})

# Test that an alternative splice site is not also reported as MX
test_that("find_mxe ignores a middle exon overlapping the candidate", {
  # tx 1 is the anchor. tx 2 drops exon 3 entirely, which is what makes
  # 21-25 a candidate at all. tx 3 keeps exon 3 but shifts its start
  # (21-25 -> 23-25): truncating an exon does not change how many exons sit
  # between the flanks, so the gap stays 2 and, without the disjointness
  # check, tx 3 is reported as an MX partner as well as an a3ss.
  base <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41),
    end = c(5, 15, 25, 35, 45),
    strand = "+",
    gene_id = 1
  )
  gr <- dplyr::bind_rows(
    base |> dplyr::mutate(tx_id = 1, estimate = -1, exon_rank = seq_len(5)),
    base[-3, ] |>
      dplyr::mutate(tx_id = 2, estimate = 1, exon_rank = seq_len(4)),
    base |>
      dplyr::mutate(
        tx_id = 3,
        estimate = 1,
        exon_rank = seq_len(5),
        start = replace(start, 3, 23)
      )
  ) |>
    plyranges::as_granges() |>
    preprocess(coef_col = "estimate")

  # tx 3's shifted exon is an alternative 3' splice site, not an MX pair
  expect_equal(length(find_mxe(gr)), 0L)

  a3ss <- find_a3ss(gr)
  expect_equal(length(a3ss), 1L)
  expect_equal(as.integer(GenomicRanges::start(a3ss)), 23L)
  expect_equal(as.integer(GenomicRanges::mcols(a3ss)$tx_id), 3L)

  # the genuine skipped exon in tx 2 is still detected
  se <- find_se(gr)
  expect_equal(length(se), 1L)
  expect_equal(as.integer(GenomicRanges::mcols(se)$event_tx_id), 2L)

  # the same holds when the donor moves instead of the acceptor
  # (21-25 -> 21-23): the check is on overlap, not on which boundary
  # shifted, so an a5ss is not reported as MX either
  gr_a5 <- dplyr::bind_rows(
    base |> dplyr::mutate(tx_id = 1, estimate = -1, exon_rank = seq_len(5)),
    base[-3, ] |>
      dplyr::mutate(tx_id = 2, estimate = 1, exon_rank = seq_len(4)),
    base |>
      dplyr::mutate(
        tx_id = 3,
        estimate = 1,
        exon_rank = seq_len(5),
        end = replace(end, 3, 23)
      )
  ) |>
    plyranges::as_granges() |>
    preprocess(coef_col = "estimate")

  expect_equal(length(find_mxe(gr_a5)), 0L)

  a5ss <- find_a5ss(gr_a5)
  expect_equal(length(a5ss), 1L)
  expect_equal(as.integer(GenomicRanges::end(a5ss)), 23L)
  expect_equal(as.integer(GenomicRanges::mcols(a5ss)$tx_id), 3L)
})

# Test for no event detected in mx
test_that("find_mxe returns empty GRanges if no events", {
  gr <- no_event_mock_data() # no_event_mock_data has no mx events
  gr <- preprocess(gr, coef_col = "estimate")
  result <- find_mxe(gr)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

# Test for calc retained introns
test_that("find_ri test", {
  gr <- create_mock_data(3, 6, 3)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_ri(gr, n_events = 3)

  result <- find_ri(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "ri"))
})

# Test for find_a3ss
test_that("find_a3ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_a3ss(gr, n_events = 3)

  result <- find_a3ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "a3ss"))
})


# Test for find_a5ss
test_that("find_a5ss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_a5ss(gr, n_events = 3)

  result <- find_a5ss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "a5ss"))
})

# Test for find_se with new mock data
test_that("find_se test with generate_se", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_se(gr, n_events = 1)

  result <- find_se(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "se"))
})

# Test for find_all_events
test_that("find_all_events returns combined GRanges with multiple event types", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_se(gr, n_events = 1)
  gr <- generate_a3ss(gr, n_events = 1)
  gr <- generate_a5ss(gr, n_events = 1)

  result <- find_all_events(gr, type = "over", verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_true(length(result) > 0L)
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true("se" %in% result$event_type)
})

test_that("find_all_events returns empty GRanges when no events exist", {
  gr <- no_event_mock_data()
  gr <- preprocess(gr, coef_col = "estimate")

  result <- find_all_events(gr, verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 0L)
})

test_that("find_se counts isoforms, not exons, when picking candidates", {
  gr <- preprocess(intron_split_mock_data(), coef_col = "estimate")

  result <- find_se(gr, type = "boundary")

  expect_s4_class(result, "GRanges")
  # exactly one event: the 100-300 exon of neg1, skipped by posC
  expect_equal(length(result), 1L)
  expect_equal(GenomicRanges::start(result), 100L)
  expect_equal(GenomicRanges::end(result), 300L)
  expect_equal(as.character(result$tx_id), "neg1")
  expect_equal(as.character(result$event_type), "se")
  # posA and posB keep the region (split in two), only posC skips it
  expect_equal(as.character(result$event_tx_id), "posC")
})
