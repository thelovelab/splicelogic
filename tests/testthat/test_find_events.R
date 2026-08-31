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

# Test that the a5ss/a3ss labels follow the strand
test_that("find_a5ss and find_a3ss respect the minus strand", {
  # the same five exons, but read right to left: exon_rank 1 is the 41-45
  # exon, so the exon end is the acceptor and the exon start is the donor,
  # the opposite of the "+" strand case above. a shifted boundary that
  # would be an a3ss on "+" is therefore an a5ss here.
  base <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41),
    end = c(5, 15, 25, 35, 45),
    strand = "-",
    gene_id = 1
  )
  # tx 1 is the anchor, tx 2 carries the shifted exon
  make_gr <- function(alt) {
    dplyr::bind_rows(
      base |> dplyr::mutate(tx_id = 1, estimate = -1, exon_rank = 5:1),
      alt |> dplyr::mutate(tx_id = 2, estimate = 1, exon_rank = 5:1)
    ) |>
      plyranges::as_granges() |>
      preprocess(coef_col = "estimate")
  }

  # the donor moves (21-25 -> 23-25): on "-" that is the exon start
  gr_a5 <- make_gr(base |> dplyr::mutate(start = replace(start, 3, 23)))
  a5ss <- find_a5ss(gr_a5)
  expect_equal(length(a5ss), 1L)
  expect_equal(as.integer(GenomicRanges::start(a5ss)), 23L)
  expect_equal(GenomicRanges::mcols(a5ss)$event_type, "a5ss")
  expect_equal(length(find_a3ss(gr_a5)), 0L)

  # the acceptor moves (21-25 -> 21-23): on "-" that is the exon end
  gr_a3 <- make_gr(base |> dplyr::mutate(end = replace(end, 3, 23)))
  a3ss <- find_a3ss(gr_a3)
  expect_equal(length(a3ss), 1L)
  expect_equal(as.integer(GenomicRanges::end(a3ss)), 23L)
  expect_equal(GenomicRanges::mcols(a3ss)$event_type, "a3ss")
  expect_equal(length(find_a5ss(gr_a3)), 0L)
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

# Test every finder against a minus-strand fixture
test_that("all find_* functions work on minus-strand data", {
  # six "-" strand genes, one engineered event each. exon_rank runs right
  # to left, so this checks that the finders key off exon_rank and strand
  # rather than off genomic order.
  gr <- preprocess(minus_strand_mock_data(), coef_col = "estimate")

  # g1: pos1 splices out 21-25, which is exon_rank 3 of neg1
  se <- find_se(gr)
  expect_s4_class(se, "GRanges")
  expect_equal(length(se), 1L)
  expect_equal(GenomicRanges::start(se), 21L)
  expect_equal(GenomicRanges::end(se), 25L)
  expect_equal(as.character(GenomicRanges::strand(se)), "-")
  expect_equal(se$tx_id, "neg1")
  expect_equal(se$exon_rank, 3L)
  expect_equal(se$event_tx_id, "pos1")
  expect_equal(se$event_type, "se")

  # g6: the roles are reversed, so 521-525 is present only in pos6
  ie <- find_ie(gr)
  expect_s4_class(ie, "GRanges")
  expect_equal(length(ie), 1L)
  expect_equal(GenomicRanges::start(ie), 521L)
  expect_equal(GenomicRanges::end(ie), 525L)
  expect_equal(ie$tx_id, "pos6")
  expect_equal(ie$event_tx_id, "neg6")
  expect_equal(ie$event_type, "ie")

  # g2: 121-125 (neg2) and 131-135 (pos2) sit between the same flanks and
  # are disjoint. results are interleaved candidate-then-partner.
  mxe <- find_mxe(gr)
  expect_s4_class(mxe, "GRanges")
  expect_equal(length(mxe), 2L)
  expect_equal(GenomicRanges::start(mxe), c(121L, 131L))
  expect_equal(GenomicRanges::end(mxe), c(125L, 135L))
  expect_equal(mxe$tx_id, c("neg2", "pos2"))
  expect_equal(mxe$event_tx_id, c("pos2", "neg2"))
  expect_equal(mxe$event_type, rep("mxe", 2L))

  # g3: pos3 retains the 216-220 intron of neg3
  ri <- find_ri(gr)
  expect_s4_class(ri, "GRanges")
  expect_equal(length(ri), 1L)
  expect_equal(GenomicRanges::start(ri), 211L)
  expect_equal(GenomicRanges::end(ri), 225L)
  expect_equal(ri$tx_id, "pos3")
  expect_equal(ri$event_tx_id, "neg3")
  expect_equal(ri$event_type, "ri")

  # g4: the exon start moves (311-315 -> 313-315). on "-" the start is the
  # donor, so this is an a5ss and must not be reported as an a3ss.
  a5ss <- find_a5ss(gr)
  expect_s4_class(a5ss, "GRanges")
  expect_equal(length(a5ss), 1L)
  expect_equal(GenomicRanges::start(a5ss), 313L)
  expect_equal(GenomicRanges::end(a5ss), 315L)
  expect_equal(a5ss$tx_id, "pos4")
  expect_equal(a5ss$event_tx_id, "neg4")
  expect_equal(a5ss$event_type, "a5ss")

  # g5: the exon end moves (411-415 -> 411-413), the acceptor on "-"
  a3ss <- find_a3ss(gr)
  expect_s4_class(a3ss, "GRanges")
  expect_equal(length(a3ss), 1L)
  expect_equal(GenomicRanges::start(a3ss), 411L)
  expect_equal(GenomicRanges::end(a3ss), 413L)
  expect_equal(a3ss$tx_id, "pos5")
  expect_equal(a3ss$event_tx_id, "neg5")
  expect_equal(a3ss$event_type, "a3ss")

  # the flank lookup labels the rank k - 1 exon "left", which on "-" is
  # genomically the right-hand one. nothing downstream depends on the
  # label -- every match type is symmetric and the filters are on
  # abs(l - r) -- so all three types must agree.
  counts <- vapply(
    c("boundary", "over", "in"),
    function(ty) {
      c(
        se = length(find_se(gr, type = ty)),
        ie = length(find_ie(gr, type = ty)),
        mxe = length(find_mxe(gr, type = ty))
      )
    },
    integer(3L)
  )
  expect_equal(
    counts,
    cbind(
      boundary = c(se = 1L, ie = 1L, mxe = 2L),
      over = c(se = 1L, ie = 1L, mxe = 2L),
      "in" = c(se = 1L, ie = 1L, mxe = 2L)
    )
  )

  # and the wrapper picks up all six types, with nothing extra
  all_events <- find_all_events(gr, verbose = FALSE)
  expect_s4_class(all_events, "GRanges")
  expect_equal(length(all_events), 7L)
  expect_equal(
    sort(unique(all_events$event_type)),
    c("a3ss", "a5ss", "ie", "mxe", "ri", "se")
  )
  expect_true(all(as.character(GenomicRanges::strand(all_events)) == "-"))
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
