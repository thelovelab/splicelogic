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

# Test that find_si is find_ri with the up and down transcripts swapped
test_that("find_si reports the intron spliced out by an up transcript", {
  # split has the 16-20 intron, fused retains it in its 11-25 exon
  make_gr <- function(split_estimate) {
    dplyr::bind_rows(
      data.frame(
        start = c(1, 11, 21, 31), end = c(5, 15, 25, 35),
        tx_id = "split", estimate = split_estimate
      ),
      data.frame(
        start = c(1, 11, 31), end = c(5, 25, 35),
        tx_id = "fused", estimate = -split_estimate
      )
    ) |>
      dplyr::mutate(seqnames = "chr1", strand = "+", gene_id = "g1") |>
      dplyr::group_by(tx_id) |>
      dplyr::mutate(exon_rank = dplyr::row_number()) |>
      dplyr::ungroup() |>
      plyranges::as_granges() |>
      preprocess(coef_col = "estimate")
  }

  # split is up: it splices out the intron that fused (down) retains
  gr <- make_gr(1)
  si <- find_si(gr)
  expect_equal(length(si), 1L)
  expect_equal(as.integer(GenomicRanges::start(si)), 11L)
  expect_equal(as.integer(GenomicRanges::end(si)), 25L)
  expect_equal(si$tx_id, "fused")
  expect_equal(si$event_tx_id, "split")
  expect_equal(si$event_type, "si")
  expect_equal(si$event_estimate, 1)
  expect_equal(length(find_ri(gr)), 0L)
  # the exons either side of the intron are not alternative splice sites
  # of the fused exon
  expect_equal(length(find_a5ss(gr)), 0L)
  expect_equal(length(find_a3ss(gr)), 0L)
  expect_equal(find_all_events(gr, verbose = FALSE)$event_type, "si")

  # split is down: the same exon is now a retained intron instead
  gr <- make_gr(-1)
  ri <- find_ri(gr)
  expect_equal(length(ri), 1L)
  expect_equal(as.integer(GenomicRanges::start(ri)), 11L)
  expect_equal(ri$tx_id, "fused")
  expect_equal(ri$event_type, "ri")
  expect_equal(length(find_si(gr)), 0L)
  expect_equal(find_all_events(gr, verbose = FALSE)$event_type, "ri")

  # flipping every sign turns each retained intron into a spliced intron
  gr <- create_mock_data(3, 6, 3) |> generate_ri(n_events = 3)
  ri <- find_ri(gr)
  si <- find_si(gr |> dplyr::mutate(estimate = -estimate))
  expect_gt(length(ri), 0L)
  expect_equal(GenomicRanges::granges(si), GenomicRanges::granges(ri))
  expect_equal(si$tx_id, ri$tx_id)
  expect_equal(si$event_tx_id, ri$event_tx_id)
  expect_true(all(si$event_type == "si"))
})

# an exon that retains an intron can also move its outer boundary, as in
# TMEM259-201 vs -209 (ENST00000333175 / ENST00000592590). only the outer
# boundaries of the group are compared, so the moved one is an a5ss / a3ss
# and the ones facing the intron are left to find_ri() / find_si().
test_that("a retained intron with a moved outer boundary is also alt ss", {
  # split has the 16-20 intron; fused retains it in 13-25, whose start
  # also moved from 11 to 13
  make_gr <- function(split_estimate, strand) {
    dplyr::bind_rows(
      data.frame(
        start = c(1, 11, 21, 31), end = c(5, 15, 25, 35),
        tx_id = "split", estimate = split_estimate
      ),
      data.frame(
        start = c(1, 13, 31), end = c(5, 25, 35),
        tx_id = "fused", estimate = -split_estimate
      )
    ) |>
      dplyr::mutate(seqnames = "chr1", strand = strand, gene_id = "g1") |>
      dplyr::group_by(tx_id) |>
      dplyr::mutate(
        exon_rank = dplyr::if_else(
          strand == "-",
          dplyr::row_number(dplyr::desc(start)),
          dplyr::row_number(start)
        )
      ) |>
      dplyr::ungroup() |>
      plyranges::as_granges() |>
      preprocess(coef_col = "estimate")
  }

  # the start is the acceptor on "+" and the donor on "-"
  for (strand in c("+", "-")) {
    moved <- if (strand == "+") "a3ss" else "a5ss"
    other <- if (strand == "+") "a5ss" else "a3ss"
    find_moved <- if (strand == "+") find_a3ss else find_a5ss
    find_other <- if (strand == "+") find_a5ss else find_a3ss

    # split is up: fused splices the intron back in (si), and split's
    # 11-15 uses a different outer site than fused's 13-25
    gr <- make_gr(1, strand)
    alt <- find_moved(gr)
    expect_equal(length(alt), 1L)
    expect_equal(as.integer(GenomicRanges::start(alt)), 11L)
    expect_equal(as.integer(GenomicRanges::end(alt)), 15L)
    expect_equal(alt$tx_id, "split")
    expect_equal(alt$event_tx_id, "fused")
    expect_equal(alt$event_type, moved)
    # 15 and 21 face the intron: not alternative splice sites
    expect_equal(length(find_other(gr)), 0L)
    expect_equal(
      sort(find_all_events(gr, verbose = FALSE)$event_type),
      sort(c("si", moved))
    )

    # fused is up: the same change from the other side
    gr <- make_gr(-1, strand)
    alt <- find_moved(gr)
    expect_equal(length(alt), 1L)
    expect_equal(as.integer(GenomicRanges::start(alt)), 13L)
    expect_equal(as.integer(GenomicRanges::end(alt)), 25L)
    expect_equal(alt$tx_id, "fused")
    expect_equal(alt$event_tx_id, "split")
    expect_equal(alt$event_type, moved)
    expect_equal(length(find_other(gr)), 0L)
    expect_equal(
      sort(find_all_events(gr, verbose = FALSE)$event_type),
      sort(c("ri", moved))
    )
  }
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

# Test that an overlapping gene cannot hide an alternative splice site
test_that("find_a5ss and find_a3ss compare exons within the same gene", {
  # gene A has a real a5ss: posA moves the donor of the 21-25 exon to 23.
  # gene B is a readthrough over A whose negative transcript reuses the
  # shifted 21-23 exon exactly. the identity check must stay within each
  # gene, or the 21-23 candidate is dropped as "shared" with gene B.
  gene_a <- data.frame(
    start = c(1, 11, 21, 31, 41),
    end = c(5, 15, 25, 35, 45),
    exon_rank = 1:5,
    gene_id = "A"
  )
  gr <- dplyr::bind_rows(
    gene_a |> dplyr::mutate(tx_id = "negA", estimate = -1),
    gene_a |>
      dplyr::mutate(
        tx_id = "posA",
        estimate = 1,
        end = replace(end, 3, 23)
      ),
    data.frame(
      start = c(21, 31, 51, 61),
      end = c(23, 35, 55, 65),
      exon_rank = 1:4,
      gene_id = "B",
      tx_id = "negB",
      estimate = -1
    )
  ) |>
    dplyr::mutate(seqnames = "chr1", strand = "+") |>
    plyranges::as_granges() |>
    preprocess(coef_col = "estimate")


  a5ss <- find_a5ss(gr)
  expect_equal(length(a5ss), 1L)
  expect_equal(as.integer(GenomicRanges::start(a5ss)), 21L)
  expect_equal(as.integer(GenomicRanges::end(a5ss)), 23L)
  expect_equal(as.character(a5ss$tx_id), "posA")
  expect_equal(as.character(a5ss$event_tx_id), "negA")
  expect_equal(length(find_a3ss(gr)), 0L)
})

# Test for find_atss
test_that("find_atss test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_atss(gr, n_events = 3)

  result <- find_atss(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "atss"))
})

# Test for find_ates
test_that("find_ates test", {
  gr <- create_mock_data(3, 3, 6)
  gr <- preprocess(gr, coef_col = "estimate")
  gr <- generate_ates(gr, n_events = 3)

  result <- find_ates(gr)

  expect_s4_class(result, "GRanges")
  expect_true("event_type" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(result$event_type == "ates"))
})

# Test the non-overlapping variant: the terminal exon is dropped, so the
# next exon along becomes terminal and does not overlap its partner
test_that("find_atss and find_ates detect mode = 'drop' events", {
  gr <- create_mock_data(3, 3, 6)
  gr <- generate_atss(gr, n_events = 3, mode = "drop")
  expect_true(any(find_atss(gr)$event_type == "atss"))

  gr <- create_mock_data(3, 3, 6)
  gr <- generate_ates(gr, n_events = 3, mode = "drop")
  expect_true(any(find_ates(gr)$event_type == "ates"))
})

# Test that the two are not confused with each other
test_that("a moved start is not reported as a moved end", {
  gr <- create_mock_data(3, 3, 6)
  gr <- generate_atss(gr, n_events = 3)
  expect_gt(length(find_atss(gr)), 0L)
  expect_equal(length(find_ates(gr)), 0L)

  gr <- create_mock_data(3, 3, 6)
  gr <- generate_ates(gr, n_events = 3)
  expect_gt(length(find_ates(gr)), 0L)
  expect_equal(length(find_atss(gr)), 0L)
})

# Test that the labels follow the strand: on "-" the first exon by rank
# is the rightmost one, so its transcription start is its genomic end
test_that("find_atss and find_ates work on the minus strand", {
  gr <- create_mock_data(3, 3, 6, strand = "-")
  gr <- generate_atss(gr, n_events = 3)
  expect_true(any(find_atss(gr)$event_type == "atss"))
  expect_equal(length(find_ates(gr)), 0L)

  gr <- create_mock_data(3, 3, 6, strand = "-")
  gr <- generate_ates(gr, n_events = 3)
  expect_true(any(find_ates(gr)$event_type == "ates"))
  expect_equal(length(find_atss(gr)), 0L)
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
  # seven "-" strand genes, one engineered event each. exon_rank runs right
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

  # g7: the roles are reversed, so neg7 retains the 616-620 intron of pos7
  si <- find_si(gr)
  expect_s4_class(si, "GRanges")
  expect_equal(length(si), 1L)
  expect_equal(GenomicRanges::start(si), 611L)
  expect_equal(GenomicRanges::end(si), 625L)
  expect_equal(si$tx_id, "neg7")
  expect_equal(si$event_tx_id, "pos7")
  expect_equal(si$event_type, "si")

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
  # genomically the right-hand one. "over" and "identical" are symmetric
  # and the filters are on abs(l - r); "boundary" compares the inner boundary,
  # which on "-" is the start of the left flank and the end of the right
  # one. all three types must agree.
  counts <- vapply(
    c("boundary", "over", "identical"),
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
      "identical" = c(se = 1L, ie = 1L, mxe = 2L)
    )
  )

  # and the wrapper picks up all seven types, with nothing extra
  all_events <- find_all_events(gr, verbose = FALSE)
  expect_s4_class(all_events, "GRanges")
  expect_equal(length(all_events), 8L)
  expect_equal(
    sort(unique(all_events$event_type)),
    c("a3ss", "a5ss", "ie", "mxe", "ri", "se", "si")
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

test_that("type = 'boundary' matches the flank's inner boundary", {
  # down skips nothing; its 21-30 exon is the candidate. up_inner shares
  # the boundaries facing the candidate (10 and 41) but not the outer ones,
  # up_outer the reverse. only up_inner is a boundary match, on both
  # strands: on "-" the rank k - 1 flank is 41-50 and its inner boundary
  # is its start, so ignoring the strand would pick up_outer instead
  make_gr <- function(tx_strand) {
    dplyr::bind_rows(
      data.frame(
        start = c(1, 21, 41), end = c(10, 30, 50),
        tx_id = "down", estimate = -1
      ),
      data.frame(
        start = c(3, 41), end = c(10, 48), tx_id = "up_inner", estimate = 1
      ),
      data.frame(
        start = c(1, 39), end = c(12, 50), tx_id = "up_outer", estimate = 1
      )
    ) |>
      dplyr::mutate(seqnames = "chr1", strand = tx_strand, gene_id = "g1") |>
      dplyr::group_by(tx_id) |>
      # exon_rank follows transcription: rank 1 is the rightmost on "-"
      dplyr::arrange(
        if (tx_strand == "-") -start else start,
        .by_group = TRUE
      ) |>
      dplyr::mutate(exon_rank = dplyr::row_number()) |>
      dplyr::ungroup() |>
      plyranges::as_granges() |>
      preprocess(coef_col = "estimate")
  }

  for (tx_strand in c("+", "-")) {
    gr <- make_gr(tx_strand)

    boundary <- find_se(gr, type = "boundary")
    expect_equal(length(boundary), 1L)
    expect_equal(as.integer(GenomicRanges::start(boundary)), 21L)
    expect_equal(boundary$event_tx_id, "up_inner")

    # "over" (the default) takes both, "identical" neither
    expect_equal(sort(find_se(gr)$event_tx_id), c("up_inner", "up_outer"))
    expect_equal(length(find_se(gr, type = "identical")), 0L)
  }
})

test_that("extra arguments are passed on to findOverlaps", {
  # up's first exon (13-15) stops 2 bp short of down's left flank (1-10),
  # so the flank only matches once maxgap covers the gap
  gr <- dplyr::bind_rows(
    data.frame(
      start = c(1, 21, 41), end = c(10, 30, 50), tx_id = "down", estimate = -1
    ),
    data.frame(start = c(13, 41), end = c(15, 50), tx_id = "up", estimate = 1)
  ) |>
    dplyr::mutate(seqnames = "chr1", strand = "+", gene_id = "g1") |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(exon_rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    plyranges::as_granges() |>
    preprocess(coef_col = "estimate")

  expect_equal(length(find_se(gr)), 0L)
  se <- find_se(gr, maxgap = 2L)
  expect_equal(length(se), 1L)
  expect_equal(as.integer(GenomicRanges::start(se)), 21L)

  # and through the wrapper
  expect_false("se" %in% find_all_events(gr, verbose = FALSE)$event_type)
  expect_true(
    "se" %in% find_all_events(gr, verbose = FALSE, maxgap = 2L)$event_type
  )
})
