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

test_that("every generator makes at most one event per gene", {
  # n_events is the number of genes that get an event, the same for every
  # generator, so asking for more than there are genes warns and caps
  set.seed(3)
  gr <- create_mock_data(n_genes = 3, n_tx_per_gene = 4, n_exons_per_tx = 7)
  bare <- gr |> dplyr::select(gene_id, tx_id, exon_rank, estimate)
  generators <- list(
    se = generate_se, ie = generate_ie, cse = generate_cse,
    cie = generate_cie, mxe = generate_mxe, ri = generate_ri,
    si = generate_si, a5ss = generate_a5ss, a3ss = generate_a3ss,
    atss = generate_atss, ates = generate_ates,
    atss_drop = function(gr, n_events) generate_atss(gr, n_events, "drop"),
    ates_drop = function(gr, n_events) generate_ates(gr, n_events, "drop")
  )
  n_event_genes <- function(out) length(unique(out$gene_id[out$sim_event]))
  for (ev in names(generators)) {
    gen <- generators[[ev]]
    expect_silent(out <- gen(gr, n_events = 2))
    expect_equal(n_event_genes(out), 2L, info = ev)

    expect_warning(out <- gen(gr, n_events = 5), "capping at 3", info = ev)
    expect_equal(n_event_genes(out), 3L, info = ev)

    # input that has not been through preprocess() is preprocessed first
    expect_true(any(gen(bare, n_events = 1)$sim_event), info = ev)
  }
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


# one row per transcript that has fewer exons in `out` than in `gr`, with
# its estimate (constant within a transcript) and how many exons it lost
exons_lost <- function(gr, out) {
  n_exons <- function(x) {
    tibble::as_tibble(x) |>
      dplyr::group_by(tx_id) |>
      dplyr::summarise(
        estimate = dplyr::first(estimate),
        n = dplyr::n(),
        .groups = "drop"
      )
  }
  n_exons(gr) |>
    dplyr::left_join(
      n_exons(out) |> dplyr::select(tx_id, n_out = n),
      by = "tx_id"
    ) |>
    dplyr::mutate(lost = n - dplyr::coalesce(n_out, 0L)) |>
    dplyr::filter(lost > 0L) |>
    dplyr::select(tx_id, estimate, lost)
}

test_that("generate_se, ie, cse and cie remove a run of internal exons", {
  # each removes n_skip consecutive internal exons from one transcript per
  # gene: an up transcript for se / cse, a down one for ie / cie. sim_event
  # marks the copies the other side keeps, so it sits on the opposite sign.
  # ranks follow transcription, so the run is the same on either strand
  for (st in c("+", "-")) {
    set.seed(5)
    gr <- create_mock_data(
      n_genes = 3, n_tx_per_gene = 3, n_exons_per_tx = 7, strand = st
    )
    cases <- list(
      se = list(out = generate_se(gr, 3), n_skip = 1L, edited = 1),
      ie = list(out = generate_ie(gr, 3), n_skip = 1L, edited = -1),
      cse = list(out = generate_cse(gr, 3), n_skip = 2L, edited = 1),
      cie = list(
        out = generate_cie(gr, 3, n_skip = 3), n_skip = 3L, edited = -1
      )
    )
    for (ev in names(cases)) {
      case <- cases[[ev]]
      info <- paste(ev, st)

      # one transcript per gene, on the edited side, n_skip exons short
      lost <- exons_lost(gr, case$out)
      expect_equal(nrow(lost), 3L, info = info)
      expect_true(all(sign(lost$estimate) == case$edited), info = info)
      expect_true(all(lost$lost == case$n_skip), info = info)

      # the run as the other side keeps it: n_skip consecutive internal
      # exons per transcript
      runs <- tibble::as_tibble(case$out) |>
        dplyr::filter(sim_event) |>
        dplyr::group_by(tx_id) |>
        dplyr::summarise(
          side = unique(sign(estimate)),
          n = dplyr::n(),
          consecutive = all(diff(sort(exon_rank)) == 1L),
          all_internal = all(internal),
          .groups = "drop"
        )
      expect_gt(nrow(runs), 0L)
      expect_true(all(runs$side == -case$edited), info = info)
      expect_true(all(runs$n == case$n_skip), info = info)
      expect_true(all(runs$consecutive & runs$all_internal), info = info)
    }
  }
})

test_that("generate_cse and generate_cie check n_skip and cap n_events", {
  gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 2, n_exons_per_tx = 4)
  # a run of one is generate_se() / generate_ie()
  expect_error(generate_cse(gr, n_skip = 1), "n_skip")
  expect_error(generate_cie(gr, n_skip = 2.5), "n_skip")
  expect_error(generate_cse(gr, n_skip = c(2, 3)), "n_skip")

  # four exons leave two internal ones: room for a run of two, not three
  expect_silent(generate_cse(gr, n_events = 2))
  expect_warning(
    out <- generate_cie(gr, n_events = 1, n_skip = 3), "capping at 0"
  )
  expect_equal(length(out), length(gr))
  expect_false(any(out$sim_event))
})

test_that("generate_ri and generate_si fuse exons on opposite sides", {
  # both fuse exons 2 and 3 of a transcript: an up one for ri, a down one
  # for si. sim_event is the fused exon itself, on the edited side
  for (st in c("+", "-")) {
    set.seed(9)
    gr <- create_mock_data(
      n_genes = 3, n_tx_per_gene = 2, n_exons_per_tx = 5, strand = st
    )
    cases <- list(
      ri = list(out = generate_ri(gr, n_events = 2), edited = 1),
      si = list(out = generate_si(gr, n_events = 2), edited = -1)
    )
    for (ev in names(cases)) {
      case <- cases[[ev]]
      info <- paste(ev, st)
      # merging toward the wrong neighbour on "-" would invert the range
      expect_true(all(GenomicRanges::width(case$out) > 0), info = info)

      lost <- exons_lost(gr, case$out)
      expect_equal(nrow(lost), 2L, info = info)
      expect_true(all(sign(lost$estimate) == case$edited), info = info)
      expect_true(all(lost$lost == 1L), info = info)

      # each fused exon spans two exons of the transcript it came from
      fused <- case$out |> dplyr::filter(sim_event)
      expect_equal(length(fused), 2L, info = info)
      expect_true(all(sign(fused$estimate) == case$edited), info = info)
      # (the two edited transcripts are in different genes, which do not
      # overlap, so each fused exon only meets its own transcript here)
      n_spanned <- plyranges::count_overlaps(
        fused, gr |> dplyr::filter(tx_id %in% fused$tx_id)
      )
      expect_true(all(n_spanned == 2L), info = info)
    }
  }
})
