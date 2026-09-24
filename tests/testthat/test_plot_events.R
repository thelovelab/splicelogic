# Tests for plot_events.R

# one "+" gene, three transcripts with exact coordinates:
# - down:  1-100, 201-300, 401-500, 601-700
# - up_se: skips 201-300, so 201-300 of down is an se
# - up_ss: 221-280 in place of 201-300, both boundaries moved, so an a5ss
#   and an a3ss on the same exon
# seqnames takes one value, or one per exon to split a transcript across
# chromosomes; tx_suffix keeps transcript ids apart when two copies are joined
plot_fixture <- function(gene_id = "g1", seqnames = "chr1", tx_suffix = "") {
  data.frame(
    seqnames = seqnames,
    start = c(1, 201, 401, 601, 1, 401, 601, 1, 221, 401, 601),
    end = c(100, 300, 500, 700, 100, 500, 700, 100, 280, 500, 700),
    strand = "+",
    gene_id = gene_id,
    tx_id = paste0(
      c(rep("down", 4), rep("up_se", 3), rep("up_ss", 4)), tx_suffix
    ),
    exon_rank = c(1:4, 1:3, 1:4),
    estimate = c(rep(-1, 4), rep(2, 3), rep(1, 4))
  ) |>
    plyranges::as_granges() |>
    preprocess(coef_col = "estimate")
}

# plot_transcripts() layers, in drawing order
layer_tx_lines <- 1
layer_exons <- 2
layer_event_exons <- 3
layer_tx_labels <- 4
layer_event_labels <- 5
# the strand arrows are always the last layer
strand_arrows <- function(p) {
  ggplot2::layer_data(p, length(p$layers))
}

test_that("event_palette covers every event type the finders return", {
  expect_setequal(
    names(event_palette()),
    c("se", "ie", "mxe", "ri", "a5ss", "a3ss", "atss", "ates")
  )
  # the up / down transcript colours are not reused for an event type
  direction_colors <- unlist(
    formals(plot_transcripts)[c("up_color", "down_color")]
  )
  expect_false(any(event_palette() %in% direction_colors))
})

test_that("plot_transcripts draws one row per transcript, by estimate", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)
  p <- plot_transcripts(events, gr)

  expect_s3_class(p, "ggplot")
  tx_lines <- ggplot2::layer_data(p, layer_tx_lines)
  expect_equal(nrow(tx_lines), 3L)

  # most up-regulated transcript at the top
  labels <- ggplot2::layer_data(p, layer_tx_labels)
  labels <- labels[order(labels$y, decreasing = TRUE), ]
  expect_equal(
    labels$label,
    c("up_se  (+2.00)", "up_ss  (+1.00)", "down  (-1.00)")
  )
})

test_that("plot_transcripts draws a strand arrow over each first exon", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_se(gr)
  p <- plot_transcripts(events, gr)

  # three points per transcript, drawn in black: the transcription start on
  # the top edge of the first exon, up half an exon box (boxes are 0.44
  # tall), then 1.1% of the panel across
  arrows <- strand_arrows(p)
  expect_equal(nrow(arrows), 6L)
  expect_true(all(arrows$colour == "black"))

  up_se <- arrows[arrows$group == min(arrows$group), ]
  row_y <- up_se$y[1] - 0.22
  expect_equal(up_se$y - row_y, c(0.22, 0.44, 0.44))
  # first exon 1-100 is drawn from 50, after the 50 bp flank, and the panel
  # runs 0-649 (last exon ends at 599, plus the flank); "+" points right
  expect_equal(up_se$x, c(50, 50, 50 + 0.011 * 649))
})

test_that("plot_transcripts highlights the event exons only", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)
  p <- plot_transcripts(events, gr)

  # 11 exons in the panel, two of them events
  expect_equal(nrow(ggplot2::layer_data(p, layer_exons)), 9L)
  event_rects <- ggplot2::layer_data(p, layer_event_exons)
  expect_equal(nrow(event_rects), 2L)

  # the a5ss and a3ss on 221-280 collapse into one exon and one label,
  # filled by the first type
  event_labels <- ggplot2::layer_data(p, layer_event_labels)
  expect_setequal(event_labels$label, c("SE", "A3SS+A5SS"))
  expect_setequal(
    event_rects$fill,
    unname(event_palette()[c("se", "a3ss")])
  )
})

test_that("plot_transcripts colours transcripts by the sign of estimate", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_se(gr)
  p <- plot_transcripts(events, gr, up_color = "blue", down_color = "red")

  tx_lines <- ggplot2::layer_data(p, layer_tx_lines)
  tx_lines <- tx_lines[order(tx_lines$y, decreasing = TRUE), ]
  # find_se() pairs down with up_se only
  expect_equal(tx_lines$colour, c("blue", "red"))
})

test_that("plot_transcripts rescales introns and keeps exon widths", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)
  p <- plot_transcripts(
    events, gr, new_intron_length = 50, flanking_length = c(50, 50)
  )
  rects <- rbind(
    ggplot2::layer_data(p, layer_exons)[, c("xmin", "xmax")],
    ggplot2::layer_data(p, layer_event_exons)[, c("xmin", "xmax")]
  )

  # union of all exons: 1-100, 201-300, 401-500, 601-700. the first starts
  # after the flank and every gap between them becomes 50 bp
  expect_setequal(unique(rects$xmin), c(50, 200, 220, 350, 500))
  # 221-280 sits 20 bp into the 201-300 union range
  expect_true(all((rects$xmax - rects$xmin) %in% c(99, 59)))

  # without rescaling the x axis is genomic
  p_genomic <- plot_transcripts(events, gr, rescale_introns = FALSE)
  event_rects <- ggplot2::layer_data(p_genomic, layer_event_exons)
  expect_setequal(event_rects$xmin, c(201, 221))
  expect_setequal(event_rects$xmax, c(300, 280))
})

test_that("plot_transcripts appends tx_annot to the transcript label", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_se(gr)
  p <- plot_transcripts(events, gr, tx_annot = c(up_se = "novel"))

  labels <- ggplot2::layer_data(p, layer_tx_labels)$label
  expect_setequal(labels, c("up_se  (+2.00)  novel", "down  (-1.00)"))
})

test_that("plot_transcripts works with numeric tx_id", {
  skip_if_not_installed("ggplot2")
  # se_mock_data() has tx_id 1, 2, 3 stored as numbers
  gr <- preprocess(se_mock_data(), coef_col = "estimate")
  events <- find_se(gr)
  p <- expect_no_warning(plot_transcripts(events, gr))

  expect_equal(nrow(ggplot2::layer_data(p, layer_tx_lines)), 3L)
  # exons 21-25 and 51-55 of tx 1
  expect_equal(nrow(ggplot2::layer_data(p, layer_event_exons)), 2L)
})

test_that("plot_transcripts highlights every event on minus-strand data", {
  skip_if_not_installed("ggplot2")
  gr <- preprocess(minus_strand_mock_data(), coef_col = "estimate")
  events <- find_all_events(gr, verbose = FALSE)

  plots <- expect_no_warning(plot_transcripts_by_gene(events, gr))
  expect_named(plots, paste0("g", 1:6), ignore.order = TRUE)

  # one event exon per gene, except the mxe pair in g2
  n_event_exons <- vapply(
    plots,
    function(p) nrow(ggplot2::layer_data(p, layer_event_exons)),
    integer(1)
  )
  expect_equal(
    n_event_exons[paste0("g", 1:6)],
    c(g1 = 1L, g2 = 2L, g3 = 1L, g4 = 1L, g5 = 1L, g6 = 1L)
  )

  # transcription runs right to left: every arrow starts on the right edge of
  # the rightmost exon and points left
  p <- plots[["g1"]]
  arrows <- strand_arrows(p)
  exon_rects <- rbind(
    ggplot2::layer_data(p, layer_exons)[, c("xmin", "xmax")],
    ggplot2::layer_data(p, layer_event_exons)[, c("xmin", "xmax")]
  )
  for (tx_arrow in split(arrows, arrows$group)) {
    expect_equal(tx_arrow$x[1], max(exon_rects$xmax))
    expect_lt(tx_arrow$x[3], tx_arrow$x[1])
  }
})

test_that("plot_transcripts warns about what it cannot draw", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)

  # a transcript missing from exons is dropped
  expect_warning(
    plot_transcripts(events, gr |> dplyr::filter(tx_id != "up_se")),
    "Dropping transcript\\(s\\) absent from 'exons': up_se"
  )

  # an event exon whose coordinates are not in exons is not highlighted
  moved <- GenomicRanges::shift(events[1], 5)
  expect_warning(
    plot_transcripts(moved, gr),
    "1 event exon\\(s\\) had no match"
  )
})

test_that("plot_transcripts refuses events from more than one gene", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()

  # a second gene on the same chromosome: nothing else would stop it, and the
  # gap between the genes would be collapsed like an intron
  same_chr <- plot_fixture(gene_id = "g2", tx_suffix = "_g2")
  expect_error(
    plot_transcripts(
      c(find_se(gr), find_se(same_chr)),
      c(gr, same_chr)
    ),
    "'events' spans 2 genes \\(g1, g2\\)"
  )

  # a second gene on another chromosome fails on the gene too, not later on
  # the chromosome
  other_chr <- plot_fixture(
    gene_id = "g2", seqnames = "chr2", tx_suffix = "_g2"
  )
  # c() warns that the two share no sequence levels, which is the point
  events_2chr <- suppressWarnings(c(find_se(gr), find_se(other_chr)))
  exons_2chr <- suppressWarnings(c(gr, other_chr))
  expect_error(
    plot_transcripts(events_2chr, exons_2chr),
    "'events' spans 2 genes"
  )

  # without gene_id on the events, the exons still give the genes away
  expect_error(
    plot_transcripts(
      c(find_se(gr), find_se(same_chr)) |>
        dplyr::select(-gene_id),
      c(gr, same_chr)
    ),
    "'events' spans 2 genes"
  )
})

test_that("plot_transcripts rejects input it cannot use", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)

  # a finder with no hits returns a bare GRanges()
  expect_error(plot_transcripts(find_ri(gr), gr), "'events' has no rows")
  expect_error(
    plot_transcripts(events |> dplyr::select(-event_estimate), gr),
    "Missing metadata columns in 'events': event_estimate"
  )
  expect_error(
    plot_transcripts(events, S4Vectors::split(gr, gr$tx_id)),
    "'exons' must be a flat GRanges"
  )
  expect_error(
    suppressWarnings(
      plot_transcripts(events, gr |> dplyr::filter(tx_id == "none"))
    ),
    "None of the event transcripts are present"
  )

  # one gene, but up_se's exons are on another chromosome
  split_chr <- plot_fixture(
    seqnames = c(rep("chr1", 4), rep("chr2", 3), rep("chr1", 4))
  )
  expect_error(
    plot_transcripts(find_se(gr), split_chr),
    "more than one chromosome"
  )
})

test_that("plot_transcripts_by_gene returns one plot per gene", {
  skip_if_not_installed("ggplot2")
  gr <- preprocess(minus_strand_mock_data(), coef_col = "estimate")
  events <- find_all_events(gr, verbose = FALSE)

  plots <- plot_transcripts_by_gene(events, gr, genes = c("g1", "g3"))
  expect_named(plots, c("g1", "g3"))
  expect_equal(plots[["g1"]]$labels$title, "g1")

  # title is a prefix, and does not clash with the per-gene title
  plots <- plot_transcripts_by_gene(events, gr, genes = "g1", title = "KD")
  expect_equal(plots[["g1"]]$labels$title, "KD: g1")

  expect_error(
    plot_transcripts_by_gene(events, gr, genes = "none"),
    "No events left to plot"
  )
  expect_error(
    plot_transcripts_by_gene(events |> dplyr::select(-gene_id), gr),
    "Missing metadata columns in 'events': gene_id"
  )
})

test_that("plot_transcripts_by_gene skips a gene it cannot draw", {
  skip_if_not_installed("ggplot2")
  gr <- preprocess(minus_strand_mock_data(), coef_col = "estimate")
  events <- find_all_events(gr, verbose = FALSE)

  # no exons for g1's transcripts: plot_transcripts() warns that it drops
  # them, then fails, and the wrapper turns that into a skip
  exons <- gr |>
    dplyr::filter(gene_id != "g1")
  expect_warning(
    expect_warning(plot_transcripts_by_gene(events, exons), "Skipping g1"),
    "Dropping transcript"
  )
  plots <- suppressWarnings(plot_transcripts_by_gene(events, exons))
  expect_named(plots, paste0("g", 2:6), ignore.order = TRUE)
})

test_that("plot_transcripts_by_pair returns one plot per unordered pair", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)

  # the se is reported from down, the a5ss / a3ss from up_ss
  plots <- plot_transcripts_by_pair(events, gr)
  expect_named(plots, c("down-up_se", "down-up_ss"))
  expect_equal(
    nrow(ggplot2::layer_data(plots[["down-up_ss"]], layer_tx_lines)),
    2L
  )

  # an existing tx_pair column is used as is
  grouped <- events |>
    dplyr::mutate(tx_pair = "all")
  expect_named(plot_transcripts_by_pair(grouped, gr), "all")
})

test_that("plot_transcript_lfc draws one bar per transcript with a value", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)
  tx_lfc <- c(down = -1.5, up_se = 3, other = 10)

  p <- plot_transcript_lfc(events, tx_lfc)
  expect_s3_class(p, "ggplot")

  # up_ss has no value, and "other" is not in the events
  bars <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(bars), 2L)
  expect_setequal(bars$xend, c(-1.5, 3))

  # values beyond lfc_max are clipped at the axis edge
  clipped <- plot_transcript_lfc(events, tx_lfc, lfc_max = 2) |>
    ggplot2::layer_data(2)
  expect_setequal(clipped$xend, c(-1.5, 2))
})

test_that("plot_transcript_lfc lines up with plot_transcripts", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)
  tx_lfc <- c(down = -1.5, up_se = 3, up_ss = 0.5)

  p_tx <- plot_transcripts(events, gr)
  p_lfc <- plot_transcript_lfc(events, tx_lfc, transcript_label = TRUE)

  y_range <- function(p) {
    ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
  }
  expect_equal(y_range(p_lfc), y_range(p_tx))

  # same row for the same transcript
  tx_rows <- ggplot2::layer_data(p_tx, layer_tx_lines)
  bars <- ggplot2::layer_data(p_lfc, 2)
  expect_setequal(bars$y, tx_rows$y)
  y_axis <- ggplot2::ggplot_build(p_lfc)$layout$panel_params[[1]]$y
  expect_setequal(y_axis$get_labels(), c("down", "up_se", "up_ss"))
})

test_that("plot_transcript_lfc checks tx_lfc", {
  skip_if_not_installed("ggplot2")
  gr <- plot_fixture()
  events <- find_all_events(gr, verbose = FALSE)

  expect_error(
    plot_transcript_lfc(events, c(1, 2)),
    "'tx_lfc' must be a numeric vector named by transcript id"
  )
  expect_warning(
    plot_transcript_lfc(events, c(other = 1)),
    "No transcript in 'events' has a value in 'tx_lfc'"
  )

  p <- plot_transcript_lfc(events, c(down = -1), lfc_ends = c("KD", "WT"))
  end_labels <- ggplot2::layer_data(p, 3)$label
  expect_equal(end_labels, c("KD", "WT"))
  expect_error(plot_transcript_lfc(events, c(down = -1), lfc_ends = "KD"))
})
