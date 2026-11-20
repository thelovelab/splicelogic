#' Create a sample GRanges with one negative coef
#' transcript and two positive coef transcripts
#' This dataset is designed to include a skipped exon
#' event at exon_rank 3 and 5 of tx_id 1
#' happening between both tx_id 2 and tx_id 3
#' @return A GRanges object with two transcripts per gene
#' @importFrom stats runif
#' @noRd
se_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41, 51, 61),
    width = 5,
    strand = "+",
    exon_rank = seq_len(7),
    gene_id = rep(1, 7),
    estimate = rep(runif(1, min = -1, max = 0), 7)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = seq_len(6),
    gene_id = rep(1, 6),
    estimate = rep(runif(1, min = 0, max = 1), 6)
  )
  df3 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = seq_len(6),
    gene_id = rep(1, 6),
    estimate = rep(runif(1, min = 0, max = 1), 6)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr3 <- plyranges::as_granges(df3)
  gr <- plyranges::bind_ranges(gr1, gr2, gr3) |>
    dplyr::mutate(
      tx_id = c(rep(1, 7), rep(2, 6), rep(3, 6))
    )

  return(gr)
}

#' Create a sample GRanges with one negative coef
#' transcript and one positive coef transcripts
#' This dataset is designed to include mutually exclusive
#' exons between exon_rank 3 of tx_id 1 and exon_rank 3
#' of tx_id 2. There is also a skipped exon event at
#' exon_rank 5 of tx_id 1.
#' @return A GRanges object with two transcripts per gene and candidate logic
#' @noRd
mx_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 41, 51, 61, 71, 91, 101),
    width = 5,
    strand = "+",
    exon_rank = seq_len(9),
    gene_id = rep(1, 9),
    estimate = rep(runif(1, min = -1, max = 0), 9)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 41, 61, 71, 81, 101),
    width = 5,
    strand = "+",
    exon_rank = seq_len(8),
    gene_id = rep(1, 8),
    estimate = rep(runif(1, min = 0, max = 1), 8)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
    dplyr::mutate(
      tx_id = c(rep(1, 9), rep(2, 8))
    )

  return(gr)
}
#' Create a sample GRanges with one negative coef
#' transcript and one positive coef transcripts
#' This dataset is designed to include no splicing events
#' @return A GRanges object with two transcripts per gene and no candidate logic
#' @noRd
no_event_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21),
    width = 5,
    strand = "+",
    exon_rank = seq_len(3),
    gene_id = rep(1, 3),
    estimate = rep(runif(1, min = -1, max = 0), 3)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21),
    width = 5,
    strand = "+",
    exon_rank = seq_len(3),
    gene_id = rep(1, 3),
    estimate = rep(runif(1, min = 0, max = 1), 3)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
    dplyr::mutate(
      tx_id = c(rep(1, 3), rep(2, 3))
    )

  return(gr)
}

#' Create a sample GRanges whose exon_rank does not start at 1
#'
#' Mimics coding-sequence (CDS) ranges rather than whole exons: a 5'
#' UTR-only exon is not part of the CDS, so the first coding range of
#' each transcript carries exon_rank 2 instead of 1. Ranks are still
#' consecutive within a transcript, they just start at an offset.
#'
#' ```
#' rank     2        3        4        5        6
#' neg1  [1-5]  [11-15]  [21-25]  [31-35]  [41-45]
#' pos2  [1-5]  [11-15]  [21-25]           [41-45]
#' rank     2        3        4                 5
#' ```
#'
#' tx 2 skips 31-35, so there is one skipped exon event at exon_rank 5
#' of tx 1. Because the ranks are offset, this fixture also covers the
#' 'internal' flag being computed relative to each transcript's own
#' rank range rather than assuming ranks run 1..nexons.
#' @return A GRanges object with two transcripts in one gene
#' @noRd
cds_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41),
    width = 5,
    strand = "+",
    exon_rank = 2:6,
    gene_id = rep(1, 5),
    estimate = rep(runif(1, min = -1, max = 0), 5)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 41),
    width = 5,
    strand = "+",
    exon_rank = 2:5,
    gene_id = rep(1, 4),
    estimate = rep(runif(1, min = 0, max = 1), 4)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
    dplyr::mutate(
      tx_id = c(rep(1, 5), rep(2, 4))
    )

  return(gr)
}

#' Create a sample GRanges with one negative coef transcript and three
#' positive coef transcripts, two of which splice an intron inside the
#' candidate exon
#'
#' ```
#' neg1  [1-50] [====100-300====] [400-450]
#' posA  [1-50] [100-150] [200-300] [400-450]
#' posB  [1-50] [100-150] [200-300] [400-450]
#' posC  [1-50]                     [400-450]
#' ```
#'
#' posC splices `1-50` straight to `400-450`, so 100-300 is a genuine
#' skipped exon relative to posC. posA and posB each contribute TWO
#' exons overlapping the candidate, so counting exons rather than
#' isoforms gives 4 against 3 positive transcripts and discards the
#' candidate before it reaches the flanking-exon stage.
#' @return A GRanges object with four transcripts in one gene
#' @noRd
intron_split_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 100, 400),
    end = c(50, 300, 450),
    strand = "+",
    exon_rank = seq_len(3),
    gene_id = "g1",
    estimate = -1
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 100, 200, 400),
    end = c(50, 150, 300, 450),
    strand = "+",
    exon_rank = seq_len(4),
    gene_id = "g1",
    estimate = 1
  )
  df3 <- df2
  df4 <- data.frame(
    seqnames = "chr1",
    start = c(1, 400),
    end = c(50, 450),
    strand = "+",
    exon_rank = seq_len(2),
    gene_id = "g1",
    estimate = 1
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr3 <- plyranges::as_granges(df3)
  gr4 <- plyranges::as_granges(df4)
  gr <- plyranges::bind_ranges(gr1, gr2, gr3, gr4) |>
    dplyr::mutate(
      tx_id = rep(c("neg1", "posA", "posB", "posC"), times = c(3, 4, 4, 2))
    )

  return(gr)
}

#' Create a sample GRanges on the minus strand covering every event type
#'
#' Seven one-gene loci, all on `-`, each engineered for a single event. What
#' distinguishes this fixture from the `+` strand ones is that exon_rank
#' runs right to left in genomic coordinates (rank 1 is the rightmost
#' exon), so the exon of rank k - 1 sits to the genomic *right* of rank k,
#' and an exon's donor / acceptor are its genomic start / end rather than
#' the other way round.
#'
#' ```
#' g1 se   neg1 [1-5]     [11-15]   [21-25]             [31-35]   [41-45]
#'         pos1 [1-5]     [11-15]                       [31-35]   [41-45]
#' g2 mxe  neg2 [100-105] [111-115] [121-125]           [141-145] [151-155]
#'         pos2 [100-105] [111-115]           [131-135] [141-145] [151-155]
#' g3 ri   neg3 [200-205] [211-215] [221-225]
#'         pos3 [200-205] [===211-225===]
#' g4 a5ss neg4 [300-305] [311-315] [321-325]
#'         pos4 [300-305]   [313-315] [321-325]
#' g5 a3ss neg5 [400-405] [411-415] [421-425]
#'         pos5 [400-405] [411-413] [421-425]
#' g6 ie   neg6 [500-505] [511-515]           [531-535] [541-545]
#'         pos6 [500-505] [511-515] [521-525] [531-535] [541-545]
#' g7 si   neg7 [600-605] [===611-625===]
#'         pos7 [600-605] [611-615] [621-625]
#' ```
#'
#' In g4 the moved boundary is the exon start, which on `-` is the donor,
#' so it is an a5ss; in g5 it is the exon end, the acceptor, so it is an
#' a3ss. Both are the mirror image of the `+` strand case. g6 and g7 are
#' g1 and g3 with the up and down roles swapped.
#' @return A GRanges object with fourteen transcripts in seven genes
#' @noRd
minus_strand_mock_data <- function() {
  tx <- function(gene_id, tx_id, start, end, estimate) {
    data.frame(
      seqnames = "chr1",
      start = start,
      end = end,
      strand = "-",
      # rank 1 is the 5'-most exon, which on "-" is the rightmost one
      exon_rank = rev(seq_along(start)),
      gene_id = gene_id,
      tx_id = tx_id,
      estimate = estimate
    )
  }

  dplyr::bind_rows(
    tx("g1", "neg1", c(1, 11, 21, 31, 41), c(5, 15, 25, 35, 45), -1),
    tx("g1", "pos1", c(1, 11, 31, 41), c(5, 15, 35, 45), 1),
    tx(
      "g2", "neg2",
      c(100, 111, 121, 141, 151), c(105, 115, 125, 145, 155), -1
    ),
    tx(
      "g2", "pos2",
      c(100, 111, 131, 141, 151), c(105, 115, 135, 145, 155), 1
    ),
    tx("g3", "neg3", c(200, 211, 221), c(205, 215, 225), -1),
    tx("g3", "pos3", c(200, 211), c(205, 225), 1),
    tx("g4", "neg4", c(300, 311, 321), c(305, 315, 325), -1),
    tx("g4", "pos4", c(300, 313, 321), c(305, 315, 325), 1),
    tx("g5", "neg5", c(400, 411, 421), c(405, 415, 425), -1),
    tx("g5", "pos5", c(400, 411, 421), c(405, 413, 425), 1),
    tx("g6", "neg6", c(500, 511, 531, 541), c(505, 515, 535, 545), -1),
    tx(
      "g6", "pos6",
      c(500, 511, 521, 531, 541), c(505, 515, 525, 535, 545), 1
    ),
    tx("g7", "neg7", c(600, 611), c(605, 625), -1),
    tx("g7", "pos7", c(600, 611, 621), c(605, 615, 625), 1)
  ) |>
    plyranges::as_granges()
}

# for create_mock_data
utils::globalVariables(c("gene_id", "strand", "tx_order", "estimate"))

#' Create mock GRanges data for splicing event testing
#' 
#' @param n_genes Number of genes to simulate
#' @param n_tx_per_gene Number of transcripts per gene. Use 2 or more:
#'   the first transcript of each gene is given a negative estimate and
#'   the second a positive one, so a single transcript per gene leaves
#'   the set all-negative and the `generate_*()` helpers with nothing to
#'   modify.
#' @param n_exons_per_tx Number of exons per transcript
#' @param coef_range Range of coefficient values to sample from
#' @param strand Strand to place every transcript on: `"+"` (default),
#'   `"-"`, or `"*"`. On `"-"` the exon ranks are reversed, so exon_rank 1
#'   is the rightmost exon in genomic coordinates.
#'
#' @return A GRanges object with simulated transcripts and exons
#'
#' @examples
#'
#' # create mock data with 2 genes, 4 transcripts
#' # per gene, and 4 exons per transcript
#' gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4)
#'
#' # the same, on the minus strand
#' gr_minus <- create_mock_data(n_genes = 2, n_exons_per_tx = 4, strand = "-")
#'
#' @export
create_mock_data <- function(
  n_genes = 1,
  n_tx_per_gene = 2,
  n_exons_per_tx = 5,
  coef_range = c(-1, 1),
  strand = c("+", "-", "*")
) {
  # bound to a distinct name so the mutate() below reads the argument
  # rather than the column it is creating
  tx_strand <- match.arg(strand)

  # the neg/pos pairing below is per gene and keys off tx_order 1 and 2,
  # so it cannot be honoured with a single transcript per gene
  if (n_tx_per_gene < 2) {
    warning(
      "n_tx_per_gene < 2: every transcript is assigned a negative ",
      "estimate, so there is no\n  positive isoform to compare against ",
      "and the generate_*() helpers have nothing to modify.",
      call. = FALSE
    )
  }

  # Generate all combinations of genes, transcripts, and exons
  data <- expand.grid(
    gene_id = seq_len(n_genes),
    tx_id = seq_len(n_tx_per_gene),
    exon_rank = seq_len(n_exons_per_tx)
  )

  data <- data |>
    dplyr::arrange(gene_id, tx_id, exon_rank)

  # gene_offset scales with n_exons_per_tx to prevent cross-gene coordinate overlap
  gene_offset <- n_exons_per_tx * 10 + 50

  # Calculate transcript IDs globally
  data <- data |>
    dplyr::mutate(
      tx_id = tx_id + (gene_id - 1) * n_tx_per_gene,
      seqnames = paste0("chr", sample(seq_len(22), 1)), # Random chromosome
      start = (exon_rank - 1) * 10 + 1 + (gene_id - 1) * gene_offset,
      width = 5, # Fixed width
      strand = tx_strand # one strand for the whole simulated set
    )

  # Reverse exon_rank for transcripts with strand == "-"
  data <- data |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      exon_rank = ifelse(strand == "-", rev(exon_rank), exon_rank)
    ) |>
    dplyr::ungroup()

  # Assign estimate per tx_id: 1 neg and 1 pos per gene when
  # n_tx_per_gene >= 2 (warned about above when it is not)
  # Create a lookup table for tx estimate:
  # first tx gets negative, second gets positive, rest random
  tx_estimate <- data |>
    dplyr::distinct(gene_id, tx_id) |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(tx_order = dplyr::row_number()) |>
    dplyr::ungroup()

  n_tx <- nrow(tx_estimate)
  # Generate all random estimate first, then override first two per gene
  tx_estimate <- tx_estimate |>
    dplyr::mutate(
      estimate = runif(n_tx, min = coef_range[1], max = coef_range[2]),
      estimate = dplyr::if_else(
        tx_order == 1,
        runif(n_tx, min = coef_range[1], max = -0.01),
        estimate
      ),
      estimate = dplyr::if_else(
        tx_order == 2,
        runif(n_tx, min = 0.01, max = coef_range[2]),
        estimate
      )
    ) |>
    dplyr::select(tx_id, estimate)

  # Join estimate back to data
  data <- data |>
    dplyr::left_join(tx_estimate, by = "tx_id")

  # Convert to GRanges
  gr <- plyranges::as_granges(data)
  gr <- preprocess(gr, coef_col = "estimate")

  return(gr)
}

#' @rdname generate_events
#' @name generate_events
#'
#' @title Generate mock splice events in a GRanges object
#'
#' @description Functions to introduce specific types of alternative splicing
#' events into mock GRanges data for testing purposes.
#'
#' Each generator edits the transcripts on one side and marks, in a
#' logical `sim_event` column, the exons the matching finder should
#' report. `generate_se()`, `generate_cse()`, `generate_ri()`,
#' `generate_a5ss()`, `generate_a3ss()`, `generate_atss()` and
#' `generate_ates()` edit up (positive `estimate`) transcripts.
#' `generate_ie()`, `generate_cie()` and `generate_si()` are their mirrors:
#' the same edit made to a down (negative `estimate`) transcript.
#' `generate_mxe()` edits one transcript on each side.
#'
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id',
#' 'tx_id', and 'estimate'.
NULL

# for sample_per_gene
utils::globalVariables(c("gene_id"))

#' Draw one candidate per gene, then up to n_events of those genes
#'
#' The sampling every generator uses, so that `n_events` means the same
#' thing for all of them: the number of genes that get an event. Two events
#' in one gene would edit the same transcripts and could break each other
#' (or one could hide the other), so a gene never gets more than one.
#'
#' @param candidates A data frame with a `gene_id` column and one row per
#'   candidate edit
#' @param n_events Number of events to generate
#' @return At most `n_events` rows of `candidates`, from distinct genes.
#'   Warns, and returns one row per gene, when fewer genes than `n_events`
#'   have a candidate.
#' @noRd
sample_per_gene <- function(candidates, n_events) {
  per_gene <- candidates |>
    dplyr::group_by(gene_id) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup()

  n_available <- nrow(per_gene)
  if (n_events > n_available) {
    warning(
      sprintf(
        paste0(
          "n_events (%d) exceeds available genes with candidates (%d); ",
          "capping at %d."
        ),
        n_events, n_available, n_available
      ),
      call. = FALSE
    )
    n_events <- n_available
  }

  per_gene |>
    dplyr::slice_sample(n = n_events)
}

# for generate_skip_block
utils::globalVariables(c(
  "estimate", "internal", "key", "gene_id", "sim_event", "coord_sig",
  "tx_id", "exon_rank", "last_rank", "run_step"
))

#' Remove a run of consecutive internal exons from one transcript per gene
#'
#' Shared implementation for generate_se() / generate_ie() /
#' generate_cse() / generate_cie(). One run start is sampled per gene,
#' up to `n_events` genes are drawn from that pool, and the `n_skip` exons
#' from each start onwards (by rank) are removed from that transcript. The
#' other transcripts of the gene keep those exons, so on the ones of the
#' opposite sign they are the event, and that is where `sim_event` is set.
#'
#' Without `inverse` an up transcript loses the run, so a down transcript's
#' copy is skipped (se / cse); with `inverse` a down transcript loses it,
#' so an up transcript's copy is included (ie / cie). Every exon of the run
#' must be internal, so that both flanks exist. Ranks follow transcription,
#' so the run is the same on either strand.
#'
#' @param gr A GRanges object with an 'estimate' column
#' @param n_events Number of events to generate, at most one per gene
#' @param n_skip Number of consecutive exons removed
#' @param inverse If TRUE, removes the run from a down transcript
#' @return A GRanges object with the events introduced
#' @noRd
generate_skip_block <- function(gr, n_events = 1, n_skip = 1L,
                                inverse = FALSE) {
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }
  # 1: the run is removed from an up transcript; -1: from a down one
  factor <- if (inverse) -1 else 1

  # candidate run starts: internal exons on the edited side whose whole run
  # stays clear of the last exon. with n_skip = 1 the run filter is implied
  # by `internal`, so this is every internal exon.
  candidates <- tibble::as_tibble(gr) |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(last_rank = max(exon_rank)) |>
    dplyr::ungroup() |>
    dplyr::filter(
      estimate * factor > 0 & internal == TRUE &
        exon_rank + n_skip - 1L < last_rank
    ) |>
    dplyr::distinct(gene_id, key, tx_id, exon_rank)

  # keys of the whole run, from each sampled start onwards by rank
  run_keys <- candidates |>
    sample_per_gene(n_events) |>
    dplyr::cross_join(tibble::tibble(run_step = seq_len(n_skip) - 1L)) |>
    dplyr::mutate(key = paste0(tx_id, "-", exon_rank + run_step)) |>
    dplyr::pull(key)

  event_coords <- gr |>
    as.data.frame() |>
    dplyr::filter(key %in% run_keys) |>
    dplyr::mutate(coord_sig = paste(start, end, gene_id)) |>
    dplyr::pull(coord_sig)

  gr <- gr |>
    dplyr::filter(!key %in% run_keys)
  gr <- rerank_exons(gr)

  # the event is on the opposite side, whose transcripts still have the run
  gr <- gr |>
    dplyr::mutate(
      sim_event = estimate * factor < 0 &
        paste(start, end, gene_id) %in% event_coords
    )
  return(gr)
}

#' Stop unless n_skip is a single whole number of at least 2
#' @param n_skip Value passed to generate_cse() / generate_cie()
#' @return `n_skip` as an integer, invisibly
#' @noRd
check_n_skip <- function(n_skip) {
  if (
    !is.numeric(n_skip) || length(n_skip) != 1L || is.na(n_skip) ||
      n_skip < 2 || n_skip %% 1 != 0
  ) {
    stop(
      "n_skip must be a single whole number of at least 2; a run of one ",
      "exon is what generate_se() / generate_ie() make.",
      call. = FALSE
    )
  }
  invisible(as.integer(n_skip))
}

#' @rdname generate_events
#' @param n_events Number of events to generate, i.e. the number of genes
#'   that get one. Every generator puts at most one event in a gene, so two
#'   events never edit the same transcripts. When fewer genes than
#'   `n_events` have a candidate, it warns and makes one in each of those.
#' @return `generate_se()`: A GRanges object with skipped exon
#' events introduced: an internal exon is removed from one up transcript
#' per gene
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_se(gr, n_events = 1)
#'
generate_se <- function(gr, n_events = 1) {
  generate_skip_block(gr, n_events, n_skip = 1L)
}

#' @rdname generate_events
#' @return `generate_ie()`: A GRanges object with included exon events
#' introduced: an internal exon is removed from one down transcript per
#' gene, so the up transcripts include it
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_ie(gr, n_events = 1)
#'
generate_ie <- function(gr, n_events = 1) {
  generate_skip_block(gr, n_events, n_skip = 1L, inverse = TRUE)
}

#' @rdname generate_events
#' @param n_skip Number of consecutive exons skipped (or included) together,
#'   at least 2. Every exon of the run is internal, so a transcript needs
#'   at least `n_skip + 2` exons to host one.
#' @return `generate_cse()`: A GRanges object with consecutive skipped exon
#' events introduced: a run of `n_skip` internal exons is removed from one
#' up transcript per gene
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
#' )
#' generate_cse(gr, n_events = 1)
#'
#' # a run of three needs at least five exons
#' generate_cse(gr, n_events = 1, n_skip = 3)
#'
generate_cse <- function(gr, n_events = 1, n_skip = 2) {
  n_skip <- check_n_skip(n_skip)
  generate_skip_block(gr, n_events, n_skip)
}

#' @rdname generate_events
#' @return `generate_cie()`: A GRanges object with consecutive included exon
#' events introduced: a run of `n_skip` internal exons is removed from one
#' down transcript per gene, so the up transcripts include it
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
#' )
#' generate_cie(gr, n_events = 1)
#'
generate_cie <- function(gr, n_events = 1, n_skip = 2) {
  n_skip <- check_n_skip(n_skip)
  generate_skip_block(gr, n_events, n_skip, inverse = TRUE)
}

# for generate_mxe
utils::globalVariables(c(
  "estimate", "internal", "next_key", "gene_id", "neg_tx_id", "key",
  "coord_sig"
))

#' @rdname generate_events
#' @return `generate_mxe()`: A GRanges object with mutually exclusive exon
#' events introduced: a down transcript loses internal exon k + 1 and an up
#' transcript of the same gene loses exon k, so each keeps the exon the
#' other lacks. One pair per gene.
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_mxe(gr, n_events = 1)
#'
generate_mxe <- function(gr, n_events = 1) {
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }

  # Find consecutive pairs of internal exons from neg transcripts
  # Both rank k and rank k+1 must be internal
  neg_internal <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate < 0 & internal == TRUE)

  pos_txs <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate > 0) |>
    dplyr::distinct(gene_id, tx_id) |>
    dplyr::rename(pos_tx_id = tx_id)

  # candidate ranks k: rank k+1 is also internal in the same neg
  # transcript, and the gene has a pos transcript to pair with
  candidates <- neg_internal |>
    dplyr::mutate(
      next_key = paste0(tx_id, "-", exon_rank + 1L)
    ) |>
    dplyr::filter(next_key %in% neg_internal$key) |>
    dplyr::semi_join(pos_txs, by = "gene_id") |>
    dplyr::select(gene_id, tx_id, exon_rank)

  # Pair each selected neg transcript with one pos transcript from same gene
  mx_pairs <- candidates |>
    sample_per_gene(n_events) |>
    dplyr::rename(neg_tx_id = tx_id) |>
    dplyr::inner_join(pos_txs, by = "gene_id") |>
    dplyr::group_by(gene_id) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup()

  # Build keys to remove (pairwise):
  # rank k from the selected pos transcript (pos loses k, keeps k+1)
  # rank k+1 from the selected neg transcript (neg keeps k, loses k+1)
  pos_remove <- paste0(mx_pairs$pos_tx_id, "-", mx_pairs$exon_rank)
  neg_remove <- paste0(mx_pairs$neg_tx_id, "-", mx_pairs$exon_rank + 1L)

  # the event is the exon each side of the pair keeps: rank k of the neg
  # transcript and rank k+1 of the pos one. recorded by transcript and
  # coordinates rather than by key, because the re-rank below changes keys
  event_keys <- c(
    paste0(mx_pairs$neg_tx_id, "-", mx_pairs$exon_rank),
    paste0(mx_pairs$pos_tx_id, "-", mx_pairs$exon_rank + 1L)
  )
  event_coords <- gr |>
    as.data.frame() |>
    dplyr::filter(key %in% event_keys) |>
    dplyr::mutate(coord_sig = paste(tx_id, start, end)) |>
    dplyr::pull(coord_sig)

  gr <- gr |>
    dplyr::filter(!key %in% c(pos_remove, neg_remove))
  gr <- rerank_exons(gr)

  # only the edited pair is marked. other isoforms of the gene that still
  # have both exons are an se / ie against the pair, not an mxe.
  gr <- gr |>
    dplyr::mutate(sim_event = paste(tx_id, start, end) %in% event_coords)
  return(gr)
}

# for generate_retained_intron
utils::globalVariables(c("estimate", "strand", "coord_sig"))

#' Fuse exons 2 and 3 of a transcript, retaining the intron between them
#'
#' Shared implementation for generate_ri() / generate_si(). One transcript
#' on the edited side is drawn per gene, then up to `n_events` of those
#' genes, and in each one
#' exon rank 2 is extended over rank 3, which is then dropped. The fused
#' exon retains the intron that the transcripts of the opposite sign still
#' splice out, and that is where `sim_event` is set.
#'
#' Without `inverse` an up transcript retains the intron (ri); with
#' `inverse` a down transcript does, so the up transcripts splice it out
#' (si). Either way the event is reported on the fused exon itself.
#'
#' @param gr A GRanges object with an 'estimate' column
#' @param n_events Number of events to generate, at most one per gene
#' @param inverse If TRUE, fuses the exons of a down transcript
#' @return A GRanges object with the events introduced
#' @noRd
generate_retained_intron <- function(gr, n_events = 1, inverse = FALSE) {
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }
  # 1: an up transcript retains the intron; -1: a down one does
  factor <- if (inverse) -1 else 1

  # generate retained introns by creating a new exon that
  # starts with exonrank x and ends at exons rank x+1 for each transcript
  # only in transcripts on the edited side
  # then remove the original exons being retained and
  # add the new retained intron exon
  # and re-rank accordingly
  ri_tx_ids <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate * factor > 0) |>
    dplyr::distinct(gene_id, tx_id) |>
    sample_per_gene(n_events) |>
    dplyr::pull(tx_id)

  exon_idx <- 2
  new_txps_with_ri <- gr |>
    dplyr::filter(tx_id %in% ri_tx_ids) |>
    #get the exons that are being retained
    # get exon 2 and merged with exon 3 to create the retained intron
    # then remove exon 2 and 3 and re rank the exons accordingly
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      # the next exon by rank sits to the genomic right on "+" but to the
      # genomic left on "-", so the intron is swallowed by pushing out the
      # end or the start accordingly. extending the end on "-" would run
      # the range backwards and error out.
      end = ifelse(
        exon_rank == exon_idx & strand != "-",
        end[exon_rank == exon_idx + 1],
        end
      ),
      start = ifelse(
        exon_rank == exon_idx & strand == "-",
        start[exon_rank == exon_idx + 1],
        start
      )
    ) |>
    dplyr::filter(!(exon_rank %in% c(exon_idx + 1))) |>
    dplyr::ungroup() |>
    plyranges::as_granges()

  ri_event_coords <- new_txps_with_ri |>
    as.data.frame() |>
    dplyr::filter(exon_rank == exon_idx) |>
    dplyr::mutate(coord_sig = paste(start, end, gene_id)) |>
    dplyr::pull(coord_sig)

  gr <- gr |>
    dplyr::filter(!tx_id %in% ri_tx_ids) |>
    plyranges::bind_ranges(new_txps_with_ri)
  gr <- rerank_exons(gr)

  # the event is the fused exon itself, on the edited side
  gr <- gr |>
    dplyr::mutate(
      sim_event = estimate * factor > 0 &
        paste(start, end, gene_id) %in% ri_event_coords
    )
  return(gr)
}

#' @rdname generate_events
#' @return `generate_ri()`: A GRanges object with retained intron
#' events introduced: exons 2 and 3 of one up transcript per gene are fused
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_ri(gr, n_events = 1)
#'
#' @export
generate_ri <- function(gr, n_events = 1) {
  generate_retained_intron(gr, n_events)
}

#' @rdname generate_events
#' @return `generate_si()`: A GRanges object with spliced intron events
#' introduced: exons 2 and 3 of one down transcript per gene are fused, so
#' the up transcripts splice out the intron it retains
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_si(gr, n_events = 1)
#'
#' @export
generate_si <- function(gr, n_events = 1) {
  generate_retained_intron(gr, n_events, inverse = TRUE)
}

# for rerank_exons
utils::globalVariables(c("strand"))

#' Re-rank exons in a GRanges object
#' @param gr A GRanges object with metadata columns: 'exon_rank'
#' @return A GRanges object with re-ranked exons
#' @noRd
rerank_exons <- function(gr) {
  gr <- gr |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      exon_rank = seq_len(dplyr::n()),
      exon_rank = dplyr::if_else(
        strand == "-",
        rev(exon_rank),
        exon_rank
      )
    ) |>
    dplyr::ungroup()
  # recalculate internal column, key and  nexons
  gr <- preprocess(gr, coef_col = "estimate")
  return(gr)
}

# for generate_a5ss
utils::globalVariables(c("estimate", "internal", "key"))

#' @rdname generate_events
#' @return `generate_a5ss()`: A GRanges object with alternative 5' splice site
#' events introduced: the donor of an internal exon is moved by 2 bp in one
#' up transcript per gene
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_a5ss(gr, n_events = 1)
#'
generate_a5ss <- function(gr, n_events = 1) {
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }

  # generate a5ss by modifying end() of random internal
  # exons in transcripts with estimate > 0, one exon per gene
  a5ss_exon_key <- gr |>
    as.data.frame() |>
    # TO DO: include first and last exons
    dplyr::filter(estimate > 0 & internal == TRUE) |>
    dplyr::distinct(gene_id, key) |>
    sample_per_gene(n_events) |>
    dplyr::pull(key)

  # the donor (5' splice site) is the exon end on "+" but the exon start
  # on "-", so move whichever boundary the strand puts downstream.
  # as.character() because an ungrouped mutate() hands back strand as an
  # Rle, which if_else() rejects as a condition. the offset is named
  # bp_shift, not shift: plyranges' data mask resolves a bare `shift` to
  # IRanges::shift(), so `end + shift` adds a function and errors.
  bp_shift <- sample(c(-2L, 2L), 1)
  on_minus <- as.character(GenomicRanges::strand(gr)) == "-"
  gr_with_a5ss <- gr |>
    dplyr::mutate(
      end = dplyr::if_else(
        key %in% a5ss_exon_key & !on_minus,
        end + bp_shift,
        end
      ),
      start = dplyr::if_else(
        key %in% a5ss_exon_key & on_minus,
        start + bp_shift,
        start
      )
    )
  gr_with_a5ss <- preprocess(gr_with_a5ss, coef_col = "estimate")
  gr_with_a5ss <- gr_with_a5ss |>
    dplyr::mutate(sim_event = key %in% a5ss_exon_key)
  return(gr_with_a5ss)
}
# for generate_a3ss
utils::globalVariables(c("estimate", "internal", "key"))

#' @rdname generate_events
#' @return `generate_a3ss()`: A GRanges object with alternative 3' splice site
#' events introduced: the acceptor of an internal exon is moved by 2 bp in
#' one up transcript per gene
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_a3ss(gr, n_events = 1)
#'
generate_a3ss <- function(gr, n_events = 1) {
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }

  # generate a3ss by modifying start() of random internal
  # exons in transcripts with estimate > 0, one exon per gene
  a3ss_exon_key <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate > 0 & internal == TRUE) |>
    dplyr::distinct(gene_id, key) |>
    sample_per_gene(n_events) |>
    dplyr::pull(key)

  # the acceptor (3' splice site) is the exon start on "+" but the exon
  # end on "-", so move whichever boundary the strand puts upstream.
  # as.character() because an ungrouped mutate() hands back strand as an
  # Rle, which if_else() rejects as a condition. the offset is named
  # bp_shift, not shift: plyranges' data mask resolves a bare `shift` to
  # IRanges::shift(), so `end + shift` adds a function and errors.
  bp_shift <- sample(c(-2L, 2L), 1)
  on_minus <- as.character(GenomicRanges::strand(gr)) == "-"
  gr_with_a3ss <- gr |>
    dplyr::mutate(
      start = dplyr::if_else(
        key %in% a3ss_exon_key & !on_minus,
        start + bp_shift,
        start
      ),
      end = dplyr::if_else(
        key %in% a3ss_exon_key & on_minus,
        end + bp_shift,
        end
      )
    )
  gr_with_a3ss <- preprocess(gr_with_a3ss, coef_col = "estimate")
  gr_with_a3ss <- gr_with_a3ss |>
    dplyr::mutate(sim_event = key %in% a3ss_exon_key)
  return(gr_with_a3ss)
}

# for generate_alt_ts
utils::globalVariables(c(
  "estimate", "key", "terminal", "exon_rank", "tx_id", "sim_event"
))

#' Introduce an alternative transcription start or end site
#'
#' Shared implementation for generate_atss() / generate_ates(). Two ways
#' of moving where a transcript begins or ends:
#'
#' - `"shift"` moves the *outer* boundary of the terminal exon (the one
#'   facing away from the rest of the transcript), leaving the inner
#'   splice site untouched. The altered exon still overlaps its partner.
#' - `"drop"` removes the terminal exon outright and re-ranks, so the
#'   transcript now begins (or ends) at what was the next exon along.
#'   That exon does not overlap the partner's terminal exon at all, which
#'   is the case an overlap-based finder could never reach.
#'
#' @param gr A GRanges object with an 'estimate' column
#' @param n_events Number of events to generate, at most one per gene
#' @param at_start If TRUE, alters the first exon; if FALSE, the last
#' @param mode Either `"shift"` or `"drop"`, as above
#' @return A GRanges object with the events introduced
#' @noRd
generate_alt_ts <- function(gr, n_events = 1, at_start = TRUE,
                            mode = c("shift", "drop")) {
  mode <- match.arg(mode)
  # if preprocessing didn't happen
  if (
    !all(
      c("key", "nexons", "internal") %in%
        names(GenomicRanges::mcols(gr))
    )
  ) {
    gr <- preprocess(gr, coef_col = "estimate")
  }

  # the first / last exon of each positive-coefficient transcript, then one
  # of those per gene. the 'internal' mcol is no help: it marks both
  # terminal exons at once. min/max rather than 1/nexons so CDS-style
  # offset ranks still work.
  alt_exon_key <- tibble::as_tibble(gr) |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      terminal = if (at_start) {
        exon_rank == min(exon_rank)
      } else {
        exon_rank == max(exon_rank)
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(estimate > 0 & terminal) |>
    dplyr::distinct(gene_id, key) |>
    sample_per_gene(n_events) |>
    dplyr::pull(key)

  if (mode == "drop") {
    # keys are rebuilt by the re-rank below, so remember the transcripts
    # rather than the keys
    alt_tx <- tibble::as_tibble(gr) |>
      dplyr::filter(key %in% alt_exon_key) |>
      dplyr::pull(tx_id)

    gr_with_alt <- gr |>
      dplyr::filter(!key %in% alt_exon_key)
    gr_with_alt <- rerank_exons(gr_with_alt)

    # the event is on whichever exon became terminal
    gr_with_alt <- gr_with_alt |>
      dplyr::group_by(tx_id) |>
      dplyr::mutate(
        sim_event = tx_id %in% alt_tx &
          if (at_start) {
            exon_rank == min(exon_rank)
          } else {
            exon_rank == max(exon_rank)
          }
      ) |>
      dplyr::ungroup()
    return(gr_with_alt)
  }

  # mode == "shift": which boundary is the outer one depends on both the
  # strand and which end of the transcript this is. a first exon begins
  # the transcript at its genomic start on "+" but at its genomic end on
  # "-", and a last exon is the reverse.
  on_minus <- as.character(GenomicRanges::strand(gr)) == "-"
  outer_is_start <- xor(at_start, on_minus)

  # always moved *inward*, shrinking the exon: outward could push the
  # first exon of the first gene below coordinate 1, and inward can never
  # collide with the neighbouring exon.
  gr_with_alt <- gr |>
    dplyr::mutate(
      start = dplyr::if_else(
        key %in% alt_exon_key & outer_is_start,
        start + 2L,
        start
      ),
      end = dplyr::if_else(
        key %in% alt_exon_key & !outer_is_start,
        end - 2L,
        end
      )
    )
  gr_with_alt <- preprocess(gr_with_alt, coef_col = "estimate")
  gr_with_alt |>
    dplyr::mutate(sim_event = key %in% alt_exon_key)
}

#' @rdname generate_events
#' @param mode How to move the site: `"shift"` moves the terminal exon's
#'   outer boundary (so it still overlaps its partner), `"drop"` removes
#'   the terminal exon and re-ranks (so the new terminal exon does not
#'   overlap its partner at all).
#' @return `generate_atss()`: A GRanges object with alternative
#' transcription start site events introduced: the first exon of one up
#' transcript per gene is shortened at its outer end (`"shift"`) or removed
#' (`"drop"`)
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_atss(gr, n_events = 1)
#'
#' # the non-overlapping variant: drop the first exon altogether
#' generate_atss(gr, n_events = 1, mode = "drop")
#'
generate_atss <- function(gr, n_events = 1, mode = c("shift", "drop")) {
  generate_alt_ts(gr, n_events, at_start = TRUE, mode = mode)
}

#' @rdname generate_events
#' @return `generate_ates()`: A GRanges object with alternative
#' transcription end site events introduced: the last exon of one up
#' transcript per gene is shortened at its outer end (`"shift"`) or removed
#' (`"drop"`)
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_ates(gr, n_events = 1)
#'
generate_ates <- function(gr, n_events = 1, mode = c("shift", "drop")) {
  generate_alt_ts(gr, n_events, at_start = FALSE, mode = mode)
}
