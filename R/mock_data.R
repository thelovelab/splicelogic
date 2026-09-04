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
#' Six one-gene loci, all on `-`, each engineered for a single event. What
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
#' ```
#'
#' In g4 the moved boundary is the exon start, which on `-` is the donor,
#' so it is an a5ss; in g5 it is the exon end, the acceptor, so it is an
#' a3ss. Both are the mirror image of the `+` strand case.
#' @return A GRanges object with twelve transcripts in six genes
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
    )
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
#' @param gr A GRanges object with metadata columns: 'exon_rank', 'gene_id',
#' 'tx_id', and 'estimate'.
NULL

# for generate_se
utils::globalVariables(c("estimate", "internal", "key", "gene_id", "sim_event"))

#' @rdname generate_events
#' @param n_events Number of events to generate
#' @return `generate_se()`: A GRanges object with skipped exon
#' events introduced
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_se(gr, n_events = 1)
#'
generate_se <- function(gr, n_events = 1) {
  # one candidate key per gene (sampled from internal positive-estimate exons),
  # then draw up to n_events genes from that pool
  candidates <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate > 0 & internal == TRUE) |>
    dplyr::distinct(gene_id, key) |>
    dplyr::group_by(gene_id) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup()

  n_available <- nrow(candidates)
  if (n_events > n_available) {
    warning(sprintf(
      "n_events (%d) exceeds available genes with candidates (%d); capping at %d.",
      n_events, n_available, n_available
    ))
    n_events <- n_available
  }

  se_exons_key <- candidates |>
    dplyr::slice_sample(n = n_events) |>
    dplyr::pull(key)

  event_coords <- gr |>
    as.data.frame() |>
    dplyr::filter(key %in% se_exons_key) |>
    dplyr::mutate(coord_sig = paste(start, end, gene_id)) |>
    dplyr::pull(coord_sig)

  gr <- gr |>
    dplyr::filter(!key %in% se_exons_key)
  gr <- rerank_exons(gr)

  gr <- gr |>
    dplyr::mutate(
      sim_event = estimate < 0 & paste(start, end, gene_id) %in% event_coords
    )
  return(gr)
}

# for generate_mxe
utils::globalVariables(c(
  "estimate", "internal", "next_key", "gene_id", "neg_tx_id", "key"
))

#' @rdname generate_events
#' @param n_events Number of events to generate
#' @return `generate_mxe()`: A GRanges object with mutually exclusive exon
#' events introduced
#' @export
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_mxe(gr, n_events = 1)
#'
generate_mxe <- function(gr, n_events = 1) {
  # Find consecutive pairs of internal exons from neg transcripts
  # Both rank k and rank k+1 must be internal
  neg_internal <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate < 0 & internal == TRUE)

  # Find rank k where rank k+1 is also internal in same transcript
  mx_candidates <- neg_internal |>
    dplyr::mutate(
      next_key = paste0(tx_id, "-", exon_rank + 1L)
    ) |>
    dplyr::filter(next_key %in% neg_internal$key) |>
    dplyr::slice_sample(n = n_events)

  # Pair each selected neg transcript with one pos transcript from same gene
  gr_df <- as.data.frame(gr)

  pos_txs <- gr_df |>
    dplyr::filter(estimate > 0) |>
    dplyr::distinct(gene_id, tx_id)

  mx_pairs <- mx_candidates |>
    dplyr::select(gene_id, neg_tx_id = tx_id, exon_rank) |>
    dplyr::inner_join(pos_txs, by = "gene_id") |>
    dplyr::rename(pos_tx_id = tx_id) |>
    dplyr::group_by(neg_tx_id, exon_rank) |>
    dplyr::slice_sample(n = 1) |>
    dplyr::ungroup()

  # Build keys to remove (pairwise):
  # rank k from the selected pos transcript (pos loses k, keeps k+1)
  # rank k+1 from the selected neg transcript (neg keeps k, loses k+1)
  pos_remove <- paste0(mx_pairs$pos_tx_id, "-", mx_pairs$exon_rank)
  neg_remove <- paste0(mx_pairs$neg_tx_id, "-", mx_pairs$exon_rank + 1L)

  # pos-removed exon coords survive in neg; neg-removed exon coords survive in pos
  pos_removed_coords <- gr |>
    as.data.frame() |>
    dplyr::filter(key %in% pos_remove) |>
    dplyr::mutate(coord_sig = paste(start, end, gene_id)) |>
    dplyr::pull(coord_sig)
  neg_removed_coords <- gr |>
    as.data.frame() |>
    dplyr::filter(key %in% neg_remove) |>
    dplyr::mutate(coord_sig = paste(start, end, gene_id)) |>
    dplyr::pull(coord_sig)

  gr <- gr |>
    dplyr::filter(!key %in% c(pos_remove, neg_remove))
  gr <- rerank_exons(gr)

  gr <- gr |>
    dplyr::mutate(
      sim_event = (estimate < 0 & paste(start, end, gene_id) %in% pos_removed_coords) |
                  (estimate > 0 & paste(start, end, gene_id) %in% neg_removed_coords)
    )
  return(gr)
}

# for generate_ri
utils::globalVariables(c("estimate", "strand"))

#' @rdname generate_events
#' @param n_events Number of events to generate
#' @return `generate_ri()`: A GRanges object with retained intron
#' events introduced
#' @examples
#'
#' gr <- create_mock_data(
#'   n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
#' )
#' generate_ri(gr, n_events = 1)
#'
#' @export
generate_ri <- function(gr, n_events = 1) {
  # generate retained introns by creating a new exon that
  # starts with exonrank x and ends at exons rank x+1 for each transcript
  # only in transcripts with estimate > 0
  # then remove the original exons being retained and
  # add the new retained intron exon with estimate = 0
  # and re-rank accordingly
  ri_tx_ids <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate > 0) |>
    dplyr::distinct(tx_id) |>
    dplyr::slice_sample(n = n_events) |>
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

  gr <- gr |>
    dplyr::mutate(
      sim_event = estimate > 0 & paste(start, end, gene_id) %in% ri_event_coords
    )
  return(gr)
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
#' @param n_events Number of events to generate
#' @return `generate_a5ss()`: A GRanges object with alternative 5' splice site
#' events introduced
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
  # exons in transcripts with estimate > 0
  a5ss_exon_key <- gr |>
    as.data.frame() |>
    # TO DO: include first and last exons
    dplyr::filter(estimate > 0 & internal == TRUE) |>
    dplyr::distinct(key) |>
    dplyr::slice_sample(n = n_events) |>
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
#' @param n_events Number of events to generate
#' @return `generate_a3ss()`: A GRanges object with alternative 3' splice site
#' events introduced
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
  # exons in transcripts with estimate > 0
  a3ss_exon_key <- gr |>
    as.data.frame() |>
    dplyr::filter(estimate > 0 & internal == TRUE) |>
    dplyr::distinct(key) |>
    dplyr::slice_sample(n = n_events) |>
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
#' @param n_events Number of events to generate
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

  # the first / last exon of each positive-coefficient transcript. the
  # 'internal' mcol is no help: it marks both terminal exons at once.
  # min/max rather than 1/nexons so CDS-style offset ranks still work.
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
    dplyr::distinct(key) |>
    dplyr::slice_sample(n = n_events) |>
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
#' @param n_events Number of events to generate
#' @param mode How to move the site: `"shift"` moves the terminal exon's
#'   outer boundary (so it still overlaps its partner), `"drop"` removes
#'   the terminal exon and re-ranks (so the new terminal exon does not
#'   overlap its partner at all).
#' @return `generate_atss()`: A GRanges object with alternative
#' transcription start site events introduced
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
#' transcription end site events introduced
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
