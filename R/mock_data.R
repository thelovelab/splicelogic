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

# for create_mock_data
utils::globalVariables(c("gene_id", "strand", "tx_order", "estimate"))

#' Create mock GRanges data for splicing event testing
#' 
#' @param n_genes Number of genes to simulate
#' @param n_tx_per_gene Number of transcripts per gene
#' @param n_exons_per_tx Number of exons per transcript
#' @param coef_range Range of coefficient values to sample from
#' 
#' @return A GRanges object with simulated transcripts and exons
#' 
#' @examples
#' 
#' # create mock data with 2 genes, 4 transcripts
#' # per gene, and 4 exons per transcript
#' gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4)
#' 
#' @export
create_mock_data <- function(
  n_genes = 1,
  n_tx_per_gene = 2,
  n_exons_per_tx = 5,
  coef_range = c(-1, 1)
) {
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
      strand = "+" # Fixed strand for simplicity TO DO
    )

  # Reverse exon_rank for transcripts with strand == "-"
  data <- data |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      exon_rank = ifelse(strand == "-", rev(exon_rank), exon_rank)
    ) |>
    dplyr::ungroup()

  # Assign estimate per tx_id ensuring at least 1 neg and 1 pos per gene
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
utils::globalVariables(c("estimate"))

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
      end = ifelse(exon_rank == exon_idx, end[exon_rank == exon_idx + 1], end)
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

  gr_with_a5ss <- gr |>
    dplyr::mutate(
      end = dplyr::if_else(
        key %in% a5ss_exon_key,
        end + sample(c(-2, 2), 1),
        end
      ) # TO DO : include - strand case (start instead of end)
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

  gr_with_a3ss <- gr |>
    dplyr::mutate(
      start = dplyr::if_else(
        key %in% a3ss_exon_key,
        start + sample(c(-2, 2), 1),
        start
      )
    ) # TO DO : include - strand case (end instead of start)
  gr_with_a3ss <- preprocess(gr_with_a3ss, coef_col = "estimate")
  gr_with_a3ss <- gr_with_a3ss |>
    dplyr::mutate(sim_event = key %in% a3ss_exon_key)
  return(gr_with_a3ss)
}
