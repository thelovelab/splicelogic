
#' Create a sample GRanges with one negative coef transcript and two positive coef transcripts
#' This dataset is designed to include a skipped exon event at exon_rank 3 and 5 of tx_id 1 
#' happening between both tx_id 2 and tx_id 3
#' 
#' @return A GRanges object with two transcripts per gene and candidate logic
#' @import GenomicRanges
#' @importFrom plyranges as_granges bind_ranges
#' @importFrom dplyr mutate case_when
#' @return A GRanges object with two transcripts per gene
#' @export
se_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41, 51, 61),
    width = 5,
    strand = "+",
    exon_rank = 1:7,
    gene_id = rep(1, 7),
    coefs = rep(runif(1, min = -1, max = 0), 7)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = 1:6,
    gene_id = rep(1, 6),
    coefs = rep(runif(1, min = 0, max = 1), 6)
  )
    df3 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  31, 41, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = 1:6,
    gene_id = rep(1, 6),
    coefs = rep(runif(1, min = 0, max = 1), 6)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr3 <- plyranges::as_granges(df3)
  gr <- plyranges::bind_ranges(gr1, gr2, gr3) |>
      plyranges::mutate(
      tx_id = c(rep(1, 7), rep(2, 6), rep(3, 6)))

return(gr)
}

#' Create a sample GRanges with one negative coef transcript and one positive coef transcripts
#' This dataset is designed to include a mutually exclusive exons between exon_rank 3 of tx_id 1
#' and exon_rank 3 of tx_id 2. There is also a skipped exon event at exon_rank 5 of tx_id 1.
#' 
#' @return A GRanges object with two transcripts per gene and candidate logic
#' @import GenomicRanges
#' @importFrom plyranges as_granges bind_ranges
#' @importFrom dplyr mutate case_when
#' @return A GRanges object with two transcripts per gene
#' @export
mx_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 41, 51, 61, 71, 91, 101),
    width = 5,
    strand = "+",
    exon_rank = 1:9,
    gene_id = rep(1, 9),
    coefs = rep(runif(1, min = -1, max = 0), 9)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 41, 61, 71, 81, 101),
    width = 5,
    strand = "+",
    exon_rank = 1:8,
    gene_id = rep(1, 8),
    coefs = rep(runif(1, min = 0, max = 1), 8)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
      plyranges::mutate(
      tx_id = c(rep(1, 9), rep(2, 8)))

return(gr)
}
#' Create a sample GRanges with one negative coef transcript and one positive coef transcripts
#' This dataset is designed to include no splicing events
#' @return A GRanges object with two transcripts per gene and no candidate logic
#' @import GenomicRanges
#' @importFrom plyranges as_granges bind_ranges
#' @importFrom dplyr mutate case_when
#' @return A GRanges object with two transcripts per gene
#' @export  
no_event_mock_data <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21),
    width = 5,
    strand = "+",
    exon_rank = 1:3,
    gene_id = rep(1, 3),
    coefs = rep(runif(1, min = -1, max = 0), 3)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11,  21),
    width = 5,
    strand = "+",
    exon_rank = 1:3,
    gene_id = rep(1, 3),
    coefs = rep(runif(1, min = 0, max = 1), 3)
  )
  gr1 <- plyranges::as_granges(df1)
  gr2 <- plyranges::as_granges(df2)
  gr <- plyranges::bind_ranges(gr1, gr2) |>
      plyranges::mutate(
      tx_id = c(rep(1, 3), rep(2, 3)))

return(gr)
}


create_mock_data <- function(n_genes = 1,
                              n_tx_per_gene = 2, 
                              n_exons_per_tx = 5, 
                              coef_range = c(-1, 1)) {
  # Generate all combinations of genes, transcripts, and exons
  data <- expand.grid(
    gene_id = seq_len(n_genes),
    tx_id = seq_len(n_tx_per_gene),
    exon_rank = seq_len(n_exons_per_tx)
  )
  # browser()

data <- data %>%
  dplyr::arrange(gene_id, tx_id, exon_rank)

  # Calculate transcript IDs globally
  data <- data |>
    dplyr::mutate(
      tx_id = tx_id + (gene_id - 1) * n_tx_per_gene,
      seqnames = paste0("chr", sample(1:22, 1)), # Random chromosome
      start = (exon_rank - 1) * 10 + 1,         # Start positions
      width = 5,                                # Fixed width
      #strand = rep(sample(c("+", "-"), n_genes * n_tx_per_gene, replace = TRUE), each = n_exons_per_tx)
      strand = "+" # Fixed strand for simplicity TO DO 
    )

  # Reverse exon_rank for transcripts with strand == "-"
  data <- data |>
    dplyr::group_by(tx_id) |>
    dplyr::mutate(
      exon_rank = ifelse(strand == "-", rev(exon_rank), exon_rank)
    ) |>
    dplyr::ungroup()
  
  # Assign random coefs per tx_id and gene_id
  data <- data |>
    dplyr::group_by(gene_id, tx_id) |>
    dplyr::mutate(
      coefs = runif(1, min = coef_range[1], max = coef_range[2]) # Same value for the group
    ) |>
    dplyr::ungroup()
  
  # Convert to GRanges
  gr <- plyranges::as_granges(data)
  return(gr)
}

generate_retained_introns <- function(gr, n_ri = 1) {
  set.seed(123) # for reproducibility
  # generate retained introns by creating a new exon that
  # starts with exonrank x and ends at exons rank x+1 for each transcript
  # only in transcripts with coefs > 0
  # then remove the original exons that are being retained and add the new retained intron exon with coefs = 0
  # and re-rank accordingly
  ri_tx_ids <- gr |>
    as.data.frame() |>
    dplyr::filter(coefs > 0) |>
    dplyr::distinct(tx_id) |>
    dplyr::slice_sample(n = n_ri) |>
    dplyr::pull(tx_id)

  exon_idx <- 2
  new_txps_with_ri <- gr |>
    plyranges::filter(tx_id %in% ri_tx_ids) |>
    #get the exons that are being retained
    # get exon 2 and merged with exon 3 to create the retained intron
    # then remove exon 2 and 3 and re rank the exons accordingly
    plyranges::group_by(tx_id) |>
    plyranges::mutate(
      end = ifelse(exon_rank == exon_idx, end[exon_rank == exon_idx + 1], end)
    ) |>
    dplyr::filter(!(exon_rank %in% c(exon_idx + 1))) |>
    # re rank the exons accordingly
    dplyr::mutate(
      exon_rank = seq_len(dplyr::n())
    ) |>
    dplyr::ungroup() |>
    plyranges::as_granges() 

  #remove old transcripts with retained introns and add the new ones
  gr <- gr |>
    plyranges::filter(!tx_id %in% ri_tx_ids) |>
    plyranges::bind_ranges(new_txps_with_ri) 
  return(gr)
}

generate_a5ss <- function(gr, n_a5ss = 1) {
    # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col ="coefs")
  }
  set.seed(123) # for reproducibility
  # generate alternative 5' splice sites by modifying the end() of random internal = TRUE exons in transcripts with coefs > 0
  a5ss_exon_key <- gr |> 
    as.data.frame() |>
    dplyr::filter(coefs > 0 & internal == TRUE) |> #TO DO : include that it can be first and last exons
    dplyr::distinct(key) |>
    dplyr::slice_sample(n = n_a5ss) |>
    dplyr::pull(key)
  
  gr_with_a5ss <- gr  |>
    plyranges::mutate(
      end = dplyr::if_else(key %in% a5ss_exon_key, end + sample(c(-2, 2), 1), end) # TO DO : include - strand case (start instead of end)
    )
    return(gr_with_a5ss)
}

generate_a3ss <- function(gr, n_a3ss = 1) {
    # if preprocessing didn't happen
  if (!all(c("key", "nexons", "internal", "event") %in% names(mcols(gr)))) {
    gr <- preprocess_input(gr, coef_col)
  }
  set.seed(123) # for reproducibility
  # generate alternative 3' splice sites by modifying the end() of random internal = TRUE exons in transcripts with coefs > 0
  a3ss_exon_key <- gr |> 
    as.data.frame() |>
    dplyr::filter(coefs > 0 & internal == TRUE) |>
    dplyr::distinct(key) |>
    dplyr::slice_sample(n = n_a3ss) |>
    dplyr::pull(key)
  
  gr_with_a3ss <- gr  |>
    plyranges::mutate(
      start = dplyr::if_else(key %in% a3ss_exon_key, start + sample(c(-2, 2), 1), start) ) # TO DO : include - strand case (end instead of start)

    return(gr_with_a3ss)
}

