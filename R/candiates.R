intersect_bp <- function(start, stop, ref_start, ref_stop) {
  # browser()
  # ref:              [----]
  # others:   [----]              # not compared
  #               [-----]
  #                   [----]
  #                      [-----]
  to_compare <- (ref_start <= stop & ref_stop >= start)
  #| (ref_start > start & ref_stop > stop)
  #  only look at positions that will overlap
  # res <- numeric(length(start))
  # uni -> union
  # int -> intersection
  to_compare <- which(to_compare)
  uni_start <- int_start <- start[to_compare]
  # ref:              [----]
  # int_start:
  #               [-->[
  #                   [
  #                      [
  int_start[int_start < ref_start] <- ref_start
  uni_stop <- int_stop <- stop[to_compare]
  # ref:              [----]
  # int_stop
  #                     ]
  #                        ]
  #                        ]<--]
  int_stop[int_stop > ref_stop] <- ref_stop
  # ref:              [----]
  # int_ranges:
  #                   [-]
  #                   [----]
  #                      [-]

  # ref:              [----]
  # uni_start:
  #               [
  #                   [
  #                   [<-[
  uni_start[uni_start > ref_start] <- ref_start
  # ref:              [----]
  # uni_stop:
  #                     ]->]
  #                        ]
  #                            ]
  uni_stop[uni_stop < ref_stop] <- ref_stop

  # ref:              [----]
  # uni_ranges:
  #               [--------]
  #                   [----]
  #                   [-------]

  # ref:              [----]
  # res:
  #   comp1:
  #                   [-]
  #               [--------]
  #   comp2:
  #                   [----]
  #                   [----]
  #   comp3:
  #                      [-]
  #                   [-------]

  # res[to_compare] <- (int_stop - int_start + 1L) / (uni_stop - uni_start + 1L)
  # # if(any(res<0))
  # #   browser()
  # # res
  # res
  list(
    index = to_compare,
    intersection = int_stop - int_start + 1L,
    union = uni_stop - uni_start + 1L
  )
}

ind_split_by_match <- function(x, ux = unique(x)) {
  matches <- match(x, ux)
  o <- order(matches)
  split(o, matches[o])
}

gr_percent_similarity <- function(gr_row, gr_col) {
  n_row <- length(gr_row)
  n_col <- length(gr_col)
  mat <- matrix(0, n_row, n_col)
  # some exons are identical,
  # dont do more work than we need
  unique_row <- unique(gr_row)
  unique_col <- unique(gr_col)
  # cache how each unique index maps back to the
  # original vector index
  uind_row <- ind_split_by_match(gr_row, unique_row)
  uind_col <- ind_split_by_match(gr_col, unique_col)
  start_row <- start(unique_row)
  stop_row <- end(unique_row)
  start_col <- start(unique_col)
  stop_col <- end(unique_col)

  row_map_lens <- vapply(uind_row, length, 1L)
  for (j in seq_len(length(unique_col))) {
    ref_start <- start_col[j]
    ref_stop <- stop_col[j]
    data <- intersect_bp(
      start_row, stop_row,
      ref_start = ref_start, ref_stop = ref_stop
    )
    # data$index maps to the indices of uind_row
    row_vec <- data$intersection / data$union
    # matrices recycle by row
    row_vec <- rep(row_vec, row_map_lens[data$index])
    row_map <- unlist(uind_row[data$index])
    mat[row_map, uind_col[[j]]] <- row_vec
  }
  mat
}

plot_gr <- function(gr) {
  data <- data.frame(
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    center = (GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2,
    width = GenomicRanges::width(gr),
    id = as.factor(paste(gr$gene_id, gr$tx_id)) |> as.integer()
  )
  ggplot(data, aes(x = start, y = id)) +
    geom_segment(
      data = \(x) {
        group_by(x, id) |> summarise(start = min(start), end = max(end))
      },
      aes(xend = end)
    ) +
    geom_rect(
      aes(
        xmin = start, xmax = end,
        ymin = id - 0.25, ymax = id + 0.25
      ),
      fill = "blue"
    )
}

se_mock_data2 <- function() {
  df1 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 21, 31, 41, 51, 61, 81),
    width = 5,
    strand = "+",
    exon_rank = 1:8,
    gene_id = rep(1, 8),
    coefs = rep(runif(1, min = -1, max = 0), 8)
  )
  df2 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 61, 71),
    width = 5,
    strand = "+",
    exon_rank = 1:5,
    gene_id = rep(1, 5),
    coefs = rep(runif(1, min = 0, max = 1), 5)
  )
  df3 <- data.frame(
    seqnames = "chr1",
    start = c(1, 11, 31, 41, 61, 71),
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
      tx_id = c(rep(1, 8), rep(2, 5), rep(3, 6))
    )

  return(gr)
}

# initialzie a matrix with values set to the intersections
# of i and j.
mat_ij <- function(i, j, nrow, ncol) {
  shift <- (j - 1L) * nrow
  out <- matrix(0, nrow = nrow, ncol = ncol)
  out[i + shift] <- 1L
  out
}

vec_shift <- function(x, n) {
  data.frame(
    head = head(x, n = -n),
    tail = tail(x, n = -n)
  )
}

left_shift_by_id <- function(id) {
  # browser()
  id <- match(id, unique(id))
  seq <- split(seq_along(id), id) |>
    unname() |>
    lapply(vec_shift, n = 1L) |>
    vctrs::list_unchop()
  n <- length(id)
  mat_ij(
    seq$tail,
    seq$head,
    nrow = n,
    ncol = n
  )
}

right_shift_by_id <- function(id) {
  # browser()
  id <- match(id, unique(id))
  seq <- split(seq_along(id), id) |>
    unname() |>
    lapply(vec_shift, n = 1L) |>
    vctrs::list_unchop()
  n <- length(id)
  mat_ij(
    seq$head,
    seq$tail,
    nrow = n,
    ncol = n
  )
}



mat_col_groups <- function(id) {
  id <- match(id, unique(id))
  mat_ij(i = seq_along(id), j = id, nrow = length(id), ncol = max(id))
}

#' Find candidates skipped exons
#' @description
#' computes similarity between positive and negative GenomicRanges
#' and then finds candidate exons in neg whose adjacent exons
#' exist in the positive set and are also adjacent
#' @examples
#' gr <- se_mock_data()
#' library(ggplot2)
#' library(GenomicRanges)
#' library(plyranges)
#' plot_gr(gr)
#' gr_pos <- filter(gr, coefs > 0)
#' gr_neg <- filter(gr, coefs < 0)
#' find_candidates(gr_pos, gr_neg, gr_pos$tx_id, gr_neg$tx_id)
#' @export
find_candidates <- function(gr_pos, gr_neg, id_pos, id_neg) {
  # returns values between [0--1]
  # candidates encodes some level of percent similarity between
  # the ranges in the positive GR (rows) and the negative GR (cols)
  candidates <- gr_percent_similarity(gr_pos, gr_neg)
  # converts this to a matrix indexed by group of gr_pos
  candidates <- t(mat_col_groups(id_pos)) %*% candidates
  # shifts columns to the left within each group.
  # i.e. the right exon (if it exists) is not in the posistion
  # of the candidate exon
  cand_right <- candidates %*% left_shift_by_id(id_neg)
  # shifts columns to the right within each group
  # # i.e. the left exon (if it exists) is not in the posistion
  # of the candidate exon
  cand_left <- candidates %*% right_shift_by_id(id_neg)
  n_pos <- nrow(candidates)
  uni_pos <- unique(id_pos)
  is_cand <- colSums(candidates > 0) < n_pos
  cand_mat <- cand_right * cand_left

  w <- which(cand_mat > 0)
  which_pos <- ((w + n_pos - 1L) %% n_pos) + 1L
  out <- rep(gr_neg[is_cand], colSums(cand_mat)[is_cand])
  out$tx_event <- uni_pos[which_pos]
  out$event <- "skipped_exon"
  out
}
