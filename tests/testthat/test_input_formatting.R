# Tests for preprocess function
test_that("preprocess adds key, nexons, internal columns", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 4)
  result <- preprocess(gr, coef_col = "estimate")
  expect_true(all(
    c("key", "nexons", "internal") %in%
      names(GenomicRanges::mcols(result))
  ))
})

test_that("preprocess doesn't wipe out original metadata", {
  gr <- se_mock_data()
  S4Vectors::metadata(gr) <- list(test=123) # add some starter metadata
  gr_pp <- preprocess(
    gr, coef_col = "estimate"
  )
  expect_equal(
    S4Vectors::metadata(gr_pp)$test, 123
  )
})

test_that("preprocess sets splicelogic_preprocessed flag", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 3)
  result <- preprocess(gr, coef_col = "estimate")
  expect_true(isTRUE(S4Vectors::metadata(result)$splicelogic_preprocessed))
})

test_that("preprocess renames coef_col to estimate", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 3)
  gr$my_coef <- gr$estimate
  result <- preprocess(gr, coef_col = "my_coef")
  expect_true("estimate" %in% names(GenomicRanges::mcols(result)))
})

test_that("preprocess stores method_string in metadata", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 3)
  result <- preprocess(gr, coef_col = "estimate", method_string = "satuRn")
  expect_equal(S4Vectors::metadata(result)$method_string, "satuRn")
})

test_that("preprocess additional_columns -> event_<col> in output", {
  gr <- se_mock_data()

  # add a per-tx, transcript-level label column
  gr$tx_label <- paste0("tx_", gr$tx_id)
  # add another one to test multiple additional columns
  gr$tx_label2 <- paste0("label2_", gr$tx_id)

  gr_pp <- preprocess(
    gr, coef_col = "estimate", additional_columns = c("tx_label", "tx_label2")
  )

  # stored under metadata for downstream consumption
  expect_equal(
    S4Vectors::metadata(gr_pp)$additional_columns, c("tx_label", "tx_label2")
  )

  result <- find_se(gr_pp)
  expect_gt(length(result), 0L)
  expect_true("event_tx_label" %in% names(GenomicRanges::mcols(result)))
  expect_true("event_tx_label2" %in% names(GenomicRanges::mcols(result)))

  # event_tx_label must come from the partner transcript (event_tx_id),
  # not the reference one
  expect_equal(
    as.character(result$event_tx_label),
    paste0("tx_", as.character(result$event_tx_id))
  )
})

test_that("preprocess missing additional_columns throws an error", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 2, n_exons_per_tx = 3)
  expect_error(
    preprocess(gr, coef_col = "estimate", additional_columns = "nonexistent_col"),
    "columns are missing in GRanges metadata: nonexistent_col"
  )
})

# ---------------------------------------------------------------
# Tests for prepare_exons
test_that("prepare_exons errors on missing dtu_table columns", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("AnnotationDbi")
  bad_table <- tibble::tibble(tx_id = "FBtr0001", effect_est = 1.0)
  expect_error(
    prepare_exons(NULL, bad_table, coef_col = "effect_est", verbose = FALSE),
    "Missing columns"
  )
})

test_that("prepare_exons returns GRanges with required columns", {
  skip_if_not_installed("AnnotationHub")
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("AnnotationDbi")
  skip_if_offline()

  ah <- suppressMessages(AnnotationHub::AnnotationHub())
  txdb <- suppressMessages(ah[["AH84134"]])

  txps <- suppressMessages(AnnotationDbi::select(
    txdb,
    keys    = AnnotationDbi::keys(txdb, "TXID"),
    columns = c("TXNAME", "GENEID"),
    keytype = "TXID"
  )) |>
    tibble::as_tibble() |>
    dplyr::select(tx_id = TXNAME, gene_id = GENEID) |>
    dplyr::filter(!is.na(gene_id)) |>
    dplyr::slice_head(n = 20) |>
    dplyr::mutate(effect_est = rnorm(dplyr::n()))

  result <- prepare_exons(txdb, txps, coef_col = "effect_est", verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_true(all(
    c("tx_id", "gene_id", "exon_rank", "effect_est") %in%
      names(GenomicRanges::mcols(result))
  ))
})

# ---------------------------------------------------------------
# Tests for prepare_exons_by_partition
test_that("prepare_exons_by_partition assigns +1 to up and -1 to down", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up   <- gr[gr$estimate > 0]
  gr_down <- gr[gr$estimate < 0]

  result <- prepare_exons_by_partition(gr_up, gr_down, verbose = FALSE)

  expect_s4_class(result, "GRanges")
  expect_true(all(result$estimate[result$estimate > 0] == 1L))
  expect_true(all(result$estimate[result$estimate < 0] == -1L))
})

test_that("prepare_exons_by_partition returns required metadata columns", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up   <- gr[gr$estimate > 0]
  gr_down <- gr[gr$estimate < 0]

  result <- prepare_exons_by_partition(gr_up, gr_down, verbose = FALSE)

  expect_true(all(
    c("tx_id", "gene_id", "exon_rank", "estimate") %in%
      names(GenomicRanges::mcols(result))
  ))
})

test_that("prepare_exons_by_partition fills missing extra columns with NA", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up   <- gr[gr$estimate > 0]
  gr_down <- gr[gr$estimate < 0]

  # add an extra column to gr_up only
  gr_up$extra_col <- "present"

  result <- prepare_exons_by_partition(gr_up, gr_down, verbose = FALSE)

  expect_true("extra_col" %in% names(GenomicRanges::mcols(result)))
  expect_true(any(is.na(result$extra_col)))
})

test_that("prepare_exons_by_partition passes through preprocess without error", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up   <- gr[gr$estimate > 0]
  gr_down <- gr[gr$estimate < 0]

  result <- prepare_exons_by_partition(gr_up, gr_down) |>
    preprocess(coef_col = "estimate")

  expect_s4_class(result, "GRanges")
  expect_true(isTRUE(S4Vectors::metadata(result)$splicelogic_preprocessed))
})

test_that("prepare_exons_by_partition mismatched input types throw an error", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up <- gr[gr$estimate > 0]

  expect_error(
    prepare_exons_by_partition(gr_up, c("ENST001", "ENST002"), verbose = FALSE),
    "both be GRanges objects or both be character vectors"
  )
})

test_that("prepare_exons_by_partition invalid input type throws an error", {
  expect_error(
    prepare_exons_by_partition(1:5, 6:10, verbose = FALSE),
    "both be GRanges objects or both be character vectors"
  )
})

test_that("prepare_exons_by_partition tx ID path without txdb throws an error", {
  expect_error(
    prepare_exons_by_partition(c("ENST001"), c("ENST002"), verbose = FALSE),
    "'txdb' is required"
  )
})

test_that("prepare_exons_by_partition with missing required columns throws an error", {
  gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
  gr_up   <- gr[gr$estimate > 0]
  gr_down <- gr[gr$estimate < 0]

  # remove a required column
  GenomicRanges::mcols(gr_up)$tx_id <- NULL

  expect_error(
    prepare_exons_by_partition(gr_up, gr_down, verbose = FALSE),
    "Missing required metadata columns"
  )
})

# --- tx ID path (requires network + AnnotationHub) ---
test_that("tx ID path returns preprocessed GRanges", {
  skip_if_not_installed("AnnotationHub")
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("AnnotationDbi")
  skip_if_offline()

  ah <- suppressMessages(AnnotationHub::AnnotationHub())
  txdb <- suppressMessages(ah[["AH84134"]])

  txps <- suppressMessages(AnnotationDbi::select(
    txdb,
    keys    = AnnotationDbi::keys(txdb, "TXID"),
    columns = c("TXNAME", "GENEID"),
    keytype = "TXID"
  )) |>
    tibble::as_tibble() |>
    dplyr::select(tx_id = TXNAME, gene_id = GENEID) |>
    dplyr::filter(!is.na(gene_id))

  test_txps <- txps |>
    dplyr::group_by(gene_id) |>
    dplyr::filter(dplyr::n() >= 4) |>
    dplyr::slice_head(n = 4) |>
    dplyr::ungroup()

  txid_up   <- test_txps$tx_id[1:2]
  txid_down <- test_txps$tx_id[3:4]

  result <- prepare_exons_by_partition(txid_up, txid_down, txdb = txdb, verbose = FALSE) |>
    preprocess(coef_col = "estimate")

  expect_s4_class(result, "GRanges")
  expect_true(isTRUE(S4Vectors::metadata(result)$splicelogic_preprocessed))
  expect_true(all(
    c("tx_id", "gene_id", "exon_rank", "estimate") %in%
      names(GenomicRanges::mcols(result))
  ))
  expect_true(any(result$estimate == 1L))
  expect_true(any(result$estimate == -1L))
})
