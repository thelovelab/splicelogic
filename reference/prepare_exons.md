# Prepare exon ranges from a TxDb and DTU results table

Extracts exon ranges from a TxDb object, merges them with differential
transcript usage (DTU) results, and returns a flat GRanges ready for
[`preprocess_input`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md).

## Usage

``` r
prepare_exons(
  txdb,
  dtu_table,
  coef_col,
  tx_id_col = "tx_id",
  gene_id_col = "gene_id",
  verbose = TRUE
)
```

## Arguments

- txdb:

  A `TxDb` object (from GenomicFeatures).

- dtu_table:

  A data.frame or tibble with DTU results. Must contain columns for
  transcript ID, gene ID, and a coefficient.

- coef_col:

  Column name in `dtu_table` with the coefficient / effect size values.

- tx_id_col:

  Column name in `dtu_table` with transcript IDs matching the TxDb
  transcript names. Default `"tx_id"`.

- gene_id_col:

  Column name in `dtu_table` with gene IDs. Default `"gene_id"`.

- verbose:

  Whether to print progress messages. Default `TRUE`.

## Value

A GRanges object with metadata columns: `gene_id`, `tx_id`, `exon_rank`,
the coefficient column, and any additional columns from `dtu_table`.
