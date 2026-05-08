# Preprocess input GRanges object for splicing event calculation

This function checks that the input is a valid GRanges object with
required metadata columns, then adds a unique key, the number of exons
per transcript, and an 'internal' flag for each exon.

## Usage

``` r
preprocess(gr, coef_col, method_string = NULL, additional_columns = NULL)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', 'coef'.

- coef_col:

  The name of the metadata column indicating upregulated (+1) and
  downregulated (-1) exons.

- method_string:

  The Differential Transcript Usage (DTU) method used to obtain the
  coef_col, for annotation purposes (optional).

- additional_columns:

  A character vector of metadata column names to record for downstream
  use. Stored in `metadata(result)$additional_columns` (optional).

## Value

A GRanges object with added 'key', 'nexons', and 'internal' columns.

## Examples

``` r

# create mock data and run preprocessing
gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4) |>
  preprocess(coef_col = "estimate", method_string = "mock_method")
```
