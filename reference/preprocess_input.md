# Preprocess input GRanges object for splicing event calculation

This function checks that the input is a valid GRanges object with
required metadata columns, then adds a unique key, the number of exons
per transcript, and an 'internal' flag for each exon. It also
initializes an 'event' column for downstream splicing event annotation.

## Usage

``` r
preprocess_input(gr, coef_col, method_string = NULL)
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

## Value

A GRanges object with added 'key', 'nexons', 'internal', and 'event'
columns.
