# Calculate alternative 5' splice sites from a GRanges object

Calculate alternative 5' splice sites from a GRanges object

## Usage

``` r
calc_a5ss(gr)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coef'. Must be preprocessed with preprocess_input().

## Value

A GRanges object with an additional 'event' metadata column indicating
a5ss events.
