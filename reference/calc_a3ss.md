# Calculate alternative 3' splice sites from a GRanges object

Calculate alternative 3' splice sites from a GRanges object

## Usage

``` r
calc_a3ss(gr)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coef'. Must be preprocessed with preprocess_input().

## Value

A GRanges object with an additional 'event' metadata column indicating
a3ss events.
