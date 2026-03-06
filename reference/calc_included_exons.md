# Calculate included exons from a GRanges object

Calculate included exons from a GRanges object

## Usage

``` r
calc_included_exons(gr, type = c("boundary", "over", "in"))
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess_input().

- type:

  The type of overlap to consider when identifying included exons.

## Value

A GRanges object with an additional 'event' metadata column indicating
included exons.
