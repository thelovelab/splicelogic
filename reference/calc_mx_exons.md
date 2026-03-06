# Calculate mutually exclusive exons from a GRanges object

Calculate mutually exclusive exons from a GRanges object

## Usage

``` r
calc_mx_exons(gr, type = c("boundary", "in", "over"))
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess_input().

- type:

  The type of overlap to consider when identifying mutually exclusive
  exons.

## Value

A GRanges object with an additional 'event' metadata column indicating
mutually exclusive exons.
