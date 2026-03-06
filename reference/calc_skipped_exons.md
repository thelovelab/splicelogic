# Calculate skipped exons from a GRanges object

Calculate skipped exons from a GRanges object

## Usage

``` r
calc_skipped_exons(gr, type = c("boundary", "over", "in"), inverse = FALSE)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess_input().

- type:

  The type of overlap to consider when identifying skipped exons.

- inverse:

  If TRUE, identifies included exons instead of skipped exons.

## Value

A GRanges object with an additional 'event' metadata column indicating
skipped exons.
