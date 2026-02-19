# Calculate skipped exons from a GRanges object

Calculate skipped exons from a GRanges object

## Usage

``` r
calc_skipped_exons(gr, coef_col, type = c("over", "in", "boundary"))
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns.

- coef_col:

  The name of the metadata column indicating upregulated (+1) and
  downregulated (-1) exons.

- type:

  The type of overlap to consider when identifying skipped exons.

## Value

A GRanges object with an additional 'event' metadata column indicating
skipped exons.
