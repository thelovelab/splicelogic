# Function to calcualte 5' and 3' alternative splice sites given a GRanges object

Function to calcualte 5' and 3' alternative splice sites given a GRanges
object

## Usage

``` r
calc_a3ss_a5ss(gr, coef_col)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coef'.

## Value

A GRanges object with an additional 'event' metadata column indicating
retained introns.
