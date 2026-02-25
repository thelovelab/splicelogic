# Function to calculate retained introns given a GRanges object

Function to calculate retained introns given a GRanges object

## Usage

``` r
calc_retained_introns(gr)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coef'.

## Value

A GRanges object with an additional 'event' metadata column indicating
retained introns.
