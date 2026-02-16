# Generate retained intron events in a GRanges object

Generate retained intron events in a GRanges object

## Usage

``` r
generate_retained_introns(gr, n_ri = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_ri:

  Number of retained intron events to generate

## Value

A GRanges object with retained intron events introduced
