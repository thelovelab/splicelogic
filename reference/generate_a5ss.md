# Generate alternative 5' splice site events in a GRanges object

Generate alternative 5' splice site events in a GRanges object

## Usage

``` r
generate_a5ss(gr, n_a5ss = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_a5ss:

  Number of alternative 5' splice site events to generate

## Value

A GRanges object with alternative 5' splice site events introduced
