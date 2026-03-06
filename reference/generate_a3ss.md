# Generate alternative 3' splice site events in a GRanges object

Generate alternative 3' splice site events in a GRanges object

## Usage

``` r
generate_a3ss(gr, n_a3ss = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_a3ss:

  Number of alternative 3' splice site events to generate

## Value

A GRanges object with alternative 3' splice site events introduced
