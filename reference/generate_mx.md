# Generate mutually exclusive exon events in a GRanges object

Generate mutually exclusive exon events in a GRanges object

## Usage

``` r
generate_mx(gr, n_mx = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_mx:

  Number of mutually exclusive exon events to generate

## Value

A GRanges object with mutually exclusive exon events introduced
