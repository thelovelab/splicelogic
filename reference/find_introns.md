# function to find introns given a GRanges object of exons

function to find introns given a GRanges object of exons

## Usage

``` r
find_introns(gr)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

## Value

A GRanges object with introns as ranges and metadata (tx_id, gene_id)
