# Generate mock splice events in a GRanges object

Functions to introduce specific types of alternative splicing events
into mock GRanges data for testing purposes.

## Usage

``` r
generate_skipped_exons_restricted(gr, n_se = 1)

generate_skipped_exons(gr, n_se = 1)

generate_mx(gr, n_mx = 1)

generate_retained_introns(gr, n_ri = 1)

generate_a5ss(gr, n_a5ss = 1)

generate_a3ss(gr, n_a3ss = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_se:

  Number of skipped exon events to generate

- n_mx:

  Number of mutually exclusive exon events to generate

- n_ri:

  Number of retained intron events to generate

- n_a5ss:

  Number of alternative 5' splice site events to generate

- n_a3ss:

  Number of alternative 3' splice site events to generate

## Value

`generate_skipped_exons_restricted()`: A GRanges object with skipped
exon events introduced

`generate_skipped_exons()`: A GRanges object with skipped exon events
introduced

`generate_mx()`: A GRanges object with mutually exclusive exon events
introduced

`generate_retained_introns()`: A GRanges object with retained intron
events introduced

`generate_a5ss()`: A GRanges object with alternative 5' splice site
events introduced

`generate_a3ss()`: A GRanges object with alternative 3' splice site
events introduced
