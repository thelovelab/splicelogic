# Generate skipped exon events in a GRanges object making all \>0 coefs transcripts have a the same exon eventfor the same gene

Generate skipped exon events in a GRanges object making all \>0 coefs
transcripts have a the same exon eventfor the same gene

## Usage

``` r
generate_skipped_exons_restricted(gr, n_se = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coefs'.

- n_se:

  Number of skipped exon events to generate

## Value

A GRanges object with skipped exon events introduced
