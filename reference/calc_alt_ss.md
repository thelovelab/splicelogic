# Function to calculate alternative splice sites

Function to calculate alternative splice sites

## Usage

``` r
calc_alt_ss(gr, by_start = TRUE)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'coef'. Must be preprocessed with preprocess_input().

- by_start:

  If TRUE, detects a5ss (same exon start, different end). If FALSE,
  detects a3ss (same end, different start).

## Value

A GRanges object with annotated events for alternative splice sites.
