# Basic function to visualize easily a gr object with transcript coefficients and events.

Basic function to visualize easily a gr object with transcript
coefficients and events.

## Usage

``` r
slPltRanges(gr)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'gene_id', 'tx_id', 'coefs',
  and optionally 'event'.

## Value

A ggplot object visualizing the transcripts and events for each gene.
