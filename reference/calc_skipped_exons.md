# Calculate skipped exons from a GRanges object

Calculate skipped exons from a GRanges object

## Usage

``` r
calc_skipped_exons(gr, type = c("boundary", "over", "in"), inverse = FALSE)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess_input().

- type:

  The type of overlap to consider when identifying skipped exons.

- inverse:

  If TRUE, identifies included exons instead of skipped exons.

## Value

A GRanges object with an additional 'event' metadata column indicating
skipped exons.

## Examples

``` r

# make some mock data and run the function
gr <- create_mock_data(n_genes = 2, n_tx = 4, n_exons = 4) |>
  preprocess_input(coef_col = "coefs") |> 
  generate_skipped_exons(1)

# this should find the skipped exon events we generated
calc_skipped_exons(gr, type = "boundary")
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr13   111-115      + |         2         5         2 -0.682818
#>   [2]    chr13   111-115      + |         2         7         2 -0.420466
#>              event  tx_event
#>        <character> <numeric>
#>   [1] skipped_exon         8
#>   [2] skipped_exon         8
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

calc_included_exons(gr, type = "boundary")
#> GRanges object with 0 ranges and 0 metadata columns:
#>    seqnames    ranges strand
#>       <Rle> <IRanges>  <Rle>
#>   -------
#>   seqinfo: no sequences
```
