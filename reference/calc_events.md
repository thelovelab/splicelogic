# Calculate splice events from a GRanges object

Functions to detect different types of alternative splicing events from
preprocessed GRanges exon data.

## Usage

``` r
calc_skipped_exons(gr, type = c("boundary", "over", "in"), inverse = FALSE)

calc_included_exons(gr, type = c("boundary", "over", "in"))

calc_mx_exons(gr, type = c("boundary", "in", "over"))

calc_retained_introns(gr)

calc_alt_ss(gr, by_start = TRUE)

calc_a5ss(gr)

calc_a3ss(gr)

calc_all_events(gr, type = c("boundary", "over", "in"), verbose = TRUE)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess_input().

- type:

  The type of overlap to consider when identifying events.

- inverse:

  If TRUE, identifies included exons instead of skipped exons.

- by_start:

  If TRUE, detects a5ss (same exon start, different end). If FALSE,
  detects a3ss (same end, different start).

- verbose:

  If TRUE, prints progress messages. Default TRUE.

## Value

`calc_skipped_exons()`: A GRanges object with an additional 'event'
metadata column indicating skipped exons.

`calc_included_exons()`: A GRanges object with an additional 'event'
metadata column indicating included exons.

`calc_mx_exons()`: A GRanges object with an additional 'event' metadata
column indicating mutually exclusive exons.

`calc_retained_introns()`: A GRanges object with an additional 'event'
metadata column indicating retained introns.

`calc_alt_ss()`: A GRanges object with annotated events for alternative
splice sites.

`calc_a5ss()`: A GRanges object with an additional 'event' metadata
column indicating a5ss events.

`calc_a3ss()`: A GRanges object with an additional 'event' metadata
column indicating a3ss events.

`calc_all_events()`: A GRanges object combining all detected events,
with an 'event' metadata column indicating the event type.

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
