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


# detect mutually exclusive exons
gr_mx <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess_input(coef_col = "coefs") |>
  generate_mx(1)
calc_mx_exons(gr_mx, type = "boundary")
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank      coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>   [1]    chr12   111-115      + |         2         8         2 -0.0419509
#>   [2]    chr12   121-125      + |         2         7         2  0.2729312
#>                    event  tx_event
#>              <character> <numeric>
#>   [1] mutually_exclusive         7
#>   [2] mutually_exclusive         8
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect retained introns
gr_ri <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess_input(coef_col = "coefs") |>
  generate_retained_introns(1)
calc_retained_introns(gr_ri)
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr22     11-25      + |         1         4         2  0.426795
#>   [2]    chr22     11-25      + |         1         4         2  0.426795
#>                 event  tx_event
#>           <character> <numeric>
#>   [1] retained_intron         1
#>   [2] retained_intron         3
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 5' splice sites
gr_a5 <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess_input(coef_col = "coefs") |>
  generate_a5ss(1)
calc_a5ss(gr_a5)
#> GRanges object with 3 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]     chr3     21-27      + |         1         2         3  0.836454
#>   [2]     chr3     21-27      + |         1         2         3  0.836454
#>   [3]     chr3     21-27      + |         1         2         3  0.836454
#>             event  tx_event
#>       <character> <numeric>
#>   [1]        a5ss         1
#>   [2]        a5ss         3
#>   [3]        a5ss         4
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 3' splice sites
gr_a3 <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess_input(coef_col = "coefs") |>
  generate_a3ss(1)
calc_a3ss(gr_a3)
#> GRanges object with 1 range and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr22   119-125      + |         2         6         3  0.128291
#>             event  tx_event
#>       <character> <numeric>
#>   [1]        a3ss         5
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect all event types at once
gr_all <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess_input(coef_col = "coefs") |>
  generate_skipped_exons(1)
calc_all_events(gr_all, type = "boundary", verbose = FALSE)
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr19     21-25      + |         1         1         3 -0.394266
#>   [2]    chr19     21-25      + |         1         4         3 -0.510564
#>              event  tx_event
#>        <character> <numeric>
#>   [1] skipped_exon         2
#>   [2] skipped_exon         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
