# Find splice events from annotated exons

Functions to find different types of alternative splicing events from
preprocessed GRanges exon data. Events include skipped exon (se),
included exon (ie), mutatualy exclusive exons (mxe), retained intron
(ri), and alternative 5' and 3' splice sites (a5ss / a3ss).

## Usage

``` r
find_se(gr, type = c("boundary", "over", "in"), inverse = FALSE)

find_ie(gr, type = c("boundary", "over", "in"))

find_mxe(gr, type = c("boundary", "in", "over"))

find_ri(gr)

find_a5ss(gr)

find_a3ss(gr)

find_all_events(gr, type = c("boundary", "over", "in"), verbose = TRUE)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess().

- type:

  The type of overlap to consider when identifying events.

- inverse:

  If TRUE, identifies included exons instead of skipped exons.

- verbose:

  If TRUE, prints progress messages. Default TRUE.

## Value

A GRanges object with an additional column `event` indicating:

`find_se()`: skipped exons

`find_ie()`: included exons

`find_mxe()`: mutually exclusive exons

`find_ri()`: retained introns

`find_a5ss()`: alternative 5' splice sites

`find_a3ss()`: : alternative 3' splice sites

`find_all_events()`: all detected events

## Examples

``` r

# make some mock data and run the function
gr <- create_mock_data(n_genes = 2, n_tx = 4, n_exons = 4) |>
  preprocess(coef_col = "estimate") |>
  generate_skipped_exons(1)

# this should find the skipped exon events we generated
find_se(gr, type = "boundary")
#> GRanges object with 1 range and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr21     11-15      + |         1         1         2 -0.525766
#>             event  tx_event
#>       <character> <numeric>
#>   [1]          se         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

find_ie(gr, type = "boundary")
#> GRanges object with 0 ranges and 0 metadata columns:
#>    seqnames    ranges strand
#>       <Rle> <IRanges>  <Rle>
#>   -------
#>   seqinfo: no sequences


# detect mutually exclusive exons
gr_mx <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_mx(1)

find_mxe(gr_mx, type = "boundary")
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr20     11-15      + |         1         3         2 -0.012726
#>   [2]    chr20     21-25      + |         1         2         2  0.285811
#>             event  tx_event
#>       <character> <numeric>
#>   [1]         mxe         2
#>   [2]         mxe         3
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect retained introns
gr_ri <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_retained_introns(1)

find_ri(gr_ri)
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr10   111-125      + |         2         7         2  0.840011
#>   [2]    chr10   111-125      + |         2         7         2  0.840011
#>             event  tx_event
#>       <character> <numeric>
#>   [1]          ri         5
#>   [2]          ri         8
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 5' splice sites
gr_a5 <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a5ss(1)

find_a5ss(gr_a5)
#> GRanges object with 1 range and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr22   121-123      + |         2         6         3  0.128291
#>             event  tx_event
#>       <character> <numeric>
#>   [1]        a5ss         5
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 3' splice sites
gr_a3 <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a3ss(1)
find_a3ss(gr_a3)
#> GRanges object with 2 ranges and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr19     23-25      + |         1         2         3    0.2278
#>   [2]    chr19     23-25      + |         1         2         3    0.2278
#>             event  tx_event
#>       <character> <numeric>
#>   [1]        a3ss         1
#>   [2]        a3ss         4
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect all event types at once
gr_all <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_skipped_exons(1)

find_all_events(gr_all, type = "boundary", verbose = FALSE)
#> GRanges object with 1 range and 6 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr13     21-25      + |         1         1         3 -0.285458
#>             event  tx_event
#>       <character> <numeric>
#>   [1]          se         2
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
