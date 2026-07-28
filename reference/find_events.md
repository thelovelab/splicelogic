# Find splice events from annotated exons

Functions to find different types of alternative splicing events from
preprocessed GRanges exon data. Events include skipped exon (se),
included exon (ie), mutatualy exclusive exons (mxe), retained intron
(ri), and alternative 5' and 3' splice sites (a5ss / a3ss).

## Usage

``` r
find_se(gr, type = c("boundary", "over", "in"), inverse = FALSE)

find_skipped_exons(gr, type = c("boundary", "over", "in"), inverse = FALSE)

find_ie(gr, type = c("boundary", "over", "in"))

find_included_exons(gr, type = c("boundary", "over", "in"))

find_mxe(gr, type = c("boundary", "in", "over"))

find_mutually_exclusive_exons(gr, type = c("boundary", "in", "over"))

find_ri(gr)

find_retained_introns(gr)

find_a5ss(gr)

find_alternative_5_prime_splice_sites(gr)

find_a3ss(gr)

find_alternative_3_prime_splice_sites(gr)

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

A GRanges object with the detected exon ranges and the following
additional metadata columns:

- `event_type`:

  The type of splicing event detected (e.g. `"se"`, `"ie"`, `"mxe"`,
  `"ri"`, `"a5ss"`, `"a3ss"`).

- `event_tx_id`:

  Transcript ID of the paired transcript involved in the event.

- `event_estimate`:

  DTU coefficient of the paired transcript.

- `event_<col>`:

  One column per name in `metadata(gr)$additional_columns`, prefixed
  with `event_`, carrying the corresponding value from the paired
  transcript.

`find_se()`: skipped exons

`find_ie()`: included exons

`find_mxe()`: mutually exclusive exons

`find_ri()`: retained introns

`find_a5ss()`: alternative 5' splice sites

`find_a3ss()`: alternative 3' splice sites

`find_all_events()`: all detected events

## Examples

``` r

# make some mock data and run the function
gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4) |>
  preprocess(coef_col = "estimate") |>
  generate_se(n_events = 1)

# this should find the skipped exon events we generated
find_se(gr, type = "boundary")
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr21   111-115      + |         2         5         3 -0.821465
#>   [2]    chr21   111-115      + |         2         7         3 -0.398338
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          se           8       0.272931
#>   [2]      TRUE          se           8       0.272931
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
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_mxe(n_events = 1)

find_mxe(gr_mx, type = "boundary")
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr22     11-15      + |         1         3         2 -0.591643
#>   [2]    chr22     21-25      + |         1         2         2  0.377796
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE         mxe           2       0.377796
#>   [2]      TRUE         mxe           3      -0.591643
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect retained introns
gr_ri <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_ri(n_events = 1)

find_ri(gr_ri)
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]     chr6   101-115      + |         2         6         2  0.409872
#>   [2]     chr6   101-115      + |         2         6         2  0.409872
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          ri           5      -0.162341
#>   [2]      TRUE          ri           7      -0.573655
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 5' splice sites
gr_a5 <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a5ss(n_events = 1)

find_a5ss(gr_a5)
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr14   101-103      + |         2         7         2  0.992694
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a5ss           5      -0.149921
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 3' splice sites
gr_a3 <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a3ss(n_events = 1)
find_a3ss(gr_a3)
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]     chr9   109-115      + |         2         6         3  0.572669
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a3ss           5      -0.765392
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect all event types at once
gr_all <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_se(n_events = 1)

find_all_events(gr_all, type = "boundary", verbose = FALSE)
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr16     11-15      + |         1         1         2  -0.14205
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          se           4       0.990246
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
