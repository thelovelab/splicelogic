# Find splice events from annotated exons

Functions to find different types of alternative splicing events from
preprocessed GRanges exon data. Events include skipped exon (se),
included exon (ie), mutatualy exclusive exons (mxe), retained intron
(ri), alternative 5' and 3' splice sites (a5ss / a3ss), and alternative
transcription start and end sites (atss / ates).

`find_a5ss()` / `find_a3ss()` do not require either boundary of an exon
to be shared with its partner, so an exon that moved both boundaries is
reported as an a5ss *and* an a3ss. A first exon has no acceptor and a
last exon has no donor, so those boundaries never produce an a3ss /
a5ss; `find_atss()` / `find_ates()` report them instead by comparing
where each transcript begins and ends. `exon_rank` runs in transcription
order, so this holds on both strands.

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

find_atss(gr)

find_alternative_transcription_start_sites(gr)

find_ates(gr)

find_alternative_transcription_end_sites(gr)

find_all_events(gr, type = c("boundary", "over", "in"), verbose = TRUE)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess().

- type:

  How an exon flanking a candidate is matched to an exon of the partner
  transcript. One of:

  `"boundary"`

  :   (default) the two exons overlap *and* share a start or an end
      coordinate, so they may still differ in length at the other end.

  `"over"`

  :   any overlap, with no coordinate in common required. The most
      permissive setting.

  `"in"`

  :   the two exons are identical (same start and same end). The
      strictest setting.

- inverse:

  If TRUE, identifies included exons instead of skipped exons.

- verbose:

  If TRUE, prints progress messages. Default TRUE.

## Value

A GRanges object with the detected exon ranges and the following
additional metadata columns:

- `event_type`:

  The type of splicing event detected (e.g. `"se"`, `"ie"`, `"mxe"`,
  `"ri"`, `"a5ss"`, `"a3ss"`, `"atss"`, `"ates"`).

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

`find_atss()`: alternative transcription start sites — first exons that
begin transcription at a different coordinate

`find_ates()`: alternative transcription end sites — last exons that end
transcription at a different coordinate

`find_all_events()`: all detected events

## Examples

``` r

# make some mock data and run the function
gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4) |>
  preprocess(coef_col = "estimate") |>
  generate_se(n_events = 1)

# this should find the skipped exon events we generated
find_se(gr, type = "boundary")
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr11   101-105      + |         2         5         2 -0.687911
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          se           8        0.53935
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
#>   [1]     chr1     11-15      + |         1         1         2 -0.980135
#>   [2]     chr1     21-25      + |         1         2         2  0.485274
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE         mxe           2       0.485274
#>   [2]      TRUE         mxe           1      -0.980135
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
#>   [1]     chr1     11-25      + |         1         2         2 0.0213647
#>   [2]     chr1     11-25      + |         1         2         2 0.0213647
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          ri           1      -0.783753
#>   [2]      TRUE          ri           4      -0.944088
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 5' splice sites
gr_a5 <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a5ss(n_events = 1)

find_a5ss(gr_a5)
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr18   101-107      + |         2         6         2  0.615734
#>   [2]    chr18   101-107      + |         2         6         2  0.615734
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a5ss           5      -0.082719
#>   [2]      TRUE        a5ss           8      -0.761029
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative 3' splice sites
gr_a3 <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_a3ss(n_events = 1)
find_a3ss(gr_a3)
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr11     13-15      + |         1         3         2  0.399659
#>   [2]    chr11     13-15      + |         1         3         2  0.399659
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a3ss           1      -0.437331
#>   [2]      TRUE        a3ss           4      -0.559999
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr_tss <- data.frame(
  seqnames = "chr1",
  start = c(1, 21, 41, 5, 21, 41),
  end = c(10, 30, 50, 10, 30, 50),
  strand = "+",
  gene_id = "g1",
  tx_id = rep(c("down", "up"), each = 3),
  exon_rank = rep(seq_len(3), 2),
  estimate = rep(c(-1, 1), each = 3)
) |>
  plyranges::as_granges() |>
  preprocess(coef_col = "estimate")

find_atss(gr_tss)
#> GRanges object with 1 range and 7 metadata columns:
#>       seqnames    ranges strand |     gene_id       tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <character> <character> <integer> <numeric>
#>   [1]     chr1      5-10      + |          g1          up         1         1
#>        event_type event_tx_id event_estimate
#>       <character> <character>      <numeric>
#>   [1]        atss        down             -1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect all event types at once
gr_all <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_se(n_events = 1)

find_all_events(gr_all, type = "boundary", verbose = FALSE)
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]     chr8     21-25      + |         1         1         3 -0.155684
#>   [2]     chr8     21-25      + |         1         3         3 -0.297502
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          se           4       0.476183
#>   [2]      TRUE          se           4       0.476183
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
