# Find splice events from annotated exons

Functions to find different types of alternative splicing events from
preprocessed GRanges exon data. Events include skipped exon (se),
included exon (ie), consecutive skipped exons (cse), consecutive
included exons (cie), mutatualy exclusive exons (mxe), retained intron
(ri), spliced intron (si), alternative 5' and 3' splice sites (a5ss /
a3ss), and alternative transcription start and end sites (atss / ates).

`find_ie()`, `find_cie()` and `find_si()` are the reverse of
`find_se()`, `find_cse()` and `find_ri()`: the same comparison with the
up and down transcripts swapped. A spliced intron is an intron of an up
transcript that a down transcript retains, i.e. an intron that the up
transcript splices out.

`find_se()` reports an exon skipped on its own: both of its neighbours
are kept by the partner transcript. Two or more consecutive exons
skipped together are reported by `find_cse()` instead, one row per exon,
when the exons on either side of the run match rank-adjacent exons of
the partner and no exon of the run overlaps the partner. The exons of
one run are the consecutive `exon_rank` values within a `tx_id` /
`event_tx_id` pair. A run that reaches the first or last exon is not a
cse; `find_atss()` / `find_ates()` cover it.

`find_a5ss()` / `find_a3ss()` do not require either boundary of an exon
to be shared with its partner, so an exon that moved both boundaries is
reported as an a5ss *and* an a3ss. A first exon has no acceptor and a
last exon has no donor, so those boundaries never produce an a3ss /
a5ss, whether the first or last exon belongs to the up or to the down
transcript of the pair; `find_atss()` / `find_ates()` report them
instead by comparing where each transcript begins and ends. The other
boundary of a first or last exon is a real splice site and is still
compared. `exon_rank` runs in transcription order, so this holds on both
strands.

When an exon overlaps two or more exons of its partner (a retained or
spliced intron), the boundaries facing the intron are left to
`find_ri()` / `find_si()`; only the outer boundaries of the group are
compared. An exon that retains an intron *and* uses a different outer
splice site is therefore reported as both an ri (or si) and an a5ss /
a3ss.

## Usage

``` r
find_se(gr, type = c("over", "boundary", "identical"), inverse = FALSE, ...)

find_skipped_exons(
  gr,
  type = c("over", "boundary", "identical"),
  inverse = FALSE,
  ...
)

find_ie(gr, type = c("over", "boundary", "identical"), ...)

find_included_exons(gr, type = c("over", "boundary", "identical"), ...)

find_cse(gr, type = c("over", "boundary", "identical"), inverse = FALSE, ...)

find_consecutive_skipped_exons(
  gr,
  type = c("over", "boundary", "identical"),
  inverse = FALSE,
  ...
)

find_cie(gr, type = c("over", "boundary", "identical"), ...)

find_consecutive_included_exons(
  gr,
  type = c("over", "boundary", "identical"),
  ...
)

find_mxe(gr, type = c("over", "boundary", "identical"), ...)

find_mutually_exclusive_exons(
  gr,
  type = c("over", "boundary", "identical"),
  ...
)

find_ri(gr, inverse = FALSE)

find_retained_introns(gr, inverse = FALSE)

find_si(gr)

find_spliced_introns(gr)

find_a5ss(gr)

find_alternative_5_prime_splice_sites(gr)

find_a3ss(gr)

find_alternative_3_prime_splice_sites(gr)

find_atss(gr)

find_alternative_transcription_start_sites(gr)

find_ates(gr)

find_alternative_transcription_end_sites(gr)

find_all_events(
  gr,
  type = c("over", "boundary", "identical"),
  verbose = TRUE,
  ...
)
```

## Arguments

- gr:

  A GRanges object with exon annotations, including 'tx_id', 'exon', and
  'coef_col' metadata columns and preprocessed with preprocess().

- type:

  How an exon flanking a candidate is matched to an exon of the partner
  transcript. One of:

  `"over"`

  : (default) any overlap, with no coordinate in common required. The
    most permissive setting.

  `"boundary"`

  : the two exons overlap *and* share the inner boundary of the flank,
    the splice site facing the candidate: the donor of the upstream
    flank and the acceptor of the downstream one. On `"+"` that is the
    end of the left exon and the start of the right one; on `"-"` the
    reverse. The outer boundaries may differ.

  `"identical"`

  : the two exons are identical (same start and same end). The strictest
    setting.

- inverse:

  If TRUE, swaps the up and down transcripts: `find_se()` then
  identifies included exons instead of skipped exons, `find_cse()`
  consecutive included exons, and `find_ri()` spliced introns instead of
  retained introns.

- ...:

  Further arguments passed to
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  when matching flanking exons, e.g. `minoverlap = 5L` to require at
  least 5 bp of overlap. Do not pass `select`, which changes the return
  type. `maxgap` can also be used, to let `"over"` and `"boundary"`
  match exons that do not overlap at all.

- verbose:

  If TRUE, prints progress messages. Default TRUE.

## Value

A GRanges object with the detected exon ranges and the following
additional metadata columns:

- `event_type`:

  The type of splicing event detected (e.g. `"se"`, `"ie"`, `"cse"`,
  `"cie"`, `"mxe"`, `"ri"`, `"si"`, `"a5ss"`, `"a3ss"`, `"atss"`,
  `"ates"`).

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

`find_cse()`: consecutive skipped exons, i.e. runs of two or more
adjacent exons skipped together, one row per exon

`find_cie()`: consecutive included exons

`find_mxe()`: mutually exclusive exons

`find_ri()`: retained introns

`find_si()`: spliced introns

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

# included exons: the same edit made to a down transcript
gr_ie <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_ie(n_events = 1)

find_ie(gr_ie, type = "boundary")
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]     chr1     11-15      + |         1         2         2  0.485274
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          ie           1      -0.980135
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# an up transcript skips a run of two consecutive exons
gr_cse <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
) |>
  preprocess(coef_col = "estimate") |>
  generate_cse(n_events = 1, n_skip = 2)

find_cse(gr_cse)
#> GRanges object with 6 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>   [1]     chr7     21-25      + |         1         1         3 -0.5858222
#>   [2]     chr7     31-35      + |         1         1         4 -0.5858222
#>   [3]     chr7     21-25      + |         1         3         3 -0.9440879
#>   [4]     chr7     31-35      + |         1         3         4 -0.9440879
#>   [5]     chr7     21-25      + |         1         4         3 -0.0612314
#>   [6]     chr7     31-35      + |         1         4         4 -0.0612314
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE         cse           2       0.889367
#>   [2]      TRUE         cse           2       0.889367
#>   [3]      TRUE         cse           2       0.889367
#>   [4]      TRUE         cse           2       0.889367
#>   [5]      TRUE         cse           2       0.889367
#>   [6]      TRUE         cse           2       0.889367
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# find_se() only reports an exon skipped on its own, so finds nothing
find_se(gr_cse)
#> GRanges object with 0 ranges and 0 metadata columns:
#>    seqnames    ranges strand
#>       <Rle> <IRanges>  <Rle>
#>   -------
#>   seqinfo: no sequences


# an up transcript includes a run of two exons that a down one skips
gr_cie <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
) |>
  preprocess(coef_col = "estimate") |>
  generate_cie(n_events = 1, n_skip = 2)

find_cie(gr_cie)
#> GRanges object with 6 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr11   111-115      + |         2         6         2  0.597618
#>   [2]    chr11   121-125      + |         2         6         3  0.597618
#>   [3]    chr11   111-115      + |         2         7         2  0.326413
#>   [4]    chr11   121-125      + |         2         7         3  0.326413
#>   [5]    chr11   111-115      + |         2         8         2  0.713150
#>   [6]    chr11   121-125      + |         2         8         3  0.713150
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE         cie           5      -0.757729
#>   [2]      TRUE         cie           5      -0.757729
#>   [3]      TRUE         cie           5      -0.757729
#>   [4]      TRUE         cie           5      -0.757729
#>   [5]      TRUE         cie           5      -0.757729
#>   [6]      TRUE         cie           5      -0.757729
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


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
#>   [1]    chr19     11-15      + |         1         4         2 -0.495589
#>   [2]    chr19     21-25      + |         1         3         2  0.136705
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE         mxe           3       0.136705
#>   [2]      TRUE         mxe           4      -0.495589
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect retained introns
gr_ri <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_ri(n_events = 1)

find_ri(gr_ri)
#> GRanges object with 3 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr12     11-25      + |         1         2         2  0.355806
#>   [2]    chr12     11-25      + |         1         2         2  0.355806
#>   [3]    chr12     11-25      + |         1         2         2  0.355806
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          ri           1     -0.6813640
#>   [2]      TRUE          ri           3     -0.7873061
#>   [3]      TRUE          ri           4     -0.0303944
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# spliced introns: the same edit made to a down transcript, so the up
# transcripts splice out the intron it retains
gr_si <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_si(n_events = 1)

find_si(gr_si)
#> GRanges object with 3 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr20   101-115      + |         2         5         2 -0.350577
#>   [2]    chr20   101-115      + |         2         5         2 -0.350577
#>   [3]    chr20   101-115      + |         2         5         2 -0.350577
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          si           6       0.667885
#>   [2]      TRUE          si           7       0.518317
#>   [3]      TRUE          si           8       0.672151
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
#>   [1]     chr7     11-13      + |         1         3         2  0.346896
#>   [2]     chr7     11-13      + |         1         3         2  0.346896
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a5ss           1      -0.462591
#>   [2]      TRUE        a5ss           4      -0.257460
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
#>   [1]    chr20   113-115      + |         2         6         3 0.0455507
#>   [2]    chr20   113-115      + |         2         6         3 0.0455507
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        a3ss           5       -0.79801
#>   [2]      TRUE        a3ss           8       -0.85807
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative transcription start sites
gr_atss <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_atss(n_events = 1)

find_atss(gr_atss)
#> GRanges object with 2 ranges and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr21       3-5      + |         1         4         1 0.0971829
#>   [2]    chr21       3-5      + |         1         4         1 0.0971829
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        atss           1      -0.985379
#>   [2]      TRUE        atss           3      -0.592199
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


# detect alternative transcription end sites
gr_ates <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
) |>
  preprocess(coef_col = "estimate") |>
  generate_ates(n_events = 1)

find_ates(gr_ates)
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr16   121-123      + |         2         8         4  0.283883
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE        ates           5     -0.0835998
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
#>       seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>   [1]     chr2   101-105      + |         2         5         2 -0.0390255
#>   [2]     chr2   101-105      + |         2         8         2 -0.7484215
#>       sim_event  event_type event_tx_id event_estimate
#>       <logical> <character>   <numeric>      <numeric>
#>   [1]      TRUE          se           6       0.857323
#>   [2]      TRUE          se           6       0.857323
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
