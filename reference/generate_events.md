# Generate mock splice events in a GRanges object

Functions to introduce specific types of alternative splicing events
into mock GRanges data for testing purposes.

Each generator edits the transcripts on one side and marks, in a logical
`sim_event` column, the exons the matching finder should report.
`generate_se()`, `generate_cse()`, `generate_ri()`, `generate_a5ss()`,
`generate_a3ss()`, `generate_atss()` and `generate_ates()` edit up
(positive `estimate`) transcripts. `generate_ie()`, `generate_cie()` and
`generate_si()` are their mirrors: the same edit made to a down
(negative `estimate`) transcript. `generate_mxe()` edits one transcript
on each side.

## Usage

``` r
generate_se(gr, n_events = 1)

generate_ie(gr, n_events = 1)

generate_cse(gr, n_events = 1, n_skip = 2)

generate_cie(gr, n_events = 1, n_skip = 2)

generate_mxe(gr, n_events = 1)

generate_ri(gr, n_events = 1)

generate_si(gr, n_events = 1)

generate_a5ss(gr, n_events = 1)

generate_a3ss(gr, n_events = 1)

generate_atss(gr, n_events = 1, mode = c("shift", "drop"))

generate_ates(gr, n_events = 1, mode = c("shift", "drop"))
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'estimate'.

- n_events:

  Number of events to generate, i.e. the number of genes that get one.
  Every generator puts at most one event in a gene, so two events never
  edit the same transcripts. When fewer genes than `n_events` have a
  candidate, it warns and makes one in each of those.

- n_skip:

  Number of consecutive exons skipped (or included) together, at
  least 2. Every exon of the run is internal, so a transcript needs at
  least `n_skip + 2` exons to host one.

- mode:

  How to move the site: `"shift"` moves the terminal exon's outer
  boundary (so it still overlaps its partner), `"drop"` removes the
  terminal exon and re-ranks (so the new terminal exon does not overlap
  its partner at all).

## Value

`generate_se()`: A GRanges object with skipped exon events introduced:
an internal exon is removed from one up transcript per gene

`generate_ie()`: A GRanges object with included exon events introduced:
an internal exon is removed from one down transcript per gene, so the up
transcripts include it

`generate_cse()`: A GRanges object with consecutive skipped exon events
introduced: a run of `n_skip` internal exons is removed from one up
transcript per gene

`generate_cie()`: A GRanges object with consecutive included exon events
introduced: a run of `n_skip` internal exons is removed from one down
transcript per gene, so the up transcripts include it

`generate_mxe()`: A GRanges object with mutually exclusive exon events
introduced: a down transcript loses internal exon k + 1 and an up
transcript of the same gene loses exon k, so each keeps the exon the
other lacks. One pair per gene.

`generate_ri()`: A GRanges object with retained intron events
introduced: exons 2 and 3 of one up transcript per gene are fused

`generate_si()`: A GRanges object with spliced intron events introduced:
exons 2 and 3 of one down transcript per gene are fused, so the up
transcripts splice out the intron it retains

`generate_a5ss()`: A GRanges object with alternative 5' splice site
events introduced: the donor of an internal exon is moved by 2 bp in one
up transcript per gene

`generate_a3ss()`: A GRanges object with alternative 3' splice site
events introduced: the acceptor of an internal exon is moved by 2 bp in
one up transcript per gene

`generate_atss()`: A GRanges object with alternative transcription start
site events introduced: the first exon of one up transcript per gene is
shortened at its outer end (`"shift"`) or removed (`"drop"`)

`generate_ates()`: A GRanges object with alternative transcription end
site events introduced: the last exon of one up transcript per gene is
shortened at its outer end (`"shift"`) or removed (`"drop"`)

## Examples

``` r

gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_se(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr7       1-5      + |         1         1         1  -0.35974
#>    [2]     chr7     11-15      + |         1         1         2  -0.35974
#>    [3]     chr7     21-25      + |         1         1         3  -0.35974
#>    [4]     chr7     31-35      + |         1         1         4  -0.35974
#>    [5]     chr7       1-5      + |         1         2         1   0.75178
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]     chr7   121-125      + |         2         7         4  0.105800
#>   [28]     chr7     91-95      + |         2         8         1 -0.561528
#>   [29]     chr7   101-105      + |         2         8         2 -0.561528
#>   [30]     chr7   111-115      + |         2         8         3 -0.561528
#>   [31]     chr7   121-125      + |         2         8         4 -0.561528
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE      TRUE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         3     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [27]         7-4         4     FALSE     FALSE
#>   [28]         8-1         4     FALSE     FALSE
#>   [29]         8-2         4      TRUE     FALSE
#>   [30]         8-3         4      TRUE     FALSE
#>   [31]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_ie(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr17       1-5      + |         1         1         1 -0.885963
#>    [2]    chr17     11-15      + |         1         1         2 -0.885963
#>    [3]    chr17     21-25      + |         1         1         3 -0.885963
#>    [4]    chr17     31-35      + |         1         1         4 -0.885963
#>    [5]    chr17       1-5      + |         1         2         1  0.534635
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr17   111-115      + |         2         7         3  0.884264
#>   [28]    chr17   121-125      + |         2         7         4  0.884264
#>   [29]    chr17     91-95      + |         2         8         1 -0.751253
#>   [30]    chr17   101-105      + |         2         8         2 -0.751253
#>   [31]    chr17   121-125      + |         2         8         3 -0.751253
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [27]         7-3         4      TRUE      TRUE
#>   [28]         7-4         4     FALSE     FALSE
#>   [29]         8-1         3     FALSE     FALSE
#>   [30]         8-2         3      TRUE     FALSE
#>   [31]         8-3         3     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
)
generate_cse(gr, n_events = 1)
#> GRanges object with 38 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.417802
#>    [2]    chr10     11-15      + |         1         1         2 -0.417802
#>    [3]    chr10     21-25      + |         1         1         3 -0.417802
#>    [4]    chr10     31-35      + |         1         1         4 -0.417802
#>    [5]    chr10     41-45      + |         1         1         5 -0.417802
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [34]    chr10   101-105      + |         2         8         1   0.78257
#>   [35]    chr10   111-115      + |         2         8         2   0.78257
#>   [36]    chr10   121-125      + |         2         8         3   0.78257
#>   [37]    chr10   131-135      + |         2         8         4   0.78257
#>   [38]    chr10   141-145      + |         2         8         5   0.78257
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         5     FALSE     FALSE
#>    [2]         1-2         5      TRUE     FALSE
#>    [3]         1-3         5      TRUE     FALSE
#>    [4]         1-4         5      TRUE     FALSE
#>    [5]         1-5         5     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [34]         8-1         5     FALSE     FALSE
#>   [35]         8-2         5      TRUE     FALSE
#>   [36]         8-3         5      TRUE     FALSE
#>   [37]         8-4         5      TRUE     FALSE
#>   [38]         8-5         5     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# a run of three needs at least five exons
generate_cse(gr, n_events = 1, n_skip = 3)
#> GRanges object with 37 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.417802
#>    [2]    chr10     11-15      + |         1         1         2 -0.417802
#>    [3]    chr10     21-25      + |         1         1         3 -0.417802
#>    [4]    chr10     31-35      + |         1         1         4 -0.417802
#>    [5]    chr10     41-45      + |         1         1         5 -0.417802
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [33]    chr10   101-105      + |         2         8         1   0.78257
#>   [34]    chr10   111-115      + |         2         8         2   0.78257
#>   [35]    chr10   121-125      + |         2         8         3   0.78257
#>   [36]    chr10   131-135      + |         2         8         4   0.78257
#>   [37]    chr10   141-145      + |         2         8         5   0.78257
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         5     FALSE     FALSE
#>    [2]         1-2         5      TRUE     FALSE
#>    [3]         1-3         5      TRUE     FALSE
#>    [4]         1-4         5      TRUE     FALSE
#>    [5]         1-5         5     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [33]         8-1         5     FALSE     FALSE
#>   [34]         8-2         5      TRUE     FALSE
#>   [35]         8-3         5      TRUE     FALSE
#>   [36]         8-4         5      TRUE     FALSE
#>   [37]         8-5         5     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 5
)
generate_cie(gr, n_events = 1)
#> GRanges object with 38 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr3       1-5      + |         1         1         1 -0.102875
#>    [2]     chr3     11-15      + |         1         1         2 -0.102875
#>    [3]     chr3     21-25      + |         1         1         3 -0.102875
#>    [4]     chr3     31-35      + |         1         1         4 -0.102875
#>    [5]     chr3     41-45      + |         1         1         5 -0.102875
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [34]     chr3   101-105      + |         2         8         1  0.843599
#>   [35]     chr3   111-115      + |         2         8         2  0.843599
#>   [36]     chr3   121-125      + |         2         8         3  0.843599
#>   [37]     chr3   131-135      + |         2         8         4  0.843599
#>   [38]     chr3   141-145      + |         2         8         5  0.843599
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         5     FALSE     FALSE
#>    [2]         1-2         5      TRUE     FALSE
#>    [3]         1-3         5      TRUE     FALSE
#>    [4]         1-4         5      TRUE     FALSE
#>    [5]         1-5         5     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [34]         8-1         5     FALSE     FALSE
#>   [35]         8-2         5      TRUE      TRUE
#>   [36]         8-3         5      TRUE      TRUE
#>   [37]         8-4         5      TRUE     FALSE
#>   [38]         8-5         5     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_mxe(gr, n_events = 1)
#> GRanges object with 30 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.745587
#>    [2]    chr10     11-15      + |         1         1         2 -0.745587
#>    [3]    chr10     21-25      + |         1         1         3 -0.745587
#>    [4]    chr10     31-35      + |         1         1         4 -0.745587
#>    [5]    chr10       1-5      + |         1         2         1  0.563896
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [26]    chr10   121-125      + |         2         7         4  0.788245
#>   [27]    chr10     91-95      + |         2         8         1 -0.592671
#>   [28]    chr10   101-105      + |         2         8         2 -0.592671
#>   [29]    chr10   111-115      + |         2         8         3 -0.592671
#>   [30]    chr10   121-125      + |         2         8         4 -0.592671
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [26]         7-4         4     FALSE     FALSE
#>   [27]         8-1         4     FALSE     FALSE
#>   [28]         8-2         4      TRUE     FALSE
#>   [29]         8-3         4      TRUE     FALSE
#>   [30]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_ri(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr11       1-5      + |         1         1         1 -0.4809627
#>    [2]    chr11     11-15      + |         1         1         2 -0.4809627
#>    [3]    chr11     21-25      + |         1         1         3 -0.4809627
#>    [4]    chr11     31-35      + |         1         1         4 -0.4809627
#>    [5]    chr11       1-5      + |         1         2         1  0.0838873
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [27]    chr11   111-115      + |         2         8         3   0.486036
#>   [28]    chr11   121-125      + |         2         8         4   0.486036
#>   [29]    chr11     91-95      + |         2         7         1   0.189638
#>   [30]    chr11   101-115      + |         2         7         2   0.189638
#>   [31]    chr11   121-125      + |         2         7         3   0.189638
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [27]         8-3         4      TRUE     FALSE
#>   [28]         8-4         4     FALSE     FALSE
#>   [29]         7-1         3     FALSE     FALSE
#>   [30]         7-2         3      TRUE      TRUE
#>   [31]         7-3         3     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_si(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr11       1-5      + |         1         1         1 -0.363362
#>    [2]    chr11     11-15      + |         1         1         2 -0.363362
#>    [3]    chr11     21-25      + |         1         1         3 -0.363362
#>    [4]    chr11     31-35      + |         1         1         4 -0.363362
#>    [5]    chr11       1-5      + |         1         2         1  0.678905
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr11   111-115      + |         2         8         3  0.655843
#>   [28]    chr11   121-125      + |         2         8         4  0.655843
#>   [29]    chr11       1-5      + |         1         3         1 -0.641931
#>   [30]    chr11     11-25      + |         1         3         2 -0.641931
#>   [31]    chr11     31-35      + |         1         3         3 -0.641931
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [27]         8-3         4      TRUE     FALSE
#>   [28]         8-4         4     FALSE     FALSE
#>   [29]         3-1         3     FALSE     FALSE
#>   [30]         3-2         3      TRUE      TRUE
#>   [31]         3-3         3     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_a5ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr16       1-5      + |         1         1         1 -0.273608
#>    [2]    chr16     11-15      + |         1         1         2 -0.273608
#>    [3]    chr16     21-25      + |         1         1         3 -0.273608
#>    [4]    chr16     31-35      + |         1         1         4 -0.273608
#>    [5]    chr16       1-5      + |         1         2         1  0.840244
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr16   121-125      + |         2         7         4  0.619392
#>   [29]    chr16     91-95      + |         2         8         1  0.118766
#>   [30]    chr16   101-105      + |         2         8         2  0.118766
#>   [31]    chr16   111-115      + |         2         8         3  0.118766
#>   [32]    chr16   121-125      + |         2         8         4  0.118766
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE     FALSE
#>   [29]         8-1         4     FALSE     FALSE
#>   [30]         8-2         4      TRUE     FALSE
#>   [31]         8-3         4      TRUE     FALSE
#>   [32]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_a3ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr4       1-5      + |         1         1         1 -0.394038
#>    [2]     chr4     11-15      + |         1         1         2 -0.394038
#>    [3]     chr4     21-25      + |         1         1         3 -0.394038
#>    [4]     chr4     31-35      + |         1         1         4 -0.394038
#>    [5]     chr4       1-5      + |         1         2         1  0.404634
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]     chr4   121-125      + |         2         7         4  0.258833
#>   [29]     chr4     91-95      + |         2         8         1  0.423911
#>   [30]     chr4   101-105      + |         2         8         2  0.423911
#>   [31]     chr4   111-115      + |         2         8         3  0.423911
#>   [32]     chr4   121-125      + |         2         8         4  0.423911
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE     FALSE
#>   [29]         8-1         4     FALSE     FALSE
#>   [30]         8-2         4      TRUE     FALSE
#>   [31]         8-3         4      TRUE     FALSE
#>   [32]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_atss(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr16       1-5      + |         1         1         1 -0.0908455
#>    [2]    chr16     11-15      + |         1         1         2 -0.0908455
#>    [3]    chr16     21-25      + |         1         1         3 -0.0908455
#>    [4]    chr16     31-35      + |         1         1         4 -0.0908455
#>    [5]    chr16       1-5      + |         1         2         1  0.7974239
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [28]    chr16   121-125      + |         2         7         4   0.316644
#>   [29]    chr16     91-95      + |         2         8         1  -0.481595
#>   [30]    chr16   101-105      + |         2         8         2  -0.481595
#>   [31]    chr16   111-115      + |         2         8         3  -0.481595
#>   [32]    chr16   121-125      + |         2         8         4  -0.481595
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE     FALSE
#>   [29]         8-1         4     FALSE     FALSE
#>   [30]         8-2         4      TRUE     FALSE
#>   [31]         8-3         4      TRUE     FALSE
#>   [32]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# the non-overlapping variant: drop the first exon altogether
generate_atss(gr, n_events = 1, mode = "drop")
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr16       1-5      + |         1         1         1 -0.0908455
#>    [2]    chr16     11-15      + |         1         1         2 -0.0908455
#>    [3]    chr16     21-25      + |         1         1         3 -0.0908455
#>    [4]    chr16     31-35      + |         1         1         4 -0.0908455
#>    [5]    chr16       1-5      + |         1         2         1  0.7974239
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [27]    chr16   121-125      + |         2         7         4   0.316644
#>   [28]    chr16     91-95      + |         2         8         1  -0.481595
#>   [29]    chr16   101-105      + |         2         8         2  -0.481595
#>   [30]    chr16   111-115      + |         2         8         3  -0.481595
#>   [31]    chr16   121-125      + |         2         8         4  -0.481595
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [27]         7-4         4     FALSE     FALSE
#>   [28]         8-1         4     FALSE     FALSE
#>   [29]         8-2         4      TRUE     FALSE
#>   [30]         8-3         4      TRUE     FALSE
#>   [31]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_ates(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr22       1-5      + |         1         1         1  -0.782176
#>    [2]    chr22     11-15      + |         1         1         2  -0.782176
#>    [3]    chr22     21-25      + |         1         1         3  -0.782176
#>    [4]    chr22     31-35      + |         1         1         4  -0.782176
#>    [5]    chr22       1-5      + |         1         2         1   0.337644
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [28]    chr22   121-125      + |         2         7         4 -0.2916265
#>   [29]    chr22     91-95      + |         2         8         1 -0.0883197
#>   [30]    chr22   101-105      + |         2         8         2 -0.0883197
#>   [31]    chr22   111-115      + |         2         8         3 -0.0883197
#>   [32]    chr22   121-125      + |         2         8         4 -0.0883197
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         4     FALSE     FALSE
#>    [2]         1-2         4      TRUE     FALSE
#>    [3]         1-3         4      TRUE     FALSE
#>    [4]         1-4         4     FALSE     FALSE
#>    [5]         2-1         4     FALSE     FALSE
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE     FALSE
#>   [29]         8-1         4     FALSE     FALSE
#>   [30]         8-2         4      TRUE     FALSE
#>   [31]         8-3         4      TRUE     FALSE
#>   [32]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
