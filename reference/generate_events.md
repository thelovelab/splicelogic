# Generate mock splice events in a GRanges object

Functions to introduce specific types of alternative splicing events
into mock GRanges data for testing purposes.

## Usage

``` r
generate_se(gr, n_events = 1)

generate_mxe(gr, n_events = 1)

generate_ri(gr, n_events = 1)

generate_a5ss(gr, n_events = 1)

generate_a3ss(gr, n_events = 1)
```

## Arguments

- gr:

  A GRanges object with metadata columns: 'exon_rank', 'gene_id',
  'tx_id', and 'estimate'.

- n_events:

  Number of events to generate

## Value

`generate_se()`: A GRanges object with skipped exon events introduced

`generate_mxe()`: A GRanges object with mutually exclusive exon events
introduced

`generate_ri()`: A GRanges object with retained intron events introduced

`generate_a5ss()`: A GRanges object with alternative 5' splice site
events introduced

`generate_a3ss()`: A GRanges object with alternative 3' splice site
events introduced

## Examples

``` r

gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_se(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr20       1-5      + |         1         1         1 -0.1333957
#>    [2]    chr20     11-15      + |         1         1         2 -0.1333957
#>    [3]    chr20     21-25      + |         1         1         3 -0.1333957
#>    [4]    chr20     31-35      + |         1         1         4 -0.1333957
#>    [5]    chr20       1-5      + |         1         2         1  0.0900365
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [27]    chr20   121-125      + |         2         7         4   0.562245
#>   [28]    chr20     91-95      + |         2         8         1  -0.423526
#>   [29]    chr20   101-105      + |         2         8         2  -0.423526
#>   [30]    chr20   111-115      + |         2         8         3  -0.423526
#>   [31]    chr20   121-125      + |         2         8         4  -0.423526
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
#>   [29]         8-2         4      TRUE      TRUE
#>   [30]         8-3         4      TRUE     FALSE
#>   [31]         8-4         4     FALSE     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_mxe(gr, n_events = 1)
#> GRanges object with 30 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr11       1-5      + |         1         1         1 -0.632443
#>    [2]    chr11     11-15      + |         1         1         2 -0.632443
#>    [3]    chr11     31-35      + |         1         1         3 -0.632443
#>    [4]    chr11       1-5      + |         1         2         1  0.729679
#>    [5]    chr11     21-25      + |         1         2         2  0.729679
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [26]    chr11   121-125      + |         2         7         4 -0.965969
#>   [27]    chr11     91-95      + |         2         8         1  0.346896
#>   [28]    chr11   101-105      + |         2         8         2  0.346896
#>   [29]    chr11   111-115      + |         2         8         3  0.346896
#>   [30]    chr11   121-125      + |         2         8         4  0.346896
#>                key    nexons  internal sim_event
#>        <character> <integer> <logical> <logical>
#>    [1]         1-1         3     FALSE     FALSE
#>    [2]         1-2         3      TRUE      TRUE
#>    [3]         1-3         3     FALSE     FALSE
#>    [4]         2-1         3     FALSE     FALSE
#>    [5]         2-2         3      TRUE      TRUE
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
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr11       1-5      + |         1         1         1 -0.331544
#>    [2]    chr11     11-15      + |         1         1         2 -0.331544
#>    [3]    chr11     21-25      + |         1         1         3 -0.331544
#>    [4]    chr11     31-35      + |         1         1         4 -0.331544
#>    [5]    chr11       1-5      + |         1         2         1  0.765840
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr11   111-115      + |         2         8         3 0.3814517
#>   [28]    chr11   121-125      + |         2         8         4 0.3814517
#>   [29]    chr11       1-5      + |         1         3         1 0.0651303
#>   [30]    chr11     11-25      + |         1         3         2 0.0651303
#>   [31]    chr11     31-35      + |         1         3         3 0.0651303
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
#>    [1]    chr10       1-5      + |         1         1         1 -0.785784
#>    [2]    chr10     11-15      + |         1         1         2 -0.785784
#>    [3]    chr10     21-25      + |         1         1         3 -0.785784
#>    [4]    chr10     31-35      + |         1         1         4 -0.785784
#>    [5]    chr10       1-5      + |         1         2         1  0.855186
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr10   121-125      + |         2         7         4 -0.747015
#>   [29]    chr10     91-95      + |         2         8         1  0.874660
#>   [30]    chr10   101-105      + |         2         8         2  0.874660
#>   [31]    chr10   111-117      + |         2         8         3  0.874660
#>   [32]    chr10   121-125      + |         2         8         4  0.874660
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
#>   [31]         8-3         4      TRUE      TRUE
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
#>    [1]    chr11       1-5      + |         1         1         1 -0.402287
#>    [2]    chr11     11-15      + |         1         1         2 -0.402287
#>    [3]    chr11     21-25      + |         1         1         3 -0.402287
#>    [4]    chr11     31-35      + |         1         1         4 -0.402287
#>    [5]    chr11       1-5      + |         1         2         1  0.334065
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr11   121-125      + |         2         7         4 -0.242718
#>   [29]    chr11     91-95      + |         2         8         1 -0.651227
#>   [30]    chr11   101-105      + |         2         8         2 -0.651227
#>   [31]    chr11   111-115      + |         2         8         3 -0.651227
#>   [32]    chr11   121-125      + |         2         8         4 -0.651227
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
