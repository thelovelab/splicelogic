# Generate mock splice events in a GRanges object

Functions to introduce specific types of alternative splicing events
into mock GRanges data for testing purposes.

## Usage

``` r
generate_skipped_exons(gr, n_events = 1)

generate_mx(gr, n_events = 1)

generate_retained_introns(gr, n_events = 1)

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

`generate_skipped_exons()`: A GRanges object with skipped exon events
introduced

`generate_mx()`: A GRanges object with mutually exclusive exon events
introduced

`generate_retained_introns()`: A GRanges object with retained intron
events introduced

`generate_a5ss()`: A GRanges object with alternative 5' splice site
events introduced

`generate_a3ss()`: A GRanges object with alternative 3' splice site
events introduced

## Examples

``` r

gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_skipped_exons(gr, n_events = 1)
#> GRanges object with 31 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.103009
#>    [2]    chr10     11-15      + |         1         1         2 -0.103009
#>    [3]    chr10     21-25      + |         1         1         3 -0.103009
#>    [4]    chr10     31-35      + |         1         1         4 -0.103009
#>    [5]    chr10       1-5      + |         1         2         1  0.855840
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr10   131-135      + |         2         7         4  0.656029
#>   [28]    chr10   101-105      + |         2         8         1 -0.798687
#>   [29]    chr10   111-115      + |         2         8         2 -0.798687
#>   [30]    chr10   121-125      + |         2         8         3 -0.798687
#>   [31]    chr10   131-135      + |         2         8         4 -0.798687
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         3     FALSE
#>    ...         ...       ...       ...
#>   [27]         7-4         4     FALSE
#>   [28]         8-1         4     FALSE
#>   [29]         8-2         4      TRUE
#>   [30]         8-3         4      TRUE
#>   [31]         8-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_mx(gr, n_events = 1)
#> GRanges object with 30 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr9       1-5      + |         1         1         1 -0.404016
#>    [2]     chr9     11-15      + |         1         1         2 -0.404016
#>    [3]     chr9     21-25      + |         1         1         3 -0.404016
#>    [4]     chr9     31-35      + |         1         1         4 -0.404016
#>    [5]     chr9       1-5      + |         1         2         1  0.395668
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [26]     chr9   111-115      + |         2         7         2 -0.838310
#>   [27]     chr9   131-135      + |         2         7         3 -0.838310
#>   [28]     chr9   101-105      + |         2         8         1  0.311965
#>   [29]     chr9   121-125      + |         2         8         2  0.311965
#>   [30]     chr9   131-135      + |         2         8         3  0.311965
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [26]         7-2         3      TRUE
#>   [27]         7-3         3     FALSE
#>   [28]         8-1         3     FALSE
#>   [29]         8-2         3      TRUE
#>   [30]         8-3         3     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_retained_introns(gr, n_events = 1)
#> GRanges object with 31 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr8       1-5      + |         1         1         1 -0.586908
#>    [2]     chr8     11-15      + |         1         1         2 -0.586908
#>    [3]     chr8     21-25      + |         1         1         3 -0.586908
#>    [4]     chr8     31-35      + |         1         1         4 -0.586908
#>    [5]     chr8       1-5      + |         1         2         1  0.134598
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]     chr8   121-125      + |         2         8         3 0.0177871
#>   [28]     chr8   131-135      + |         2         8         4 0.0177871
#>   [29]     chr8       1-5      + |         1         4         1 0.3303049
#>   [30]     chr8     11-25      + |         1         4         2 0.3303049
#>   [31]     chr8     31-35      + |         1         4         3 0.3303049
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [27]         8-3         4      TRUE
#>   [28]         8-4         4     FALSE
#>   [29]         4-1         3     FALSE
#>   [30]         4-2         3      TRUE
#>   [31]         4-3         3     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_a5ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr20       1-5      + |         1         1         1 -0.478438
#>    [2]    chr20     11-15      + |         1         1         2 -0.478438
#>    [3]    chr20     21-25      + |         1         1         3 -0.478438
#>    [4]    chr20     31-35      + |         1         1         4 -0.478438
#>    [5]    chr20       1-5      + |         1         2         1  0.886854
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr20   131-135      + |         2         7         4  0.629192
#>   [29]    chr20   101-105      + |         2         8         1 -0.858070
#>   [30]    chr20   111-115      + |         2         8         2 -0.858070
#>   [31]    chr20   121-125      + |         2         8         3 -0.858070
#>   [32]    chr20   131-135      + |         2         8         4 -0.858070
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [28]         7-4         4     FALSE
#>   [29]         8-1         4     FALSE
#>   [30]         8-2         4      TRUE
#>   [31]         8-3         4      TRUE
#>   [32]         8-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_a3ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr14       1-5      + |         1         1         1 -0.741426
#>    [2]    chr14     11-15      + |         1         1         2 -0.741426
#>    [3]    chr14     21-25      + |         1         1         3 -0.741426
#>    [4]    chr14     31-35      + |         1         1         4 -0.741426
#>    [5]    chr14       1-5      + |         1         2         1  0.367170
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr14   131-135      + |         2         7         4  0.655003
#>   [29]    chr14   101-105      + |         2         8         1 -0.744885
#>   [30]    chr14   111-115      + |         2         8         2 -0.744885
#>   [31]    chr14   121-125      + |         2         8         3 -0.744885
#>   [32]    chr14   131-135      + |         2         8         4 -0.744885
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [28]         7-4         4     FALSE
#>   [29]         8-1         4     FALSE
#>   [30]         8-2         4      TRUE
#>   [31]         8-3         4      TRUE
#>   [32]         8-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
