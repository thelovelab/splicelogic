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
#> GRanges object with 31 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr4       1-5      + |         1         1         1 -0.620463
#>    [2]     chr4     11-15      + |         1         1         2 -0.620463
#>    [3]     chr4     21-25      + |         1         1         3 -0.620463
#>    [4]     chr4     31-35      + |         1         1         4 -0.620463
#>    [5]     chr4       1-5      + |         1         2         1  0.218208
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]     chr4   121-125      + |         2         7         4  0.812103
#>   [28]     chr4     91-95      + |         2         8         1  0.545461
#>   [29]     chr4   101-105      + |         2         8         2  0.545461
#>   [30]     chr4   111-115      + |         2         8         3  0.545461
#>   [31]     chr4   121-125      + |         2         8         4  0.545461
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [27]         7-4         4     FALSE
#>   [28]         8-1         4     FALSE
#>   [29]         8-2         4      TRUE
#>   [30]         8-3         4      TRUE
#>   [31]         8-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_mxe(gr, n_events = 1)
#> GRanges object with 30 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.0303205
#>    [2]    chr10     11-15      + |         1         1         2 -0.0303205
#>    [3]    chr10     21-25      + |         1         1         3 -0.0303205
#>    [4]    chr10     31-35      + |         1         1         4 -0.0303205
#>    [5]    chr10       1-5      + |         1         2         1  0.6024902
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [26]    chr10   121-125      + |         2         7         4   0.313992
#>   [27]    chr10     91-95      + |         2         8         1  -0.341366
#>   [28]    chr10   101-105      + |         2         8         2  -0.341366
#>   [29]    chr10   111-115      + |         2         8         3  -0.341366
#>   [30]    chr10   121-125      + |         2         8         4  -0.341366
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         1-1         4     FALSE
#>    [2]         1-2         4      TRUE
#>    [3]         1-3         4      TRUE
#>    [4]         1-4         4     FALSE
#>    [5]         2-1         4     FALSE
#>    ...         ...       ...       ...
#>   [26]         7-4         4     FALSE
#>   [27]         8-1         4     FALSE
#>   [28]         8-2         4      TRUE
#>   [29]         8-3         4      TRUE
#>   [30]         8-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_ri(gr, n_events = 1)
#> GRanges object with 31 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr16       1-5      + |         1         1         1 -0.607554
#>    [2]    chr16     11-15      + |         1         1         2 -0.607554
#>    [3]    chr16     21-25      + |         1         1         3 -0.607554
#>    [4]    chr16     31-35      + |         1         1         4 -0.607554
#>    [5]    chr16       1-5      + |         1         2         1  0.760473
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr16   111-115      + |         2         8         3  0.275371
#>   [28]    chr16   121-125      + |         2         8         4  0.275371
#>   [29]    chr16     91-95      + |         2         6         1  0.394466
#>   [30]    chr16   101-115      + |         2         6         2  0.394466
#>   [31]    chr16   121-125      + |         2         6         3  0.394466
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
#>   [29]         6-1         3     FALSE
#>   [30]         6-2         3      TRUE
#>   [31]         6-3         3     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_a5ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr1       1-5      + |         1         1         1 -0.568967
#>    [2]     chr1     11-15      + |         1         1         2 -0.568967
#>    [3]     chr1     21-25      + |         1         1         3 -0.568967
#>    [4]     chr1     31-35      + |         1         1         4 -0.568967
#>    [5]     chr1       1-5      + |         1         2         1  0.609315
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]     chr1   121-125      + |         2         7         4 0.0536606
#>   [29]     chr1     91-95      + |         2         8         1 0.5269497
#>   [30]     chr1   101-105      + |         2         8         2 0.5269497
#>   [31]     chr1   111-115      + |         2         8         3 0.5269497
#>   [32]     chr1   121-125      + |         2         8         4 0.5269497
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
  n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4
)
generate_a3ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank   estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer>  <numeric>
#>    [1]     chr5       1-5      + |         1         1         1 -0.1548135
#>    [2]     chr5     11-15      + |         1         1         2 -0.1548135
#>    [3]     chr5     21-25      + |         1         1         3 -0.1548135
#>    [4]     chr5     31-35      + |         1         1         4 -0.1548135
#>    [5]     chr5       1-5      + |         1         2         1  0.0406147
#>    ...      ...       ...    ... .       ...       ...       ...        ...
#>   [28]     chr5   121-125      + |         2         7         4  -0.423166
#>   [29]     chr5     91-95      + |         2         8         1  -0.970462
#>   [30]     chr5   101-105      + |         2         8         2  -0.970462
#>   [31]     chr5   111-115      + |         2         8         3  -0.970462
#>   [32]     chr5   121-125      + |         2         8         4  -0.970462
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
