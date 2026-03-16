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
  'tx_id', and 'coefs'.

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
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr21       1-5      + |         1         1         1 -0.121154
#>    [2]    chr21     11-15      + |         1         1         2 -0.121154
#>    [3]    chr21     21-25      + |         1         1         3 -0.121154
#>    [4]    chr21     31-35      + |         1         1         4 -0.121154
#>    [5]    chr21       1-5      + |         1         2         1  0.223939
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]    chr21   131-135      + |         2         7         4 -0.356291
#>   [28]    chr21   101-105      + |         2         8         1 -0.276932
#>   [29]    chr21   111-115      + |         2         8         2 -0.276932
#>   [30]    chr21   121-125      + |         2         8         3 -0.276932
#>   [31]    chr21   131-135      + |         2         8         4 -0.276932
#>                key    nexons  internal estimates
#>        <character> <integer> <logical> <numeric>
#>    [1]         1-1         4     FALSE -0.121154
#>    [2]         1-2         4      TRUE -0.121154
#>    [3]         1-3         4      TRUE -0.121154
#>    [4]         1-4         4     FALSE -0.121154
#>    [5]         2-1         4     FALSE  0.223939
#>    ...         ...       ...       ...       ...
#>   [27]         7-4         4     FALSE -0.356291
#>   [28]         8-1         4     FALSE -0.276932
#>   [29]         8-2         4      TRUE -0.276932
#>   [30]         8-3         4      TRUE -0.276932
#>   [31]         8-4         4     FALSE -0.276932
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_mx(gr, n_events = 1)
#> GRanges object with 30 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr20       1-5      + |         1         1         1 -0.244809
#>    [2]    chr20     11-15      + |         1         1         2 -0.244809
#>    [3]    chr20     21-25      + |         1         1         3 -0.244809
#>    [4]    chr20     31-35      + |         1         1         4 -0.244809
#>    [5]    chr20       1-5      + |         1         2         1  0.718034
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [26]    chr20   131-135      + |         2         7         3  0.518317
#>   [27]    chr20   101-105      + |         2         8         1  0.672151
#>   [28]    chr20   111-115      + |         2         8         2  0.672151
#>   [29]    chr20   121-125      + |         2         8         3  0.672151
#>   [30]    chr20   131-135      + |         2         8         4  0.672151
#>                key    nexons  internal estimates
#>        <character> <integer> <logical> <numeric>
#>    [1]         1-1         4     FALSE -0.244809
#>    [2]         1-2         4      TRUE -0.244809
#>    [3]         1-3         4      TRUE -0.244809
#>    [4]         1-4         4     FALSE -0.244809
#>    [5]         2-1         4     FALSE  0.718034
#>    ...         ...       ...       ...       ...
#>   [26]         7-3         3     FALSE  0.518317
#>   [27]         8-1         4     FALSE  0.672151
#>   [28]         8-2         4      TRUE  0.672151
#>   [29]         8-3         4      TRUE  0.672151
#>   [30]         8-4         4     FALSE  0.672151
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_retained_introns(gr, n_events = 1)
#> GRanges object with 31 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]     chr8       1-5      + |         1         1         1 -0.341499
#>    [2]     chr8     11-15      + |         1         1         2 -0.341499
#>    [3]     chr8     21-25      + |         1         1         3 -0.341499
#>    [4]     chr8     31-35      + |         1         1         4 -0.341499
#>    [5]     chr8       1-5      + |         1         3         1  0.958587
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [27]     chr8   121-125      + |         2         8         3  0.355956
#>   [28]     chr8   131-135      + |         2         8         4  0.355956
#>   [29]     chr8       1-5      + |         1         2         1  0.959888
#>   [30]     chr8     11-25      + |         1         2         2  0.959888
#>   [31]     chr8     31-35      + |         1         2         3  0.959888
#>                key    nexons  internal estimates
#>        <character> <integer> <logical> <numeric>
#>    [1]         1-1         4     FALSE -0.341499
#>    [2]         1-2         4      TRUE -0.341499
#>    [3]         1-3         4      TRUE -0.341499
#>    [4]         1-4         4     FALSE -0.341499
#>    [5]         3-1         4     FALSE  0.958587
#>    ...         ...       ...       ...       ...
#>   [27]         8-3         4      TRUE  0.355956
#>   [28]         8-4         4     FALSE  0.355956
#>   [29]         2-1         3     FALSE  0.959888
#>   [30]         2-2         3      TRUE  0.959888
#>   [31]         2-3         3     FALSE  0.959888
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_a5ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr16       1-5      + |         1         1         1 -0.805743
#>    [2]    chr16     11-15      + |         1         1         2 -0.805743
#>    [3]    chr16     21-25      + |         1         1         3 -0.805743
#>    [4]    chr16     31-35      + |         1         1         4 -0.805743
#>    [5]    chr16       1-5      + |         1         2         1  0.556948
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr16   131-135      + |         2         7         4  0.350417
#>   [29]    chr16   101-105      + |         2         8         1  0.892590
#>   [30]    chr16   111-115      + |         2         8         2  0.892590
#>   [31]    chr16   121-125      + |         2         8         3  0.892590
#>   [32]    chr16   131-135      + |         2         8         4  0.892590
#>                key    nexons  internal estimates
#>        <character> <integer> <logical> <numeric>
#>    [1]         1-1         4     FALSE -0.805743
#>    [2]         1-2         4      TRUE -0.805743
#>    [3]         1-3         4      TRUE -0.805743
#>    [4]         1-4         4     FALSE -0.805743
#>    [5]         2-1         4     FALSE  0.556948
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE  0.350417
#>   [29]         8-1         4     FALSE  0.892590
#>   [30]         8-2         4      TRUE  0.892590
#>   [31]         8-3         4      TRUE  0.892590
#>   [32]         8-4         4     FALSE  0.892590
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths


gr <- create_mock_data(
  n_genes = 2, n_tx = 4, n_exons = 4
)
generate_a3ss(gr, n_events = 1)
#> GRanges object with 32 ranges and 8 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank     coefs
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>    [1]    chr10       1-5      + |         1         1         1 -0.785784
#>    [2]    chr10     11-15      + |         1         1         2 -0.785784
#>    [3]    chr10     21-25      + |         1         1         3 -0.785784
#>    [4]    chr10     31-35      + |         1         1         4 -0.785784
#>    [5]    chr10       1-5      + |         1         2         1  0.855186
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [28]    chr10   131-135      + |         2         7         4 -0.747015
#>   [29]    chr10   101-105      + |         2         8         1  0.874660
#>   [30]    chr10   111-115      + |         2         8         2  0.874660
#>   [31]    chr10   123-125      + |         2         8         3  0.874660
#>   [32]    chr10   131-135      + |         2         8         4  0.874660
#>                key    nexons  internal estimates
#>        <character> <integer> <logical> <numeric>
#>    [1]         1-1         4     FALSE -0.785784
#>    [2]         1-2         4      TRUE -0.785784
#>    [3]         1-3         4      TRUE -0.785784
#>    [4]         1-4         4     FALSE -0.785784
#>    [5]         2-1         4     FALSE  0.855186
#>    ...         ...       ...       ...       ...
#>   [28]         7-4         4     FALSE -0.747015
#>   [29]         8-1         4     FALSE  0.874660
#>   [30]         8-2         4      TRUE  0.874660
#>   [31]         8-3         4      TRUE  0.874660
#>   [32]         8-4         4     FALSE  0.874660
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
