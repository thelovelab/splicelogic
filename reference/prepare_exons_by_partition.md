# Prepare exons from two transcript partitions

Combines two transcript partitions (up- and down-regulated) and assigns
an `estimate` coefficient: `+1` to `up` and `-1` to `down`. Accepts
either GRanges objects or character vectors of transcript IDs (in which
case `txdb` is required to look up exon coordinates). The result is
ready to pass to
[`preprocess`](https://thelovelab.github.io/splicelogic/reference/preprocess.md)
with `coef_col = "estimate"`.

## Usage

``` r
prepare_exons_by_partition(
  up,
  down,
  txdb = NULL,
  tx_id_col = "TXNAME",
  verbose = TRUE
)
```

## Arguments

- up:

  A GRanges object or character vector of transcript IDs for the
  upregulated partition (assigned `estimate = +1`).

- down:

  A GRanges object or character vector of transcript IDs for the
  downregulated partition (assigned `estimate = -1`).

- txdb:

  A `TxDb` object (from GenomicFeatures). Required when `up` and `down`
  are character vectors of transcript IDs.

- tx_id_col:

  The keytype in `txdb` matching the transcript IDs in `up` and `down`.
  Default `"TXNAME"`. Only used when `up` and `down` are character
  vectors. See `AnnotationDbi::keytypes(txdb)` for available options.

- verbose:

  Whether to print progress messages. Default `TRUE`.

## Value

A combined GRanges object with an `estimate` column (`+1` for `up`, `-1`
for `down`), ready for
[`preprocess`](https://thelovelab.github.io/splicelogic/reference/preprocess.md).

## Details

When `up` and `down` are GRanges, both must have `exon_rank`, `gene_id`,
and `tx_id` metadata columns. Extra columns are kept; if one object
lacks a column present in the other, those entries receive `NA`.

## Examples

``` r

# GRanges input
gr <- create_mock_data(n_genes = 1, n_tx_per_gene = 4, n_exons_per_tx = 4)
gr <- generate_se(gr, n_events = 1)
gr_down <- gr[gr$estimate < 0]
gr_up <- gr[gr$estimate > 0]
prepare_exons_by_partition(gr_up, gr_down) |>
  preprocess(coef_col = "estimate")
#> GRanges object with 15 ranges and 7 metadata columns:
#>        seqnames    ranges strand |   gene_id     tx_id exon_rank  estimate
#>           <Rle> <IRanges>  <Rle> | <integer> <numeric> <integer> <integer>
#>    [1]    chr17       1-5      + |         1         2         1         1
#>    [2]    chr17     21-25      + |         1         2         2         1
#>    [3]    chr17     31-35      + |         1         2         3         1
#>    [4]    chr17       1-5      + |         1         1         1        -1
#>    [5]    chr17     11-15      + |         1         1         2        -1
#>    ...      ...       ...    ... .       ...       ...       ...       ...
#>   [11]    chr17     31-35      + |         1         3         4        -1
#>   [12]    chr17       1-5      + |         1         4         1        -1
#>   [13]    chr17     11-15      + |         1         4         2        -1
#>   [14]    chr17     21-25      + |         1         4         3        -1
#>   [15]    chr17     31-35      + |         1         4         4        -1
#>                key    nexons  internal
#>        <character> <integer> <logical>
#>    [1]         2-1         3     FALSE
#>    [2]         2-2         3      TRUE
#>    [3]         2-3         3     FALSE
#>    [4]         1-1         4     FALSE
#>    [5]         1-2         4      TRUE
#>    ...         ...       ...       ...
#>   [11]         3-4         4     FALSE
#>   [12]         4-1         4     FALSE
#>   [13]         4-2         4      TRUE
#>   [14]         4-3         4      TRUE
#>   [15]         4-4         4     FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
