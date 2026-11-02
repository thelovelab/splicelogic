# Plot the transcripts involved in splicelogic events, one plot per gene

Splits `events` by `gene_id` and draws each gene with
[`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md).

## Usage

``` r
plot_transcripts_by_gene(events, exons, genes = NULL, title = NULL, ...)
```

## Arguments

- events:

  A GRanges returned by
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  or any `find_*()` function, with a `gene_id` column.

- exons:

  A flat GRanges of exons with a `tx_id` column holding every transcript
  in `events`, e.g. the preprocessed GRanges the events were found in.

- genes:

  Optional character vector of gene ids to restrict to.

- title:

  Optional prefix for the plot titles, which are otherwise the gene id.

- ...:

  Passed to
  [`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md).

## Value

A named list of `ggplot` objects, one per gene.

## Examples

``` r
exons <- data.frame(
  seqnames = "chr1",
  start = c(1, 201, 401, 1, 401),
  end = c(100, 300, 500, 100, 500),
  strand = "+",
  gene_id = "g1",
  tx_id = c("down", "down", "down", "up", "up"),
  exon_rank = c(1, 2, 3, 1, 2),
  estimate = c(-1, -1, -1, 1, 1)
) |>
  plyranges::as_granges() |>
  preprocess(coef_col = "estimate")

events <- find_all_events(exons, verbose = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plots <- plot_transcripts_by_gene(events, exons)
  plots[["g1"]]
}

```
