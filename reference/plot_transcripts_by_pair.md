# Plot the transcripts involved in splicelogic events, one plot per pair

Splits `events` by transcript pair, so each plot shows two transcripts
and every event between them. The pair is unordered, so events reported
from either side land in the same plot.

## Usage

``` r
plot_transcripts_by_pair(events, exons, title = NULL, ...)
```

## Arguments

- events:

  A GRanges returned by
  [`find_all_events()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  or any `find_*()` function. A `tx_pair` column, if present, is used as
  the grouping as is.

- exons:

  A flat GRanges of exons with a `tx_id` column holding every transcript
  in `events`, e.g. the preprocessed GRanges the events were found in.

- title:

  Optional prefix for the plot titles, which are otherwise the
  transcript pair.

- ...:

  Passed to
  [`plot_transcripts()`](https://thelovelab.github.io/splicelogic/reference/plot_transcripts.md).

## Value

A named list of `ggplot` objects, one per transcript pair.

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
  plots <- plot_transcripts_by_pair(events, exons)
  plots[["down-up"]]
}

```
