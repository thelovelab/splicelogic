# Calculate all splicing events from a GRanges object

Runs all event detection functions and returns a single concatenated
GRanges with results from each.

## Usage

``` r
calc_all_events(gr, type = c("boundary", "over", "in"))
```

## Arguments

- gr:

  A GRanges object with exon annotations, preprocessed with
  preprocess_input().

- type:

  The type of overlap to consider for skipped exons, included exons, and
  mutually exclusive exons.

## Value

A GRanges object combining all detected events, with an 'event' metadata
column indicating the event type.
