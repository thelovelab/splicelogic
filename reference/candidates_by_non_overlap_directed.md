# Filter candidates that do not overlap any pos_exons, and are internal Then get the left and right exons for each candidate Return a named list with three GRanges objects: candidates, left_exons, right_exons

Filter candidates that do not overlap any pos_exons, and are internal
Then get the left and right exons for each candidate Return a named list
with three GRanges objects: candidates, left_exons, right_exons

## Usage

``` r
candidates_by_non_overlap_directed(neg_exons, pos_exons, gr, type)
```

## Arguments

- neg_exons:

  A GRanges object with candidate negative exons (neg_exons)

- pos_exons:

  A GRanges object with positive exons (pos_exons)

- gr:

  The original GRanges object with all exons (for looking up left/right)

- type:

  The type of overlap to consider when identifying matches.

## Value

A named list with three GRanges objects: candidates, left_exons,
right_exons candidates, left_exons, right_exons are all from neg_exons
set
