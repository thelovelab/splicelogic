# Batch find overlaps between query and subject GRanges based on match type

Batch find overlaps between query and subject GRanges based on match
type

## Usage

``` r
find_matches_batch(query, subject, type = c("boundary", "over", "in"))
```

## Arguments

- query:

  A GRanges object (e.g. pos_exons)

- subject:

  A GRanges object (e.g. all left or right exons from the candidates
  (neg exons set))

- type:

  Match type: "over" for any overlap, "in" for exact match, "boundary"
  for shared start or end coordinate

## Value

A Hits object with queryHits and subjectHits indices
