# Find candidates and build flanking match tables for skipped/MX exon detection

Shared setup for calc_skipped_exons and calc_mx_exons: splits pos/neg
exons, finds candidates, runs batch overlap matching, and gene-filters
the results.

## Usage

``` r
find_candidates_and_flanks(gr, type, factor)
```

## Arguments

- gr:

  A preprocessed GRanges object

- type:

  Match type passed to find_matches_batch

- factor:

  Sign multiplier: 1 for normal, -1 for inverse (included exons)

## Value

A named list with pos_tbl, cand_tbl, left_match_tbl, right_match_tbl
tibbles for downstream processing. Returns NULL if no candidates are
found.
