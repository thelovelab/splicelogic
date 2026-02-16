# For a given candidate exon, find matching left and right exons in pos_exons

For a given candidate exon, find matching left and right exons in
pos_exons

## Usage

``` r
match_left_right(pos_exons, left_exon, right_exon, type)
```

## Arguments

- pos_exons:

  A GRanges object with positive exons (pos_exons)

- left_exon:

  A GRanges object with the left exon to match

- right_exon:

  A GRanges object with the right exon to match

- type:

  The type of overlap to consider when identifying matches.

## Value

A list with two tibbles: left_tbl and right_tbl
