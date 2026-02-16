# Takes a GRanges object and left/right exons to compute matches. i.e takes the pos_exons (GRanges) and left and right exons and returns a tibble that is the same as pos_exons but with two additional columns 'match_left' and 'match_right' indicating whether each exon matches the left and right exons respectively. The matching is done based on the 'type' parameter which can be "in", "over", or "boundary".

Takes a GRanges object and left/right exons to compute matches. i.e
takes the pos_exons (GRanges) and left and right exons and returns a
tibble that is the same as pos_exons but with two additional columns
'match_left' and 'match_right' indicating whether each exon matches the
left and right exons respectively. The matching is done based on the
'type' parameter which can be "in", "over", or "boundary".

## Usage

``` r
compute_matches(gr, left_exon, right_exon, type = c("in", "over", "boundary"))
```

## Arguments

- gr:

  A GRanges object with exon annotations

- left_exon:

  A GRanges object with left exon(s) to match

- right_exon:

  A GRanges object with right exon(s) to match

- type:

  The type of overlap to consider when identifying matches.

## Value

A tibble created from gr with two additional columns: 'match_left' and
'match_right' indicating the number of overlaps with the left and right
exons, respectively.
