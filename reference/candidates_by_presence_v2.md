# Filter candidates based on their presence in transcripts Then get the left and right exons for each candidate Return a named list with three tibbles: candidates, left_exons, right_exons

Filter candidates based on their presence in transcripts Then get the
left and right exons for each candidate Return a named list with three
tibbles: candidates, left_exons, right_exons

## Usage

``` r
candidates_by_presence_v2(gr, neg_exons, pos_exons)
```

## Arguments

- gr:

  A GRanges object with all exons

## Value

A named list with three granges: candidates, left_exons, right_exons
candidates, left_exons, right_exons are all from neg_exons set
