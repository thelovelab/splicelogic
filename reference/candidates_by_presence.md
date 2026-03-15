# Filter candidates based on their presence in transcripts Then get the left and right exons for each candidate Return a named list with three GRanges objects: candidates, left_exons, right

Filter candidates based on their presence in transcripts Then get the
left and right exons for each candidate Return a named list with three
GRanges objects: candidates, left_exons, right

## Usage

``` r
candidates_by_presence(gr, coef_col)
```

## Arguments

- gr:

  A GRanges object with all exons

- coef_col:

  The name of the coefficient column in gr

## Value

A named list with three GRanges objects: candidates, left_exons,
right_exons. All are from the neg_exons set.
