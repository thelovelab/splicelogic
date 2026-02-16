# Combine two GRanges objects into one for splicing event analysis This function takes two GRanges objects, typically representing positive and negative sets of exons, and combines them into a single GRanges object.

Combine two GRanges objects into one for splicing event analysis This
function takes two GRanges objects, typically representing positive and
negative sets of exons, and combines them into a single GRanges object.

## Usage

``` r
combine_gr_input(gr1, gr2, coef_col)
```

## Arguments

- gr1:

  A GRanges object (e.g., positive set)

- gr2:

  A GRanges object (e.g., negative set)

- coef_col:

  Name of the coefficient metadata column (string)

## Value

A combined GRanges object with appropriate coef metadata
