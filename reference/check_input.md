# Check that input is a valid GRanges object with required metadata columns exon_rank", "gene_id", "tx_id", coef_col

Check that input is a valid GRanges object with required metadata
columns exon_rank", "gene_id", "tx_id", coef_col

## Usage

``` r
check_input(gr, coef_col)
```

## Arguments

- gr:

  A GRanges object

- coef_col:

  Name of the coefficient metadata column (string)

## Value

TRUE if input is valid, otherwise throws an error
