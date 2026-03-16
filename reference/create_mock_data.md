# Create mock GRanges data for splicing event testing

Create mock GRanges data for splicing event testing

## Usage

``` r
create_mock_data(
  n_genes = 1,
  n_tx_per_gene = 2,
  n_exons_per_tx = 5,
  coef_range = c(-1, 1)
)
```

## Arguments

- n_genes:

  Number of genes to simulate

- n_tx_per_gene:

  Number of transcripts per gene

- n_exons_per_tx:

  Number of exons per transcript

- coef_range:

  Range of coefficient values to sample from

## Value

A GRanges object with simulated transcripts and exons

## Examples

``` r

# create mock data with 2 genes, 4 transcripts
# per gene, and 4 exons per transcript
gr <- create_mock_data(n_genes = 2, n_tx_per_gene = 4, n_exons_per_tx = 4)
```
