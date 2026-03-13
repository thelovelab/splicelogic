# Create a sample GRanges with one negative coef transcript and one positive coef transcripts This dataset is designed to include a mutually exclusive exons between exon_rank 3 of tx_id 1 and exon_rank 3 of tx_id 2. There is also a skipped exon event at exon_rank 5 of tx_id 1.

Create a sample GRanges with one negative coef transcript and one
positive coef transcripts This dataset is designed to include a mutually
exclusive exons between exon_rank 3 of tx_id 1 and exon_rank 3 of tx_id
2. There is also a skipped exon event at exon_rank 5 of tx_id 1.

## Usage

``` r
mx_mock_data()
```

## Value

A GRanges object with two transcripts per gene and candidate logic
