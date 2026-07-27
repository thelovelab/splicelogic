# Extract sequences for GRanges

Given a GRanges, look up the DNA or RNA sequence of each range in a
reference genome. This is a simple helper function that rearranges the
`getSeq` arguments and combines it with other helper functions to refer
to genomes by name.

## Usage

``` r
get_seq(gr, genome, as_rna = FALSE)
```

## Arguments

- gr:

  A GRanges object.

- genome:

  Either a character string naming an installed BSgenome package (e.g.
  `"hg38"`, or `"BSgenome.Hsapiens.UCSC.hg38"`) or a BSgenome object.
  See `?getBSgenome` for details.

- as_rna:

  If `TRUE`, return an `RNAStringSet` (`T -> U`). Default `FALSE`
  returns a `DNAStringSet`.

## Value

A `DNAStringSet` (or `RNAStringSet` if `as_rna = TRUE`), one entry per
range in `gr`, in the same order. Assign it onto `gr` (e.g.
`gr$seq <- ...`) to keep it as a metadata column.

## Examples

``` r

  gr <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(1e6, 1.1e6), width = 50)
  )

  suppressPackageStartupMessages(
    library(Biostrings)
  )

  get_seq(gr, "hg38")
#> DNAStringSet object of length 2:
#>     width seq
#> [1]    50 GGTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGA
#> [2]    50 CCAGCTATTCTGGAGACTGAGGCAGGAGGATCACTTGAGCCCAGGAGTTT
  get_seq(gr, "hg38", as_rna = TRUE)
#> RNAStringSet object of length 2:
#>     width seq
#> [1]    50 GGUGGAGCGCGCCGCCACGGACCACGGGCGGGCUGGCGGGCGAGCGGCGA
#> [2]    50 CCAGCUAUUCUGGAGACUGAGGCAGGAGGAUCACUUGAGCCCAGGAGUUU

  set.seed(123)
  gr <- create_mock_data(
    n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6
  ) |>
    generate_se(n_events = 1) |>
    GenomicRanges::shift(50e6) # move mock data
  
  skipped <- find_se(gr)

  suppressPackageStartupMessages(
    library(magrittr)
  )
  # magrittr pipe needed for the `.` placeholder below
  skipped %>% 
    plyranges::flank_upstream(100) %>%
    dplyr::mutate(seq = get_seq(., "hg38"))
#> GRanges object with 1 range and 8 metadata columns:
#>       seqnames            ranges strand |   gene_id     tx_id exon_rank
#>          <Rle>         <IRanges>  <Rle> | <integer> <numeric> <integer>
#>   [1]    chr15 49999911-50000010      + |         1         1         2
#>        estimate  event_type event_tx_id event_estimate                     seq
#>       <numeric> <character>   <numeric>      <numeric>          <DNAStringSet>
#>   [1] -0.454079          se           2       0.900827 TTTTCAAAAG...TCATCCAGAG
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
