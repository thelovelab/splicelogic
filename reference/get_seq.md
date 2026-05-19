# Extract sequences for a GRanges

Given a GRanges (e.g. the output of
[`find_se`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
[`find_a5ss`](https://thelovelab.github.io/splicelogic/reference/find_events.md),
etc.), look up the DNA or RNA sequence of each range from a BSgenome
reference.

## Usage

``` r
get_seq(gr, genome, as_rna = FALSE)
```

## Arguments

- gr:

  A GRanges object.

- genome:

  Either a character string naming an installed BSgenome package (e.g.
  `"BSgenome.Hsapiens.UCSC.hg38"`) or a BSgenome object.

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
  get_seq(gr, "BSgenome.Hsapiens.UCSC.hg38")
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://p3m.dev/cran/__linux__/noble/latest
#> Error in .stopOnAvailablePkg(genome): BSgenome.Hsapiens.UCSC.hg38 package is not currently installed.
#>   You first need to install it, which you can do with:
#>       library(BiocManager)
#>       install("BSgenome.Hsapiens.UCSC.hg38")
  get_seq(gr, "BSgenome.Hsapiens.UCSC.hg38", as_rna = TRUE)
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://p3m.dev/cran/__linux__/noble/latest
#> Error in .stopOnAvailablePkg(genome): BSgenome.Hsapiens.UCSC.hg38 package is not currently installed.
#>   You first need to install it, which you can do with:
#>       library(BiocManager)
#>       install("BSgenome.Hsapiens.UCSC.hg38")

# On splicelogic output: upstream-flank RNA of a skipped exon.
  gr <- create_mock_data(
    n_genes = 2, n_tx_per_gene = 3, n_exons_per_tx = 6
  ) |>
    generate_se(n_events = 1)
  skipped <- find_se(gr)
  # shift mock coords into a real hg38 region before lookup
  skipped <- GenomicRanges::shift(skipped, 1e6)
  flanked <- plyranges::flank_upstream(skipped, 100)
  flanked$seq <- get_seq(flanked, genome = "hg38", as_rna = TRUE)
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://p3m.dev/cran/__linux__/noble/latest
#> Error in .stopOnAvailablePkg(genome): BSgenome.Hsapiens.UCSC.hg38 package is not currently installed.
#>   You first need to install it, which you can do with:
#>       library(BiocManager)
#>       install("BSgenome.Hsapiens.UCSC.hg38")
  flanked
#> GRanges object with 1 range and 7 metadata columns:
#>       seqnames         ranges strand |   gene_id     tx_id exon_rank  estimate
#>          <Rle>      <IRanges>  <Rle> | <integer> <numeric> <integer> <numeric>
#>   [1]    chr15 999911-1000010      + |         1         1         2 -0.984157
#>        event_type event_tx_id event_estimate
#>       <character>   <numeric>      <numeric>
#>   [1]          se           2       0.700679
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
