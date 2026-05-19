# Prepare exon ranges from a TxDb and DTU results table

Extracts exon ranges from a TxDb object, merges them with differential
transcript usage (DTU) results, and returns a flat GRanges ready for
[`preprocess`](https://thelovelab.github.io/splicelogic/reference/preprocess.md).

## Usage

``` r
prepare_exons(
  txdb,
  dtu_table,
  coef_col,
  tx_id_col = "tx_id",
  gene_id_col = "gene_id",
  verbose = TRUE
)
```

## Arguments

- txdb:

  A `TxDb` object (from GenomicFeatures).

- dtu_table:

  A data.frame or tibble with DTU results. Must contain columns for
  transcript ID, gene ID, and a coefficient.

- coef_col:

  Column name in `dtu_table` with the coefficient / effect size values.

- tx_id_col:

  Column name in `dtu_table` with transcript IDs matching the TxDb
  transcript names. Default `"tx_id"`.

- gene_id_col:

  Column name in `dtu_table` with gene IDs. Default `"gene_id"`.

- verbose:

  Whether to print progress messages. Default `TRUE`.

## Value

A GRanges object with metadata columns: `gene_id`, `tx_id`, `exon_rank`,
the coefficient column, and any additional columns from `dtu_table`.

## Examples

``` r

library(AnnotationHub)
#> Loading required package: BiocFileCache
#> Loading required package: dbplyr
#> 
#> Attaching package: ‘AnnotationHub’
#> The following object is masked from ‘package:rtracklayer’:
#> 
#>     hubUrl
library(AnnotationDbi)
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:AnnotationHub’:
#> 
#>     cache
library(GenomicFeatures)
library(tibble)

ah <- AnnotationHub()
txdb <- ah[["AH84134"]] # fly TxDb (Drosophila melanogaster)
#> loading from cache

# build a simulated DTU table from the TxDb transcripts
txps <- txdb |>
  AnnotationDbi::select(
    keys(txdb, "TXID"), c("TXNAME", "GENEID"), "TXID"
  ) |>
  tibble::as_tibble() |>
  dplyr::select(tx_id = TXNAME, gene_id = GENEID)|>
  dplyr::filter(!is.na(gene_id))
#> 'select()' returned 1:1 mapping between keys and columns

sim_dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )

fly_exons <- prepare_exons(
  txdb, sim_dtu_table, coef_col = "effect_est", verbose = TRUE
)
#> Extracting exons from TxDb...
#> 'select()' returned 1:1 mapping between keys and columns
#> Mapping transcript IDs...
#> Merging DTU results onto exons...
#> Done. Returned 188169 exon ranges from 34920 unique transcripts.
```
