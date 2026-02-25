# Introduction to splicelogic

## Introduction

`splicelogic` is an R/Bioconductor package for detecting alternative
splicing events from exon-level data stored as `GRanges` objects. Given
a set of exons annotated with a coefficient column indicating
differential transcript usage (DTU), `splicelogic` identifies the
following types of splicing events:

- **Skipped exons (SE)** – exons present in one isoform but absent in
  another.
- **Mutually exclusive exons (MXE)** – pairs of exons where one is
  included at the expense of the other.
- **Retained introns (RI)** – intronic regions that are retained as part
  of an exon in an alternative isoform.
- **Alternative 5’ (A5SS)** – exons that share the same 3’ splice site
  but differ at the 5’ end.
  - **Alternative 3’ (A3SS)** – exons that share the same 5’ splice site
    but differ at the 3’ end.

## Quick start

With DTU results attached to a *GRanges* of the exons from significant
transcripts, one can use the following code to identify splice events:

``` r

exons <- preprocess_input(exons, coef_col = "coefs")
skipped <- exons |> calc_skipped_exons()
mut_exc <- exons |> calc_mutually_exclusive()
# etc.
```

## Input data

`splicelogic` assumes the user has run some kind of DTU (e.g. isoform
switching) statistical analysis providing error bound (FDR) and effect
estimate (e.g. model coefficient or deltaPSI) (see [upstream
methods](#upstream-methods)). It also assumes the user can obtain ranges
representing the exon structure of each transcript being analyzed (see
[obtaining exon ranges](#obtaining-exon-ranges). The exons should be in
a flat *GRanges* object (one range per exon), containing exon-level
metadata such as the gene ID, transcript ID, and rank in the transcript.
Information from the DTU analysis should also be present in the
metadata. Code below demonstrates generating this from a GTF file and
attaching metdata from a DTU results table.

| Column      | Description                                       |
|-------------|---------------------------------------------------|
| `gene_id`   | Gene identifier                                   |
| `tx_id`     | Transcript identifier                             |
| `exon_rank` | Exon rank within the transcript                   |
| *coef_col*  | DTU effect estimate (column name is user-defined) |

The coefficient input in *coef_col* indicates the differential
transcript usage (DTU) of the specific transcript containing the exon.
All exons from the same transcript share the same coefficient value.
This values comes from a prior differential transcript usage analysis,
such as that performed by
[satuRn](https://bioconductor.org/packages/satuRn). Positive coefficient
values indicate upregulated exons and negative values indicate
downregulated exons.

### Obtaining exon ranges

Here we demontrate reading in exon ranges from a GTF and adding example
DTU results.

``` r

# first we get an example TxDB for our exon ranges
suppressPackageStartupMessages({
  library(AnnotationHub)
  library(GenomicFeatures)
  library(tibble)
})
# here we will use the transcript database for GENCODE v32 from AHub
# typically, you should supply your own GTF to makeTxDbFromGFF()
ah <- AnnotationHub()
txdb <- ah[["AH75191"]]
```

    ## loading from cache

The following extracts the exons from the *TxDb*:

``` r

# extract exons as a GRangesList
ebt <- GenomicFeatures::exonsBy(txdb, by="tx") # exon id, name, and rank
```

Here we build a transcript table, likely this is already provided by an
upstream DTU method.

``` r

# extract transcript table as tibble
txps <- txdb |>
  AnnotationDbi::select(keys(txdb, "TXID"), c("TXNAME","GENEID"), "TXID") |>
  tibble::as_tibble() |>
  dplyr::select(tx_num = TXID, tx_id = TXNAME, gene_id = GENEID)
```

    ## 'select()' returned 1:1 mapping between keys and columns

Suppose we already have DTU results for transcripts:

``` r

# here just simulated results
txps <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )
# show the top 3
txps[1:3,]
```

    ## # A tibble: 3 × 5
    ##   tx_num tx_id             gene_id             padj effect_est
    ##    <int> <chr>             <chr>              <dbl>      <dbl>
    ## 1      1 ENST00000456328.2 ENSG00000223972.5 0.0808     -0.497
    ## 2      2 ENST00000450305.2 ENSG00000223972.5 0.834       0.938
    ## 3      3 ENST00000473358.1 ENSG00000243485.5 0.601      -1.84

``` r

# check the concordance...
all.equal(txps$tx_num, as.integer(names(ebt)))
```

    ## [1] TRUE

``` r

# and if aligned, set new names for exons
exons <- ebt
names(exons) <- txps$tx_id
```

Next flattening the exons:

``` r

exons <- unlist(exons)
exons$tx_id <- names(exons)
names(exons) <- exons$exon_name
exons[1:3] # show the top 3
```

    ## GRanges object with 3 ranges and 4 metadata columns:
    ##                     seqnames      ranges strand |   exon_id         exon_name
    ##                        <Rle>   <IRanges>  <Rle> | <integer>       <character>
    ##   ENSE00002234944.1     chr1 11869-12227      + |         1 ENSE00002234944.1
    ##   ENSE00003582793.1     chr1 12613-12721      + |         5 ENSE00003582793.1
    ##   ENSE00002312635.1     chr1 13221-14409      + |         8 ENSE00002312635.1
    ##                     exon_rank             tx_id
    ##                     <integer>       <character>
    ##   ENSE00002234944.1         1 ENST00000456328.2
    ##   ENSE00003582793.1         2 ENST00000456328.2
    ##   ENSE00002312635.1         3 ENST00000456328.2
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

Adding DTU results and gene ID:

``` r

# arrange `txps` as in `exons`, including duplicates
txp_idx <- match(exons$tx_id, txps$tx_id)
add_columns <- txps[txp_idx,] |> dplyr::select(-c(tx_id, tx_num))
merged_DF <- cbind(mcols(exons), add_columns)
mcols(exons) <- merged_DF
```

We end up with all the exons-per-transcript in a flat *GRanges* object
with the following columns:

``` r

exons[1:3] |> mcols()
```

    ## DataFrame with 3 rows and 7 columns
    ##                     exon_id         exon_name exon_rank             tx_id
    ##                   <integer>       <character> <integer>       <character>
    ## ENSE00002234944.1         1 ENSE00002234944.1         1 ENST00000456328.2
    ## ENSE00003582793.1         5 ENSE00003582793.1         2 ENST00000456328.2
    ## ENSE00002312635.1         8 ENSE00002312635.1         3 ENST00000456328.2
    ##                             gene_id      padj effect_est
    ##                         <character> <numeric>  <numeric>
    ## ENSE00002234944.1 ENSG00000223972.5 0.0807501  -0.496625
    ## ENSE00003582793.1 ENSG00000223972.5 0.0807501  -0.496625
    ## ENSE00002312635.1 ENSG00000223972.5 0.0807501  -0.496625

### Upstream methods

Something about what methods can be used for upstream DTU or switching
analysis.

## How to start detecting splicing events

``` r

library(splicelogic)
```

### Preprocessing

We can create mock data for demonstration purposes using the
[`create_mock_data()`](https://thelovelab.github.io/splicelogic/reference/create_mock_data.md)
function. This function generates a `GRanges` object with the required
metadata columns and random coefficient values for each transcript,
making sure that at least one transcript is upregulated (coef\>0) and
one is downregulated (coef\<0) for each gene.

``` r

gr <- create_mock_data( n_genes = 2, n_tx = 4, n_exons = 6 )
mcols(gr)
```

    ## DataFrame with 48 rows and 7 columns
    ##       gene_id     tx_id exon_rank     coefs         key    nexons  internal
    ##     <integer> <numeric> <integer> <numeric> <character> <integer> <logical>
    ## 1           1         1         1 -0.821507         1-1         6     FALSE
    ## 2           1         1         2 -0.821507         1-2         6      TRUE
    ## 3           1         1         3 -0.821507         1-3         6      TRUE
    ## 4           1         1         4 -0.821507         1-4         6      TRUE
    ## 5           1         1         5 -0.821507         1-5         6      TRUE
    ## ...       ...       ...       ...       ...         ...       ...       ...
    ## 44          2         8         2 -0.707009         8-2         6      TRUE
    ## 45          2         8         3 -0.707009         8-3         6      TRUE
    ## 46          2         8         4 -0.707009         8-4         6      TRUE
    ## 47          2         8         5 -0.707009         8-5         6      TRUE
    ## 48          2         8         6 -0.707009         8-6         6     FALSE

### Skipped exons

[`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_skipped_exons.md)
identifies exons that are present in downregulated transcripts but
absent (skipped) in upregulated ones.

Using the gr object generated by
[`create_mock_data()`](https://thelovelab.github.io/splicelogic/reference/create_mock_data.md),
we can generate skiped exon events with the
[`generate_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/generate_skipped_exons.md)
function, which modifies the input `GRanges` object to include skipped
exon events. Then we can run the detection function
[`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_skipped_exons.md)
on the modified `GRanges` object:

``` r

gr_se <- generate_skipped_exons(gr, n_se = 2)
gr_se <- preprocess_input(gr_se, coef_col = "coefs")
se_result <- calc_skipped_exons(gr_se, coef_col = "coefs")
se_result
```

The result is a `GRanges` object containing only the exons flagged as
skipped, with an `event` column set to `"skipped_exon"`. The `"tx_id"`
refers to the transcript that includes the skipped exon, `"tx_event"`
indicates the transcripts in respect to which the exon is skipped
(i.e. the downregulated transcripts that do not include the exon).

``` r

mcols(se_result)[, c("gene_id", "tx_id", "n_txp", "event", "tx_event")]
```

### Mutually exclusive exons

[`calc_mutually_exclusive()`](https://thelovelab.github.io/splicelogic/reference/calc_mutually_exclusive.md)
detects pairs of exons where one exon is included in one isoform and a
different exon occupies the equivalent position in another isoform.

The
[`mx_mock_data()`](https://thelovelab.github.io/splicelogic/reference/mx_mock_data.md)
dataset is designed to contain mutually exclusive exon events:

``` r

gr_mx <- mx_mock_data()
gr_mx <- preprocess_input(gr_mx, coef_col = "coefs")
mx_result <- calc_mutually_exclusive(gr_mx, coef_col = "coefs")
mx_result
```

``` r

mcols(mx_result)[, c("gene_id", "tx_id", "exon_rank", "event")]
```

### Retained introns

[`calc_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/calc_retained_introns.md)
finds exons in upregulated transcripts that overlap intronic regions of
downregulated transcripts, indicating intron retention.

``` r

gr_ri <- create_mock_data(3, 6, 3)
gr_ri <- preprocess_input(gr_ri, coef_col = "coefs")
gr_ri <- generate_retained_introns(gr_ri, n_ri = 3)
ri_result <- calc_retained_introns(gr_ri, coef_col = "coefs")
ri_result
```

``` r

if (length(ri_result) > 0) {
    mcols(ri_result)[, c("gene_id", "tx_id", "exon_rank", "event")]
}
```

### Alternative 3’ and 5’ splice sites

**NOTE: I think it’s more common to list 5 then 3**

[`calc_a3ss_a5ss()`](https://thelovelab.github.io/splicelogic/reference/calc_a3ss_a5ss.md)
detects exons that share one boundary (start or end) with an exon in
another isoform but differ at the opposite end, indicating alternative
splice site usage.

``` r

gr_ss <- create_mock_data(3, 3, 6)
gr_ss <- preprocess_input(gr_ss, coef_col = "coefs")
gr_ss <- generate_a3ss(gr_ss, n_a3ss = 3)
ss_result <- calc_a3ss_a5ss(gr_ss, coef_col = "coefs")
ss_result
```

``` r

if (length(ss_result) > 0) {
    mcols(ss_result)[, c("gene_id", "tx_id", "exon_rank", "event")]
}
```

## Handling cases with no events

When the input data does not contain any splicing events, each detection
function returns an empty `GRanges` object:

``` r

gr_none <- no_event_mock_data()
gr_none <- preprocess_input(gr_none, coef_col = "coefs")

calc_skipped_exons(gr_none, coef_col = "coefs")
calc_mutually_exclusive(gr_none, coef_col = "coefs")
```

## Session info

``` r

sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] splicelogic_0.0.67     tibble_3.3.1           GenomicFeatures_1.62.0
    ##  [4] AnnotationDbi_1.72.0   Biobase_2.70.0         GenomicRanges_1.62.1  
    ##  [7] Seqinfo_1.0.0          IRanges_2.44.0         S4Vectors_0.48.0      
    ## [10] AnnotationHub_4.0.0    BiocFileCache_3.0.0    dbplyr_2.5.2          
    ## [13] BiocGenerics_0.56.0    generics_0.1.4        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1            dplyr_1.2.0                
    ##  [3] blob_1.3.0                  filelock_1.0.3             
    ##  [5] Biostrings_2.78.0           bitops_1.0-9               
    ##  [7] fastmap_1.2.0               RCurl_1.98-1.17            
    ##  [9] GenomicAlignments_1.46.0    XML_3.99-0.22              
    ## [11] digest_0.6.39               lifecycle_1.0.5            
    ## [13] plyranges_1.30.1            KEGGREST_1.50.0            
    ## [15] RSQLite_2.4.6               magrittr_2.0.4             
    ## [17] compiler_4.5.2              rlang_1.1.7                
    ## [19] sass_0.4.10                 tools_4.5.2                
    ## [21] utf8_1.2.6                  yaml_2.3.12                
    ## [23] rtracklayer_1.70.1          knitr_1.51                 
    ## [25] S4Arrays_1.10.1             htmlwidgets_1.6.4          
    ## [27] bit_4.6.0                   curl_7.0.0                 
    ## [29] DelayedArray_0.36.0         abind_1.4-8                
    ## [31] BiocParallel_1.44.0         withr_3.0.2                
    ## [33] purrr_1.2.1                 desc_1.4.3                 
    ## [35] grid_4.5.2                  SummarizedExperiment_1.40.0
    ## [37] cli_3.6.5                   rmarkdown_2.30             
    ## [39] crayon_1.5.3                ragg_1.5.0                 
    ## [41] otel_0.2.0                  httr_1.4.8                 
    ## [43] rjson_0.2.23                DBI_1.2.3                  
    ## [45] cachem_1.1.0                parallel_4.5.2             
    ## [47] BiocManager_1.30.27         XVector_0.50.0             
    ## [49] restfulr_0.0.16             matrixStats_1.5.0          
    ## [51] vctrs_0.7.1                 Matrix_1.7-4               
    ## [53] jsonlite_2.0.0              bit64_4.6.0-1              
    ## [55] systemfonts_1.3.1           jquerylib_0.1.4            
    ## [57] glue_1.8.0                  pkgdown_2.2.0              
    ## [59] codetools_0.2-20            BiocVersion_3.22.0         
    ## [61] BiocIO_1.20.0               pillar_1.11.1              
    ## [63] rappdirs_0.3.4              htmltools_0.5.9            
    ## [65] R6_2.6.1                    httr2_1.2.2                
    ## [67] textshaping_1.0.4           evaluate_1.0.5             
    ## [69] lattice_0.22-9              png_0.1-8                  
    ## [71] Rsamtools_2.26.0            cigarillo_1.0.0            
    ## [73] memoise_2.0.1               bslib_0.10.0               
    ## [75] SparseArray_1.10.8          xfun_0.56                  
    ## [77] fs_1.6.6                    MatrixGenerics_1.22.0      
    ## [79] pkgconfig_2.0.3
