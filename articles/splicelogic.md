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
mut_exc <- exons |> calc_mx_exons()
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

### Jones et al mouse long read dataset

We will use a published long read dataset and it’s reported splicing
changes to demonstrate some of the functionality in *splicelogic*.

The citation is:

> Emma F. Jones, Timothy C. Howton, Victoria L. Flanary, Amanda D.
> Clark, Brittany N. Lasseigne Long-read RNA sequencing identifies
> region- and sex-specific C57BL/6J mouse brain mRNA isoform expression
> and usage **Mol Brain** 17, 40 (2024). [doi:
> 10.1186/s13041-024-01112-7](https://doi.org/10.1186/s13041-024-01112-7)

And information about the paper, including code and publicly available
data can be found at this URL:

<https://github.com/lasseignelab/230227_EJ_MouseBrainIsoDiv>

In the abstract, Jones *et al.* describe the experiment:

> To assess differences in AS across the cerebellum, cortex,
> hippocampus, and striatum by sex, we generated and analyzed Oxford
> Nanopore Technologies (ONT) long-read RNA sequencing (lrRNA-Seq)
> C57BL/6J mouse brain cDNA libraries. From \> 85 million reads that
> passed quality control metrics, we calculated differential gene
> expression (DGE), differential transcript expression (DTE), and
> differential transcript usage (DTU) across brain regions and by sex.

``` r

library(readr)
# load DTU results 
dir <- system.file("extdata", package="splicelogic")
dtu_table <- readr::read_delim(file.path(dir, "dtu_table.tsv"))
```

    ## Rows: 49 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): tx_id, gene_id, gene_name
    ## dbl (2): estimate, padj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r

library(plyranges)
```

    ## Loading required package: BiocGenerics
    ## Loading required package: generics
    ## 
    ## Attaching package: 'generics'
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union
    ## 
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min
    ## 
    ## Loading required package: IRanges
    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: GenomicRanges
    ## Loading required package: Seqinfo
    ## Loading required package: dplyr
    ## 
    ## Attaching package: 'dplyr'
    ## 
    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following object is masked from 'package:Seqinfo':
    ## 
    ##     intersect
    ## 
    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union
    ## 
    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union
    ## 
    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union
    ## 
    ## The following object is masked from 'package:generics':
    ## 
    ##     explain
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union
    ## 
    ## 
    ## Attaching package: 'plyranges'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, n, n_distinct

``` r

# load exon ranges and metadata, and set seqinfo
exons <- read_bed(file.path(dir, "exons_M31.bed.gz"))
mcols(exons) <- DataFrame(readr::read_delim(file.path(dir, "exons_mcols.tsv.gz")))
```

    ## Rows: 601 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): exon_name, tx_id
    ## dbl (2): exon_id, exon_rank
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r

# M31 is paired with GRCm39 / mm39
si <- Seqinfo::Seqinfo(genome="mm39")
seqlevels(exons) <- seqlevels(si)
seqinfo(exons) <- si
```

next we join the …

``` r

# insert dtu results into exon metadata
txp_idx <- match(exons$tx_id, dtu_table$tx_id)
add_columns <- dtu_table[txp_idx, ] |>
  dplyr::select(-c(tx_id))
merged_DF <- cbind(mcols(exons), add_columns)
mcols(exons) <- merged_DF
```

As of Bioc 3.23 and *plyranges* version 1.32, the above can be replaced
with:

``` r

exons <- exons |> 
  plyranges::join_mcols_left(dtu_table, by = "tx_id")
```

splicelogic works on exons from significantly changed transcripts, so we
first filter out transcripts that were not signficant:

``` r

exons <- exons |> filter(padj < .05)
```

``` r

library(splicelogic)

processed_exons <- preprocess_input(exons, coef_col = "estimate")
calc_all_events(processed_exons)
```

    ## GRanges object with 16 ranges and 10 metadata columns:
    ##        seqnames              ranges strand |   exon_id            exon_name
    ##           <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##    [1]    chr17   66647479-66647535      - |    480827 ENSMUSE00000443570.7
    ##    [2]    chr14   20517526-20517591      - |    408079 ENSMUSE00001423050.2
    ##    [3]     chr8 120884207-120884236      + |    250989 ENSMUSE00001243257.2
    ##    [4]    chr12   91799829-91799996      - |    374021 ENSMUSE00001304078.2
    ##    [5]    chr12   91798541-91798558      - |    374018 ENSMUSE00001473756.2
    ##    ...      ...                 ...    ... .       ...                  ...
    ##   [12]     chr8 112437109-112438026      - |    262490 ENSMUSE00001391518.2
    ##   [13]     chr4 101504990-101505022      + |    107521 ENSMUSE00000631777.3
    ##   [14]     chr4 101513375-101513492      + |    107524 ENSMUSE00000671573.2
    ##   [15]     chr9   21849570-21849860      + |    266124 ENSMUSE00001334761.2
    ##   [16]    chr17   66643977-66645149      - |    480823 ENSMUSE00000791759.2
    ##        exon_rank                 tx_id               gene_id     gene_name
    ##        <numeric>           <character>           <character>   <character>
    ##    [1]        14 ENSMUST00000097291.10 ENSMUSG00000052105.18         Mtcl1
    ##    [2]         6  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##    [3]        10 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##    [4]         4 ENSMUST00000021347.12 ENSMUSG00000020964.15         Sel1l
    ##    [5]         4  ENSMUST00000178462.8 ENSMUSG00000020964.15         Sel1l
    ##    ...       ...                   ...                   ...           ...
    ##   [12]         7  ENSMUST00000212349.2 ENSMUSG00000031955.11         Bcar1
    ##   [13]         1  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [14]         3  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [15]         1  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [16]        14 ENSMUST00000086693.12 ENSMUSG00000052105.18         Mtcl1
    ##         estimate        padj              event              tx_event
    ##        <numeric>   <numeric>        <character>           <character>
    ##    [1]  -2.97320 0.009701213       skipped_exon ENSMUST00000086693.12
    ##    [2]   3.02406 0.018472554      included_exon ENSMUST00000065504.17
    ##    [3]   4.14230 0.000995041      included_exon  ENSMUST00000108951.8
    ##    [4]  -3.28535 0.001719345 mutually_exclusive  ENSMUST00000178462.8
    ##    [5]   2.99384 0.001719345 mutually_exclusive ENSMUST00000021347.12
    ##    ...       ...         ...                ...                   ...
    ##   [12]   3.85659  0.02700941               a3ss  ENSMUST00000166232.4
    ##   [13]   9.13055  0.01847255               a3ss ENSMUST00000030254.15
    ##   [14]   9.13055  0.01847255               a3ss ENSMUST00000030254.15
    ##   [15]   6.33671  0.00607823               a3ss ENSMUST00000046371.13
    ##   [16]   2.88524  0.01396713               a3ss ENSMUST00000097291.10
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Obtaining exon ranges with prepare_exons

[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
extracts exon ranges from a *TxDb* object and merges them with your DTU
results table in a single call. It returns a flat *GRanges* ready for
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md).

``` r

library(AnnotationHub)
ah <- AnnotationHub()
txdb <- ah[["AH75191"]]
# extract transcript table as tibble
txps <- txdb |>
  AnnotationDbi::select(keys(txdb, "TXID"), c("TXNAME","GENEID"), "TXID") |>
  tibble::as_tibble() |>
  dplyr::select(tx_num = TXID, tx_id = TXNAME, gene_id = GENEID)
# simulate DTU results
dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )

exons <- prepare_exons(txdb, dtu_table, coef_col = "effect_est", verbose = TRUE)
gr <- preprocess_input(exons, coef_col = "effect_est")
```

For a step-by-step breakdown of what
[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
does internally, see [Obtaining exon ranges
manually](#obtaining-exon-ranges-manually).

### Upstream methods

Something about what methods can be used for upstream DTU or switching
analysis.

## How to start detecting splicing events

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

### Skipped exons

[`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
identifies exons that are present in downregulated transcripts but
absent (skipped) in upregulated ones.

Using the gr object generated by
[`create_mock_data()`](https://thelovelab.github.io/splicelogic/reference/create_mock_data.md),
we can generate skiped exon events with the
[`generate_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/generate_events.md)
function, which modifies the input `GRanges` object to include skipped
exon events. Then we can run the detection function
[`calc_skipped_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
on the modified `GRanges` object:

``` r

gr_se <- generate_skipped_exons(gr, n_se = 2)
gr_se <- preprocess_input(gr_se, coef_col = "coefs")
se_result <- calc_skipped_exons(gr_se)
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

[`calc_mx_exons()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
detects pairs of exons where one exon is included in one isoform and a
different exon occupies the equivalent position in another isoform.

The
[`mx_mock_data()`](https://thelovelab.github.io/splicelogic/reference/mx_mock_data.md)
dataset is designed to contain mutually exclusive exon events:

``` r

gr_mx <- mx_mock_data()
gr_mx <- preprocess_input(gr_mx, coef_col = "coefs")
mx_result <- calc_mx_exons(gr_mx)
mx_result
```

``` r

mcols(mx_result)[, c("gene_id", "tx_id", "exon_rank", "event")]
```

### Retained introns

[`calc_retained_introns()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
finds exons in upregulated transcripts that overlap intronic regions of
downregulated transcripts, indicating intron retention.

``` r

gr_ri <- create_mock_data(3, 6, 3)
gr_ri <- preprocess_input(gr_ri, coef_col = "coefs")
gr_ri <- generate_retained_introns(gr_ri, n_ri = 3)
ri_result <- calc_retained_introns(gr_ri)
ri_result
```

``` r

if (length(ri_result) > 0) {
    mcols(ri_result)[, c("gene_id", "tx_id", "exon_rank", "event")]
}
```

### Alternative 3’ and 5’ splice sites

**NOTE: I think it’s more common to list 5 then 3**

[`calc_a3ss()`](https://thelovelab.github.io/splicelogic/reference/calc_events.md)
detects exons that share one boundary (start or end) with an exon in
another isoform but differ at the opposite end, indicating alternative
splice site usage.

``` r

gr_ss <- create_mock_data(3, 3, 6)
gr_ss <- preprocess_input(gr_ss, coef_col = "coefs")
gr_ss <- generate_a3ss(gr_ss, n_a3ss = 3)
ss_result <- calc_a3ss(gr_ss)
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

calc_skipped_exons(gr_none)
calc_mx_exons(gr_none)
```

## Obtaining exon ranges manually

This section walks through the steps that
[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
performs internally. This is useful if you need more control over the
process or want to understand how exon ranges are built from a *TxDb*.

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

Suppose we already have DTU results for transcripts:

``` r

# here just simulated results
dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )
```

``` r

# check the concordance...
all.equal(dtu_table$tx_num, as.integer(names(ebt)))
# and if aligned, set new names for exons
exons <- ebt
names(exons) <- dtu_table$tx_id
```

Next flattening the exons:

``` r

exons <- unlist(exons)
exons$tx_id <- names(exons)
names(exons) <- exons$exon_name
```

Adding DTU results and gene ID:

``` r

# arrange `dtu_table` as in `exons`, including duplicates
txp_idx <- match(exons$tx_id, dtu_table$tx_id)
add_columns <- dtu_table[txp_idx,] |>
  dplyr::select(-c(tx_id, tx_num))
merged_DF <- cbind(mcols(exons), add_columns)
mcols(exons) <- merged_DF
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
    ##  [1] splicelogic_0.0.67   plyranges_1.30.1     dplyr_1.2.0         
    ##  [4] GenomicRanges_1.62.1 Seqinfo_1.0.0        IRanges_2.44.0      
    ##  [7] S4Vectors_0.48.0     BiocGenerics_0.56.0  generics_0.1.4      
    ## [10] readr_2.2.0         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] SummarizedExperiment_1.40.0 rjson_0.2.23               
    ##  [3] xfun_0.56                   bslib_0.10.0               
    ##  [5] htmlwidgets_1.6.4           Biobase_2.70.0             
    ##  [7] lattice_0.22-9              tzdb_0.5.0                 
    ##  [9] vctrs_0.7.1                 tools_4.5.2                
    ## [11] bitops_1.0-9                curl_7.0.0                 
    ## [13] parallel_4.5.2              tibble_3.3.1               
    ## [15] pkgconfig_2.0.3             Matrix_1.7-4               
    ## [17] desc_1.4.3                  cigarillo_1.0.0            
    ## [19] lifecycle_1.0.5             compiler_4.5.2             
    ## [21] Rsamtools_2.26.0            textshaping_1.0.5          
    ## [23] Biostrings_2.78.0           codetools_0.2-20           
    ## [25] GenomeInfoDb_1.46.2         htmltools_0.5.9            
    ## [27] sass_0.4.10                 RCurl_1.98-1.17            
    ## [29] yaml_2.3.12                 pillar_1.11.1              
    ## [31] pkgdown_2.2.0               crayon_1.5.3               
    ## [33] jquerylib_0.1.4             BiocParallel_1.44.0        
    ## [35] cachem_1.1.0                DelayedArray_0.36.0        
    ## [37] abind_1.4-8                 tidyselect_1.2.1           
    ## [39] digest_0.6.39               restfulr_0.0.16            
    ## [41] fastmap_1.2.0               grid_4.5.2                 
    ## [43] cli_3.6.5                   SparseArray_1.10.9         
    ## [45] magrittr_2.0.4              S4Arrays_1.10.1            
    ## [47] XML_3.99-0.22               withr_3.0.2                
    ## [49] UCSC.utils_1.6.1            bit64_4.6.0-1              
    ## [51] httr_1.4.8                  rmarkdown_2.30             
    ## [53] XVector_0.50.0              matrixStats_1.5.0          
    ## [55] bit_4.6.0                   otel_0.2.0                 
    ## [57] ragg_1.5.1                  hms_1.1.4                  
    ## [59] evaluate_1.0.5              knitr_1.51                 
    ## [61] BiocIO_1.20.0               rtracklayer_1.70.1         
    ## [63] rlang_1.1.7                 glue_1.8.0                 
    ## [65] vroom_1.7.0                 jsonlite_2.0.0             
    ## [67] R6_2.6.1                    MatrixGenerics_1.22.0      
    ## [69] GenomicAlignments_1.46.0    systemfonts_1.3.2          
    ## [71] fs_1.6.7
