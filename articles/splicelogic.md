# Introduction to splicelogic

## Introduction

`splicelogic` is an R/Bioconductor package for detecting alternative
splicing events from exon-level data stored as `GRanges` objects. Given
a set of exons annotated with a coefficient column indicating
differential transcript usage (DTU), `splicelogic` identifies the
following types of splicing events:

- **Skipped exons (SE)** – exons present in one isoform but absent in
  another.
- **Retained exons (RE)** – exons present in one isoform but absent in
  another.
- **Mutually exclusive exons (MXE)** – pairs of exons where one is
  included at the expense of the other.
- **Retained introns (RI)** – intronic regions that are retained as part
  of an exon in an alternative isoform.
- **Alternative 5’ (A5SS)** – exons that share the same 3’ splice site
  but differ at the 5’ splice site.
- **Alternative 3’ (A3SS)** – exons that share the same 5’ splice site
  but differ at the 3’ splice site.

## Quick start

With DTU results attached to a *GRanges* of the exons from significant
transcripts, one can use the following code to identify splice events:

``` r

exons <- preprocess_input(exons, coef_col = "estimates")
skipped <- exons |> calc_skipped_exons()
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

### Loading data

``` r

library(readr)
# load DTU results 
dir <- system.file("extdata", package="splicelogic")
dtu_table <- readr::read_delim(file.path(dir, "dtu_table.tsv"))

library(plyranges)
# load exon ranges and metadata, and set seqinfo
exons <- read_bed(file.path(dir, "exons_M31.bed.gz"))
mcols(exons) <- DataFrame(
  readr::read_delim(file.path(dir, "exons_mcols.tsv.gz"))
)
# M31 is paired with GRCm39 / mm39
si <- Seqinfo::Seqinfo(genome="mm39")
seqlevels(exons) <- seqlevels(si)
seqinfo(exons) <- si
```

Next, we join the DTU results into the exon GRanges metadata. This is
necessary for the splicelogic functions to know which transcripts are up
or downregulated, and to which gene they belong. The code below matches
the transcript IDs in the exon metadata with those in the DTU results
table, and adds the relevant columns from the DTU table to the exon
metadata.

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

## How to start detecting splicing events

`splicelogic` works on exons from significantly changed transcripts, so
we first filter out transcripts that were not signficant:

``` r

exons <- exons |> filter(padj < .05)
```

### Preprocessing input data

The first step is to run
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md),
which prepares the exon data for event detection. This function checks
the input data, ensures that the necessary

``` r

library(splicelogic)
processed_exons <- preprocess_input(exons, coef_col = "estimate")
```

### Calculate skipped exons

After preprocessing, you can run the various functions to calculate
different types of splicing events. To calculate skipped exons:

``` r

skipped_exons <- calc_skipped_exons(processed_exons)
skipped_exons
```

    ## GRanges object with 1 range and 10 metadata columns:
    ##       seqnames            ranges strand |   exon_id            exon_name
    ##          <Rle>         <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr17 66647479-66647535      - |    480827 ENSMUSE00000443570.7
    ##       exon_rank                 tx_id               gene_id   gene_name
    ##       <numeric>           <character>           <character> <character>
    ##   [1]        14 ENSMUST00000097291.10 ENSMUSG00000052105.18       Mtcl1
    ##        estimate       padj        event              tx_event
    ##       <numeric>  <numeric>  <character>           <character>
    ##   [1]   -2.9732 0.00970121 skipped_exon ENSMUST00000086693.12
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### calculate included exons

``` r

included_exons <- calc_included_exons(processed_exons)
included_exons
```

    ## GRanges object with 2 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr14   20517526-20517591      - |    408079 ENSMUSE00001423050.2
    ##   [2]     chr8 120884207-120884236      + |    250989 ENSMUSE00001243257.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]         6  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [2]        10 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##        estimate        padj         event              tx_event
    ##       <numeric>   <numeric>   <character>           <character>
    ##   [1]   3.02406 0.018472554 included_exon ENSMUST00000065504.17
    ##   [2]   4.14230 0.000995041 included_exon  ENSMUST00000108951.8
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Calculate mutually exclusive exons

``` r

mx_exons <- calc_mx_exons(processed_exons)
mx_exons
```

    ## GRanges object with 2 ranges and 10 metadata columns:
    ##       seqnames            ranges strand |   exon_id            exon_name
    ##          <Rle>         <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr12 91799829-91799996      - |    374021 ENSMUSE00001304078.2
    ##   [2]    chr12 91798541-91798558      - |    374018 ENSMUSE00001473756.2
    ##       exon_rank                 tx_id               gene_id   gene_name
    ##       <numeric>           <character>           <character> <character>
    ##   [1]         4 ENSMUST00000021347.12 ENSMUSG00000020964.15       Sel1l
    ##   [2]         4  ENSMUST00000178462.8 ENSMUSG00000020964.15       Sel1l
    ##        estimate       padj              event              tx_event
    ##       <numeric>  <numeric>        <character>           <character>
    ##   [1]  -3.28535 0.00171934 mutually_exclusive  ENSMUST00000178462.8
    ##   [2]   2.99384 0.00171934 mutually_exclusive ENSMUST00000021347.12
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Calculate retained introns

``` r

retained_introns <- calc_retained_introns(processed_exons)
retained_introns
```

    ## GRanges object with 0 ranges and 0 metadata columns:
    ##    seqnames    ranges strand
    ##       <Rle> <IRanges>  <Rle>
    ##   -------
    ##   seqinfo: no sequences

### Calculate alternative 5’ splice sites

``` r

a5ss <- calc_a5ss(processed_exons)
a5ss
```

    ## GRanges object with 4 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr14   20529963-20530189      - |    408097 ENSMUSE00000901772.3
    ##   [2]     chr8 120887954-120892045      + |    250998 ENSMUSE00000446870.6
    ##   [3]     chr9   21858242-21858348      + |    266133 ENSMUSE00001322549.2
    ##   [4]     chr9   21858900-21860203      + |    266139 ENSMUSE00001327764.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]         1  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [2]        13 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##   [3]         7  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [4]         9  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##        estimate        padj       event              tx_event
    ##       <numeric>   <numeric> <character>           <character>
    ##   [1]   3.02406 0.018472554        a5ss ENSMUST00000065504.17
    ##   [2]   4.14230 0.000995041        a5ss  ENSMUST00000108951.8
    ##   [3]   6.33671 0.006078234        a5ss ENSMUST00000046371.13
    ##   [4]   6.33671 0.006078234        a5ss ENSMUST00000046371.13
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Calculate alternative 3’ splice sites

``` r

a3ss <- calc_a3ss(processed_exons)
a3ss
```

    ## GRanges object with 7 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr14   20505329-20506669      - |    408059 ENSMUSE00000564348.6
    ##   [2]     chr8 120840891-120841056      + |    250964 ENSMUSE00000678589.2
    ##   [3]     chr8 112437109-112438026      - |    262490 ENSMUSE00001391518.2
    ##   [4]     chr4 101504990-101505022      + |    107521 ENSMUSE00000631777.3
    ##   [5]     chr4 101513375-101513492      + |    107524 ENSMUSE00000671573.2
    ##   [6]     chr9   21849570-21849860      + |    266124 ENSMUSE00001334761.2
    ##   [7]    chr17   66643977-66645149      - |    480823 ENSMUSE00000791759.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]        14  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [2]         1 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##   [3]         7  ENSMUST00000212349.2 ENSMUSG00000031955.11         Bcar1
    ##   [4]         1  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [5]         3  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [6]         1  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [7]        14 ENSMUST00000086693.12 ENSMUSG00000052105.18         Mtcl1
    ##        estimate        padj       event              tx_event
    ##       <numeric>   <numeric> <character>           <character>
    ##   [1]   3.02406 0.018472554        a3ss ENSMUST00000065504.17
    ##   [2]   4.14230 0.000995041        a3ss  ENSMUST00000108951.8
    ##   [3]   3.85659 0.027009405        a3ss  ENSMUST00000166232.4
    ##   [4]   9.13055 0.018472554        a3ss ENSMUST00000030254.15
    ##   [5]   9.13055 0.018472554        a3ss ENSMUST00000030254.15
    ##   [6]   6.33671 0.006078234        a3ss ENSMUST00000046371.13
    ##   [7]   2.88524 0.013967132        a3ss ENSMUST00000097291.10
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Calculate all events

``` r

all_events <- calc_all_events(processed_exons)
```

    ## Calculating skipped exon events...

    ## Calculating included exon events...

    ## Calculating mutually exclusive exon events...

    ## Calculating retained intron events...

    ## Calculating alternative 5' splice site events...

    ## Calculating alternative 3' splice site events...

    ## Done! 16 events detected.

``` r

all_events
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

## Upstream methods

Something about what methods can be used for upstream DTU or switching
analysis.

## Obtaining exon ranges with prepare_exons

[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
extracts exon ranges from a *TxDb* object and merges them with your DTU
results table in a single call. It returns a flat *GRanges* ready for
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md).

``` r

library(AnnotationHub)
```

    ## Loading required package: BiocFileCache

    ## Loading required package: dbplyr

    ## 
    ## Attaching package: 'dbplyr'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     ident, sql

``` r

ah <- AnnotationHub()
txdb <- ah[["AH75191"]]
```

    ## loading from cache

    ## Loading required package: GenomicFeatures

    ## Loading required package: AnnotationDbi

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:AnnotationHub':
    ## 
    ##     cache

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r

# extract transcript table as tibble
txps <- txdb |>
  AnnotationDbi::select(keys(txdb, "TXID"), c("TXNAME","GENEID"), "TXID") |>
  tibble::as_tibble() |>
  dplyr::select(tx_num = TXID, tx_id = TXNAME, gene_id = GENEID)
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r

# simulate DTU results
dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )

exons <- prepare_exons(txdb, dtu_table, coef_col = "effect_est", verbose = TRUE)
```

    ## Extracting exons from TxDb...

    ## Mapping transcript IDs...

    ## Merging DTU results onto exons...

    ## Done. Returned 1372308 exon ranges from 227462unique transcripts.

``` r

gr <- preprocess_input(exons, coef_col = "effect_est")
```

For a step-by-step breakdown of what
[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
does internally, see [Obtaining exon ranges
manually](#obtaining-exon-ranges-manually).

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
dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )
```

``` r

# check the concordance...
all.equal(dtu_table$tx_num, as.integer(names(ebt)))
```

    ## [1] TRUE

``` r

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
    ##  [1] tibble_3.3.1           GenomicFeatures_1.62.0 AnnotationDbi_1.72.0  
    ##  [4] Biobase_2.70.0         AnnotationHub_4.0.0    BiocFileCache_3.0.0   
    ##  [7] dbplyr_2.5.2           splicelogic_0.0.67     plyranges_1.30.1      
    ## [10] dplyr_1.2.0            GenomicRanges_1.62.1   Seqinfo_1.0.0         
    ## [13] IRanges_2.44.0         S4Vectors_0.48.0       BiocGenerics_0.56.0   
    ## [16] generics_0.1.4         readr_2.2.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1            blob_1.3.0                 
    ##  [3] filelock_1.0.3              Biostrings_2.78.0          
    ##  [5] bitops_1.0-9                fastmap_1.2.0              
    ##  [7] RCurl_1.98-1.17             GenomicAlignments_1.46.0   
    ##  [9] XML_3.99-0.22               digest_0.6.39              
    ## [11] lifecycle_1.0.5             KEGGREST_1.50.0            
    ## [13] RSQLite_2.4.6               magrittr_2.0.4             
    ## [15] compiler_4.5.2              rlang_1.1.7                
    ## [17] sass_0.4.10                 tools_4.5.2                
    ## [19] yaml_2.3.12                 rtracklayer_1.70.1         
    ## [21] knitr_1.51                  S4Arrays_1.10.1            
    ## [23] htmlwidgets_1.6.4           bit_4.6.0                  
    ## [25] curl_7.0.0                  DelayedArray_0.36.0        
    ## [27] abind_1.4-8                 BiocParallel_1.44.0        
    ## [29] purrr_1.2.1                 withr_3.0.2                
    ## [31] desc_1.4.3                  grid_4.5.2                 
    ## [33] SummarizedExperiment_1.40.0 cli_3.6.5                  
    ## [35] rmarkdown_2.30              crayon_1.5.3               
    ## [37] ragg_1.5.1                  otel_0.2.0                 
    ## [39] httr_1.4.8                  tzdb_0.5.0                 
    ## [41] rjson_0.2.23                DBI_1.3.0                  
    ## [43] cachem_1.1.0                parallel_4.5.2             
    ## [45] BiocManager_1.30.27         XVector_0.50.0             
    ## [47] restfulr_0.0.16             matrixStats_1.5.0          
    ## [49] vctrs_0.7.1                 Matrix_1.7-4               
    ## [51] jsonlite_2.0.0              hms_1.1.4                  
    ## [53] bit64_4.6.0-1               systemfonts_1.3.2          
    ## [55] jquerylib_0.1.4             glue_1.8.0                 
    ## [57] pkgdown_2.2.0               codetools_0.2-20           
    ## [59] BiocVersion_3.22.0          GenomeInfoDb_1.46.2        
    ## [61] BiocIO_1.20.0               UCSC.utils_1.6.1           
    ## [63] pillar_1.11.1               rappdirs_0.3.4             
    ## [65] htmltools_0.5.9             R6_2.6.1                   
    ## [67] httr2_1.2.2                 textshaping_1.0.5          
    ## [69] vroom_1.7.0                 evaluate_1.0.5             
    ## [71] lattice_0.22-9              png_0.1-9                  
    ## [73] Rsamtools_2.26.0            cigarillo_1.0.0            
    ## [75] memoise_2.0.1               bslib_0.10.0               
    ## [77] SparseArray_1.10.9          xfun_0.56                  
    ## [79] fs_1.6.7                    MatrixGenerics_1.22.0      
    ## [81] pkgconfig_2.0.3
