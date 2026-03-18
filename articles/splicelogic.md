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

`splicelogic` assumes the user has run a differential transcript usage
(DTU) or differential splicing analysis providing an error bound
(e.g. adjusted p-value / FDR) and and effect estimate with direction
(e.g. GLM estimated coefficient or deltaPSI) (see [upstream
methods](#upstream-methods)).

It also assumes the user has obtained genomic ranges representing the
exon structure of each transcript being analyzed (see [obtaining exon
ranges](#obtaining-exon-ranges)).

The `exons` should be provided in a flat *GRanges* object (one range per
exon), containing exon-level metadata in `mcols(exons)` such as the gene
ID, transcript ID, rank in the transcript, and the estimated change or
at least direction of change from the differential analysis.

Required columns for exons:

| Column | Description |
|----|----|
| `gene_id` | Gene identifier |
| `tx_id` | Transcript identifier |
| `exon_rank` | Exon rank within the transcript |
| `coef_col = "<user supplied>"` | differential effect estimate for the transcript |

The differential effect estimate column is specified by the user in
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md).
This should be the differential effect estimate associated with the
specific transcript containing the exon, and minimally could indicate
the direction of effect with `+1/-1`. Positive values indicate
up-regulated exons and negative values indicate down-regulated exons.
All exons from the same transcript will share the same value for this
column.

### Jones et al mouse long read dataset

We will use a published long read dataset and its reported splicing
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

We load the differential transcript usage analysis from the Jones *et
al.* paper. Specifically, we use the cortex-specific sex comparison (F
vs M) that identified the following transcripts as DTU:

<https://github.com/lasseignelab/230227_EJ_MouseBrainIsoDiv/blob/main/results/dtu_transcripts/cortex_sex_sig_features.csv>

The DTU table for this example, and the exons for the differential
transcripts have been saved in the `extdata` directory of *splicelogic*.

``` r

library(readr)
# load DTU results 
dir <- system.file("extdata", package="splicelogic")
dtu_table <- readr::read_delim(file.path(dir, "dtu_table.tsv"))
```

We load the 601 exons from GENCODE M31 annotation for the 49
differential transcripts, and attach the metadata columns.

``` r

library(plyranges)
exons_file <- "exons_M31.bed.gz"
exons_mcols_file <- "exons_mcols.tsv.gz"
exons <- plyranges::read_bed(file.path(dir, exons_file))
mcols(exons) <- DataFrame(readr::read_delim(file.path(dir, exons_mcols_file)))
```

Finally, we add `Seqinfo` for mm49, the reference mouse genome for
GENCODE M31.

``` r

si <- Seqinfo::Seqinfo(genome="mm39")
seqlevels(exons) <- seqlevels(si)
seqinfo(exons) <- si
```

Next, we join the DTU results onto the `exons` GRanges. This is
necessary for the *splicelogic* functions to know which transcripts are
up- or down-regulated, and to which gene the exons and transcripts
belong.

The code below matches the transcript IDs in the `exons` metadata with
the DTU results table, and then binds the two tables before assigning
back to `exons`.

``` r

txp_idx <- match(exons$tx_id, dtu_table$tx_id)
cols_to_add <- dtu_table[txp_idx, ] |> dplyr::select(-tx_id)
merged_DF <- S4Vectors::cbind(mcols(exons), cols_to_add)
mcols(exons) <- merged_DF
```

As of Bioc 3.23 and *plyranges* version 1.32, the above code to add
metadata columns from an additional table can be replaced with the
following:

``` r

exons <- exons |> 
  plyranges::join_mcols_left(dtu_table, by = "tx_id")
```

*splicelogic* is designed to take as input the exons from significantly
changed transcripts, so we first filter out transcripts that were not
signficant at FDR 10%.

``` r

exons <- exons |> filter(padj < .1)
```

## Detecting splicing events

### Preprocessing input data

The first step is to run
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md),
which prepares the `exons` for event detection. This function checks the
input data, ensures that the necessary columns are present:

``` r

library(splicelogic)
exons <- exons |> 
  preprocess_input(coef_col = "estimate")
```

### Individual events

Next we can run the various functions for calculating different types of
splicing events. E.g. to calculate exons which are skipped in
up-regulated transcripts relative to down-regulated transcripts, across
all genes:

**Skipped exons**

``` r

skipped <- exons |> calc_skipped_exons()
skipped
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

**Included exons**

``` r

included <- exons |> calc_included_exons()
included
```

    ## GRanges object with 3 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr14   20517526-20517591      - |    408079 ENSMUSE00001423050.2
    ##   [2]     chr1 163739641-163739706      + |     12966 ENSMUSE00000368805.4
    ##   [3]     chr8 120884207-120884236      + |    250989 ENSMUSE00001243257.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]         6  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [2]        20 ENSMUST00000077642.12 ENSMUSG00000026585.14        Kifap3
    ##   [3]        10 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##        estimate        padj         event              tx_event
    ##       <numeric>   <numeric>   <character>           <character>
    ##   [1]   3.02406 0.018472554 included_exon ENSMUST00000065504.17
    ##   [2]   1.04389 0.010395185 included_exon  ENSMUST00000027877.7
    ##   [3]   4.14230 0.000995041 included_exon  ENSMUST00000108951.8
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

**Mutually exclusive exons**

``` r

mx <- exons |> calc_mx_exons()
mx
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

**Retained introns**

``` r

ri <- exons |> calc_retained_introns()
ri
```

    ## GRanges object with 0 ranges and 0 metadata columns:
    ##    seqnames    ranges strand
    ##       <Rle> <IRanges>  <Rle>
    ##   -------
    ##   seqinfo: no sequences

**Alternative 5’ and 3’ splice sites**

``` r

a5ss <- exons |> calc_a5ss()
a5ss
```

    ## GRanges object with 5 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr10   88081618-88081868      + |    304334 ENSMUSE00001310024.2
    ##   [2]    chr14   20529963-20530189      - |    408097 ENSMUSE00000901772.3
    ##   [3]     chr8 120887954-120892045      + |    250998 ENSMUSE00000446870.6
    ##   [4]     chr9   21858242-21858348      + |    266133 ENSMUSE00001322549.2
    ##   [5]     chr9   21858900-21860203      + |    266139 ENSMUSE00001327764.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]         7  ENSMUST00000182183.8 ENSMUSG00000020056.17        Washc3
    ##   [2]         1  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [3]        13 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##   [4]         7  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [5]         9  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##        estimate        padj       event              tx_event
    ##       <numeric>   <numeric> <character>           <character>
    ##   [1]   9.19059 0.086737576        a5ss ENSMUST00000020248.16
    ##   [2]   3.02406 0.018472554        a5ss ENSMUST00000065504.17
    ##   [3]   4.14230 0.000995041        a5ss  ENSMUST00000108951.8
    ##   [4]   6.33671 0.006078234        a5ss ENSMUST00000046371.13
    ##   [5]   6.33671 0.006078234        a5ss ENSMUST00000046371.13
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

``` r

a3ss <- exons |> calc_a3ss()
a3ss
```

    ## GRanges object with 9 ranges and 10 metadata columns:
    ##       seqnames              ranges strand |   exon_id            exon_name
    ##          <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##   [1]    chr10   88037014-88037154      + |    304312 ENSMUSE00001309977.2
    ##   [2]    chr10   88055115-88055222      + |    304327 ENSMUSE00001223513.2
    ##   [3]    chr14   20505329-20506669      - |    408059 ENSMUSE00000564348.6
    ##   [4]     chr8 120840891-120841056      + |    250964 ENSMUSE00000678589.2
    ##   [5]     chr8 112437109-112438026      - |    262490 ENSMUSE00001391518.2
    ##   [6]     chr4 101504990-101505022      + |    107521 ENSMUSE00000631777.3
    ##   [7]     chr4 101513375-101513492      + |    107524 ENSMUSE00000671573.2
    ##   [8]     chr9   21849570-21849860      + |    266124 ENSMUSE00001334761.2
    ##   [9]    chr17   66643977-66645149      - |    480823 ENSMUSE00000791759.2
    ##       exon_rank                 tx_id               gene_id     gene_name
    ##       <numeric>           <character>           <character>   <character>
    ##   [1]         1  ENSMUST00000182183.8 ENSMUSG00000020056.17        Washc3
    ##   [2]         5  ENSMUST00000182183.8 ENSMUSG00000020056.17        Washc3
    ##   [3]        14  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##   [4]         1 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##   [5]         7  ENSMUST00000212349.2 ENSMUSG00000031955.11         Bcar1
    ##   [6]         1  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [7]         3  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [8]         1  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [9]        14 ENSMUST00000086693.12 ENSMUSG00000052105.18         Mtcl1
    ##        estimate        padj       event              tx_event
    ##       <numeric>   <numeric> <character>           <character>
    ##   [1]   9.19059 0.086737576        a3ss ENSMUST00000020248.16
    ##   [2]   9.19059 0.086737576        a3ss ENSMUST00000020248.16
    ##   [3]   3.02406 0.018472554        a3ss ENSMUST00000065504.17
    ##   [4]   4.14230 0.000995041        a3ss  ENSMUST00000108951.8
    ##   [5]   3.85659 0.027009405        a3ss  ENSMUST00000166232.4
    ##   [6]   9.13055 0.018472554        a3ss ENSMUST00000030254.15
    ##   [7]   9.13055 0.018472554        a3ss ENSMUST00000030254.15
    ##   [8]   6.33671 0.006078234        a3ss ENSMUST00000046371.13
    ##   [9]   2.88524 0.013967132        a3ss ENSMUST00000097291.10
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

### Calculate all events

``` r

all_events <- exons |> calc_all_events()
```

    ## Calculating skipped exon events...

    ## Calculating included exon events...

    ## Calculating mutually exclusive exon events...

    ## Calculating retained intron events...

    ## Calculating alternative 5' splice site events...

    ## Calculating alternative 3' splice site events...

    ## Done! 20 events detected.

``` r

all_events
```

    ## GRanges object with 20 ranges and 10 metadata columns:
    ##        seqnames              ranges strand |   exon_id            exon_name
    ##           <Rle>           <IRanges>  <Rle> | <numeric>          <character>
    ##    [1]    chr17   66647479-66647535      - |    480827 ENSMUSE00000443570.7
    ##    [2]    chr14   20517526-20517591      - |    408079 ENSMUSE00001423050.2
    ##    [3]     chr1 163739641-163739706      + |     12966 ENSMUSE00000368805.4
    ##    [4]     chr8 120884207-120884236      + |    250989 ENSMUSE00001243257.2
    ##    [5]    chr12   91799829-91799996      - |    374021 ENSMUSE00001304078.2
    ##    ...      ...                 ...    ... .       ...                  ...
    ##   [16]     chr8 112437109-112438026      - |    262490 ENSMUSE00001391518.2
    ##   [17]     chr4 101504990-101505022      + |    107521 ENSMUSE00000631777.3
    ##   [18]     chr4 101513375-101513492      + |    107524 ENSMUSE00000671573.2
    ##   [19]     chr9   21849570-21849860      + |    266124 ENSMUSE00001334761.2
    ##   [20]    chr17   66643977-66645149      - |    480823 ENSMUSE00000791759.2
    ##        exon_rank                 tx_id               gene_id     gene_name
    ##        <numeric>           <character>           <character>   <character>
    ##    [1]        14 ENSMUST00000097291.10 ENSMUSG00000052105.18         Mtcl1
    ##    [2]         6  ENSMUST00000100844.6 ENSMUSG00000021814.18         Anxa7
    ##    [3]        20 ENSMUST00000077642.12 ENSMUSG00000026585.14        Kifap3
    ##    [4]        10 ENSMUST00000034281.13 ENSMUSG00000031824.16 6430548M08Rik
    ##    [5]         4 ENSMUST00000021347.12 ENSMUSG00000020964.15         Sel1l
    ##    ...       ...                   ...                   ...           ...
    ##   [16]         7  ENSMUST00000212349.2 ENSMUSG00000031955.11         Bcar1
    ##   [17]         1  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [18]         3  ENSMUST00000106927.2 ENSMUSG00000035212.15        Leprot
    ##   [19]         1  ENSMUST00000190387.7 ENSMUSG00000040563.14        Plppr2
    ##   [20]        14 ENSMUST00000086693.12 ENSMUSG00000052105.18         Mtcl1
    ##         estimate        padj              event              tx_event
    ##        <numeric>   <numeric>        <character>           <character>
    ##    [1]  -2.97320 0.009701213       skipped_exon ENSMUST00000086693.12
    ##    [2]   3.02406 0.018472554      included_exon ENSMUST00000065504.17
    ##    [3]   1.04389 0.010395185      included_exon  ENSMUST00000027877.7
    ##    [4]   4.14230 0.000995041      included_exon  ENSMUST00000108951.8
    ##    [5]  -3.28535 0.001719345 mutually_exclusive  ENSMUST00000178462.8
    ##    ...       ...         ...                ...                   ...
    ##   [16]   3.85659  0.02700941               a3ss  ENSMUST00000166232.4
    ##   [17]   9.13055  0.01847255               a3ss ENSMUST00000030254.15
    ##   [18]   9.13055  0.01847255               a3ss ENSMUST00000030254.15
    ##   [19]   6.33671  0.00607823               a3ss ENSMUST00000046371.13
    ##   [20]   2.88524  0.01396713               a3ss ENSMUST00000097291.10
    ##   -------
    ##   seqinfo: 61 sequences (1 circular) from mm39 genome

``` r

barplot(table(all_events$event))
```

![](splicelogic_files/figure-html/unnamed-chunk-3-1.png)

## Upstream methods

Something about what methods can be used for upstream DTU or switching
analysis.

## Obtaining exon ranges with prepare_exons

[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
extracts exon ranges from a *TxDb* object and merges them with your DTU
results table in a single call. It returns a flat *GRanges* ready for
[`preprocess_input()`](https://thelovelab.github.io/splicelogic/reference/preprocess_input.md).
We demonstrate using GENCODE v32 (human genes).

The first step is to load a *TxDb* object. Typically, a user would
supply their own GTF to txdbmaker::makeTxDbFromGFF() to generate this.
For this demonstration we will load a pre-constructed *TxDb* from
Bioconductor’s AnnotationHub.

``` r

suppressPackageStartupMessages({
  library(AnnotationHub)
})
ah <- AnnotationHub()
txdb <- ah[["AH75191"]] # GENCODE v32 (human)
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

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(tibble)
})
txps <- txdb |>
  AnnotationDbi::select(keys(txdb, "TXID"), c("TXNAME","GENEID"), "TXID") |>
  tibble::as_tibble() |>
  dplyr::select(tx_num = TXID, tx_id = TXNAME, gene_id = GENEID)
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r

txps
```

    ## # A tibble: 227,462 × 3
    ##    tx_num tx_id             gene_id          
    ##     <int> <chr>             <chr>            
    ##  1      1 ENST00000456328.2 ENSG00000223972.5
    ##  2      2 ENST00000450305.2 ENSG00000223972.5
    ##  3      3 ENST00000473358.1 ENSG00000243485.5
    ##  4      4 ENST00000469289.1 ENSG00000243485.5
    ##  5      5 ENST00000607096.1 ENSG00000284332.1
    ##  6      6 ENST00000606857.1 ENSG00000268020.3
    ##  7      7 ENST00000642116.1 ENSG00000240361.2
    ##  8      8 ENST00000492842.2 ENSG00000240361.2
    ##  9      9 ENST00000641515.2 ENSG00000186092.6
    ## 10     10 ENST00000335137.4 ENSG00000186092.6
    ## # ℹ 227,452 more rows

``` r

# simulate DTU results
dtu_table <- txps |>
  dplyr::mutate(
    padj = runif(dplyr::n()),
    effect_est = rnorm(dplyr::n())
  )
```

``` r

exons <- prepare_exons(txdb, dtu_table, coef_col = "effect_est", verbose = TRUE)
```

    ## Extracting exons from TxDb...

    ## Mapping transcript IDs...

    ## Merging DTU results onto exons...

    ## Done. Returned 1372308 exon ranges from 227462unique transcripts.

``` r

exons |> preprocess_input(coef_col = "effect_est")
```

    ## GRanges object with 1372308 ranges and 12 metadata columns:
    ##                     seqnames      ranges strand |   exon_id         exon_name
    ##                        <Rle>   <IRanges>  <Rle> | <integer>       <character>
    ##   ENSE00002234944.1     chr1 11869-12227      + |         1 ENSE00002234944.1
    ##   ENSE00003582793.1     chr1 12613-12721      + |         5 ENSE00003582793.1
    ##   ENSE00002312635.1     chr1 13221-14409      + |         8 ENSE00002312635.1
    ##   ENSE00001948541.1     chr1 12010-12057      + |         2 ENSE00001948541.1
    ##   ENSE00001671638.2     chr1 12179-12227      + |         3 ENSE00001671638.2
    ##                 ...      ...         ...    ... .       ...               ...
    ##   ENSE00001544488.1     chrM   5826-5891      - |    745885 ENSE00001544488.1
    ##   ENSE00001544487.2     chrM   7446-7514      - |    745886 ENSE00001544487.2
    ##   ENSE00001434974.2     chrM 14149-14673      - |    745887 ENSE00001434974.2
    ##   ENSE00001544476.1     chrM 14674-14742      - |    745888 ENSE00001544476.1
    ##   ENSE00001544473.2     chrM 15956-16023      - |    745889 ENSE00001544473.2
    ##                     exon_rank             tx_id    tx_num           gene_id
    ##                     <integer>       <character> <integer>       <character>
    ##   ENSE00002234944.1         1 ENST00000456328.2         1 ENSG00000223972.5
    ##   ENSE00003582793.1         2 ENST00000456328.2         1 ENSG00000223972.5
    ##   ENSE00002312635.1         3 ENST00000456328.2         1 ENSG00000223972.5
    ##   ENSE00001948541.1         1 ENST00000450305.2         2 ENSG00000223972.5
    ##   ENSE00001671638.2         2 ENST00000450305.2         2 ENSG00000223972.5
    ##                 ...       ...               ...       ...               ...
    ##   ENSE00001544488.1         1 ENST00000387409.1    227458 ENSG00000210144.1
    ##   ENSE00001544487.2         1 ENST00000387416.2    227459 ENSG00000210151.2
    ##   ENSE00001434974.2         1 ENST00000361681.2    227460 ENSG00000198695.2
    ##   ENSE00001544476.1         1 ENST00000387459.1    227461 ENSG00000210194.1
    ##   ENSE00001544473.2         1 ENST00000387461.2    227462 ENSG00000210196.2
    ##                          padj effect_est                 key    nexons
    ##                     <numeric>  <numeric>         <character> <integer>
    ##   ENSE00002234944.1 0.0807501  -0.496625 ENST00000456328.2-1         3
    ##   ENSE00003582793.1 0.0807501  -0.496625 ENST00000456328.2-2         3
    ##   ENSE00002312635.1 0.0807501  -0.496625 ENST00000456328.2-3         3
    ##   ENSE00001948541.1 0.8343330   0.938183 ENST00000450305.2-1         6
    ##   ENSE00001671638.2 0.8343330   0.938183 ENST00000450305.2-2         6
    ##                 ...       ...        ...                 ...       ...
    ##   ENSE00001544488.1  0.917140   0.886557 ENST00000387409.1-1         1
    ##   ENSE00001544487.2  0.712942  -0.383043 ENST00000387416.2-1         1
    ##   ENSE00001434974.2  0.746749  -2.515144 ENST00000361681.2-1         1
    ##   ENSE00001544476.1  0.118248  -1.629253 ENST00000387459.1-1         1
    ##   ENSE00001544473.2  0.398895   0.376717 ENST00000387461.2-1         1
    ##                      internal estimates
    ##                     <logical> <numeric>
    ##   ENSE00002234944.1     FALSE -0.496625
    ##   ENSE00003582793.1      TRUE -0.496625
    ##   ENSE00002312635.1     FALSE -0.496625
    ##   ENSE00001948541.1     FALSE  0.938183
    ##   ENSE00001671638.2      TRUE  0.938183
    ##                 ...       ...       ...
    ##   ENSE00001544488.1     FALSE  0.886557
    ##   ENSE00001544487.2     FALSE -0.383043
    ##   ENSE00001434974.2     FALSE -2.515144
    ##   ENSE00001544476.1     FALSE -1.629253
    ##   ENSE00001544473.2     FALSE  0.376717
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

``` r

# followed by calculation of splicing events
```

The following section details what
[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
does internally.

## Obtaining exon ranges manually

This section walks through the steps that
[`prepare_exons()`](https://thelovelab.github.io/splicelogic/reference/prepare_exons.md)
performs internally. This is useful if you need more control over the
process or want to understand how exon ranges are built from a *TxDb*.

The following extracts the exons grouped by transcript from the *TxDb*:

``` r

suppressPackageStartupMessages({
  library(GenomicFeatures)
})
# extract exons as a GRangesList
exons_list <- GenomicFeatures::exonsBy(txdb, by="tx") # exon id, name, and rank
# dtu_table aligns with txps, which aligns with the names of the GRangesList.
# prepare_exons() handles this with alignment checks
names(exons_list) <- dtu_table$tx_id
```

Next flattening the exons:

``` r

exons <- unlist(exons_list)
# swap tx_id with exon_name as the names of the GRanges
exons$tx_id <- names(exons) # store transcript ids
names(exons) <- exons$exon_name
```

Adding DTU results and gene ID:

``` r

txp_idx <- match(exons$tx_id, dtu_table$tx_id)
cols_to_add <- dtu_table[txp_idx,] |> dplyr::select(-c(tx_id, tx_num))
merged_DF <- cbind(mcols(exons), cols_to_add)
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
    ##  [7] dbplyr_2.5.2           splicelogic_0.99.0     plyranges_1.30.1      
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
    ## [19] utf8_1.2.6                  yaml_2.3.12                
    ## [21] rtracklayer_1.70.1          knitr_1.51                 
    ## [23] S4Arrays_1.10.1             htmlwidgets_1.6.4          
    ## [25] bit_4.6.0                   curl_7.0.0                 
    ## [27] DelayedArray_0.36.0         abind_1.4-8                
    ## [29] BiocParallel_1.44.0         purrr_1.2.1                
    ## [31] withr_3.0.2                 desc_1.4.3                 
    ## [33] grid_4.5.2                  SummarizedExperiment_1.40.0
    ## [35] cli_3.6.5                   rmarkdown_2.30             
    ## [37] crayon_1.5.3                ragg_1.5.1                 
    ## [39] otel_0.2.0                  httr_1.4.8                 
    ## [41] tzdb_0.5.0                  rjson_0.2.23               
    ## [43] DBI_1.3.0                   cachem_1.1.0               
    ## [45] parallel_4.5.2              BiocManager_1.30.27        
    ## [47] XVector_0.50.0              restfulr_0.0.16            
    ## [49] matrixStats_1.5.0           vctrs_0.7.1                
    ## [51] Matrix_1.7-4                jsonlite_2.0.0             
    ## [53] hms_1.1.4                   bit64_4.6.0-1              
    ## [55] systemfonts_1.3.2           jquerylib_0.1.4            
    ## [57] glue_1.8.0                  pkgdown_2.2.0              
    ## [59] codetools_0.2-20            BiocVersion_3.22.0         
    ## [61] GenomeInfoDb_1.46.2         BiocIO_1.20.0              
    ## [63] UCSC.utils_1.6.1            pillar_1.11.1              
    ## [65] rappdirs_0.3.4              htmltools_0.5.9            
    ## [67] R6_2.6.1                    httr2_1.2.2                
    ## [69] textshaping_1.0.5           vroom_1.7.0                
    ## [71] evaluate_1.0.5              lattice_0.22-9             
    ## [73] png_0.1-9                   Rsamtools_2.26.0           
    ## [75] cigarillo_1.0.0             memoise_2.0.1              
    ## [77] bslib_0.10.0                SparseArray_1.10.9         
    ## [79] xfun_0.56                   fs_1.6.7                   
    ## [81] MatrixGenerics_1.22.0       pkgconfig_2.0.3
