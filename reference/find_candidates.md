# Find candidates skipped exons

computes similarity between positive and negative GenomicRanges and then
finds candidate exons in neg whose adjacent exons exist in the positive
set and are also adjacent

## Usage

``` r
find_candidates(gr_pos, gr_neg, id_pos, id_neg)
```

## Examples

``` r
gr <- se_mock_data()
library(ggplot2)
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
library(plyranges)
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:GenomicRanges’:
#> 
#>     intersect, setdiff, union
#> The following object is masked from ‘package:Seqinfo’:
#> 
#>     intersect
#> The following objects are masked from ‘package:IRanges’:
#> 
#>     collapse, desc, intersect, setdiff, slice, union
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     first, intersect, rename, setdiff, setequal, union
#> The following objects are masked from ‘package:BiocGenerics’:
#> 
#>     combine, intersect, setdiff, setequal, union
#> The following object is masked from ‘package:generics’:
#> 
#>     explain
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: ‘plyranges’
#> The following objects are masked from ‘package:dplyr’:
#> 
#>     between, n, n_distinct
plot_gr(gr)
#> Error in plot_gr(gr): could not find function "plot_gr"
gr_pos <- filter(gr, coefs > 0)
gr_neg <- filter(gr, coefs < 0)
find_candidates(gr_pos, gr_neg, gr_pos$tx_id, gr_neg$tx_id)
#> GRanges object with 4 ranges and 6 metadata columns:
#>       seqnames    ranges strand | exon_rank   gene_id     coefs     tx_id
#>          <Rle> <IRanges>  <Rle> | <integer> <numeric> <numeric> <numeric>
#>   [1]     chr1     21-25      + |         3         1 -0.399239         1
#>   [2]     chr1     21-25      + |         3         1 -0.399239         1
#>   [3]     chr1     51-55      + |         6         1 -0.399239         1
#>   [4]     chr1     51-55      + |         6         1 -0.399239         1
#>        tx_event        event
#>       <numeric>  <character>
#>   [1]         2 skipped_exon
#>   [2]         3 skipped_exon
#>   [3]         2 skipped_exon
#>   [4]         3 skipped_exon
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
