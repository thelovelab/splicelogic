# splicelogic: splicing events from transcript sets

*splicelogic* turns sets of transcripts into discrete splicing events.
Unlike event-based tools that work at the junction level, *splicelogic*
operates on whole transcript structures: within each gene it compares
two groups of transcripts and all of their exons, so events are derived
with full isoform context. It detects skipped exons, included exons,
mutually exclusive exons, retained introns, alternative 5’ and 3’ splice
sites, and alternative transcription start and end sites.

*splicelogic* allows the two groups of transcripts to be compared to be
defined upstream in one of two ways. They can come from an explicit
partition of the transcripts into two sets, e.g. a set of novel
transcripts against a set of annotated reference transcripts, or
transcripts found in disease samples relative to those found in control
samples. Or they can come from a differential transcript usage (DTU)
analysis, which gives each transcript an effect estimate whose sign
defines the comparison. Because it takes transcript-level results tables
as input, *splicelogic* is compatible with any upstream DTU method
(including DRIMSeq, DEXSeq, satuRn, and edgeR), supporting flexible
experimental designs.

*splicelogic* operates on exon-level data stored as *GRanges* objects
within R/Bioconductor. Given a set of exons carrying the gene ID,
transcript ID and exon rank, *splicelogic* can be used to identify a
variety of splicing events. See the
[vignette](https://thelovelab.github.io/splicelogic/articles/splicelogic.html)
for more details. *splicelogic* works well within the
[tidyomics](https://github.com/tidyomics/) framework. Event results
including the exons and introns underlying each event can be further
manipulated downstream using
[plyranges](https://tidyomics.github.io/plyranges).

# How to install

`splicelogic` is [available from
Bioconductor](https://bioconductor.org/packages/splicelogic):

    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("splicelogic")

The most recent development version can be installed from the `devel`
branch on GitHub:

    devtools::install_github("thelovelab/splicelogic")

# Quick start

## Finding events from transcript sets

    # prepare exons from two sets of transcripts
    exons <- prepare_exons_by_partition(
      up = <GRanges OF EXONS, OR tx IDs>,
      down = <GRanges OF EXONS, OR tx IDs>
    ) |>
      preprocess(coef_col = "estimate")

    # find skipped exons
    skipped <- exons |> find_se()

    # find all splicing events
    all_events <- exons |> find_all_events()

## Finding events from DTU results

    # prepare exons from a TxDb and DTU results
    exons <- prepare_exons(
      txdb = <A TxDb OBJECT>,
      dtu_table = <DTU_TABLE>,
      coef_col = "estimate"
    )

    # preprocess for further analysis
    exons <- preprocess(exons, coef_col = "estimate")

    # find skipped exons
    skipped <- exons |> find_se()

    # find all splicing events
    all_events <- exons |> find_all_events()

# Future directions

- Support detection of additional event types, such as consecutive
  skipped exons or loss of retained introns.
- Extraction and labelling of the specific splice junctions associated
  with each event, adding metadata columns such as the donor–acceptor
  dinucleotide sequence (e.g. AG-GT) and a logical indicating whether
  the junction is canonical, for downstream interpretation.
- Custom plotting functions to easily visualize the results, showing the
  exons involved in each event alongside the transcripts they come from
  and the event annotation.
- Support detection of alternative UTR events (alternative 5’ and 3’
  UTRs), when the reference annotation includes UTR coordinates
  (e.g. GENCODE), building on the alternative transcription start and
  end sites reported by
  [`find_atss()`](https://thelovelab.github.io/splicelogic/reference/find_events.md)
  and
  [`find_ates()`](https://thelovelab.github.io/splicelogic/reference/find_events.md).
- Facilitating RNA-binding protein (RBP) motif detection.
- Facilitating interpretation of downstream structural consequences.

# Feedback

We would love to hear your feedback. Please post to an [Issue on
GitHub](https://github.com/thelovelab/splicelogic/issues/new).

# Funding

*splicelogic* is supported by NHGRI R01-HG009937, and the Wellcome Trust
as part of the EOSS program.
