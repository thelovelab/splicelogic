# splicelogic: DTU to splice events

*splicelogic* allows users to find alternative splicing events after
performing differential transcript usage (DTU) analysis. Unlike
event-based tools that work at the junction level, *splicelogic*
operates on whole transcript structures: each transcript and all its
exons are annotated with a DTU effect estimate, allowing splicing events
to be derived directly from transcript quantification with full isoform
context. By comparing up- and down-regulated transcripts, *splicelogic*
can detect skipped exons, included exons, mutually exclusive exons,
retained introns, and alternative 5’ and 3’splice sites. Because it
takes transcript-level effect estimates as input, it is compatible with
any upstream DTU method (including DRIMSeq, DEXSeq, satuRn, and edgeR),
supporting flexible experimental designs.

*splicelogic* is included as a R/Bioconductor package for detecting
alternative splicing events from exon-level data stored as *GRanges*
objects. Given a set of exons annotated with a coefficient column
indicating differential transcript usage (DTU), *splicelogic* can be
used to identify a variety of splicing events. See the
[vignette](https://thelovelab.github.io/splicelogic/articles/splicelogic.html)
for more details.

# How to install

`splicelogic` will be submitted to Bioconductor. For now you can test it
by installing from GitHub:

    devtools::install_github("thelovelab/splicelogic")

# Quick start

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

- Support detection of alternative UTR events (alternative 5’ and 3’
  UTRs), when the reference annotation includes UTR coordinates
  (e.g. GENCODE).
- Support detection of additional event types, such as consecutive
  skipped exons or loss of retained introns.
- Extraction and labelling of the specific splice junctions associated
  with each event, adding metadata columns such as the donor–acceptor
  dinucleotide sequence (e.g. AG-GT) and a logical indicating whether
  the junction is canonical, for downstream interpretation.
- Facilitating RNA-binding protein (RBP) motif detection.
- Facilitating interpretation of downstream structural consequences.

# Feedback

We would love to hear your feedback. Please post to an [Issue on
GitHub](https://github.com/thelovelab/splicelogic/issues/new).

# Funding

*splicelogic* is supported by NHGRI R01-HG009937, and the Wellcome Trust
as part of the EOSS program.
