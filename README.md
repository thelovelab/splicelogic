# splicelogic: splicing events from transcript sets

_splicelogic_ turns sets of transcripts into discrete splicing events.
Unlike event-based tools that work at the junction level, _splicelogic_ operates on whole transcript structures: within each gene it compares two groups of transcripts and all of their exons, so events are derived with full isoform context. It detects skipped exons, included exons, mutually exclusive exons, retained introns, and alternative 5' and 3' splice sites.

_splicelogic_ does not decide what the two groups are; that is settled upstream, in one of two ways. They can come from an explicit partition of the transcripts into two sets, with no differential analysis behind it — a reference annotation against a set of novel transcripts, say, or a control sample against a disease sample. Or they can come from a differential transcript usage (DTU) analysis, which gives each transcript an effect estimate whose sign defines the comparison; because it takes transcript-level effect estimates as input, _splicelogic_ is compatible with any upstream DTU method (including DRIMSeq, DEXSeq, satuRn, and edgeR), supporting flexible experimental designs.

*splicelogic* operates on exon-level data stored as *GRanges* objects within R/Bioconductor.
Given a set of exons carrying the gene ID, transcript ID and exon rank, _splicelogic_ can be used to identify a variety of splicing events. See the [vignette](https://thelovelab.github.io/splicelogic/articles/splicelogic.html) for more details.

# How to install

`splicelogic` will be submitted to Bioconductor. For now you can test it by
installing from GitHub:

```
devtools::install_github("thelovelab/splicelogic")
```

# Quick start

## Finding events from transcript sets

```
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
```

## Finding events from DTU results

```
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
```

# Future directions

- Support detection of alternative UTR events (alternative 5' and 3' UTRs), when the reference annotation includes UTR coordinates (e.g. GENCODE).
- Support detection of additional event types, such as consecutive skipped exons or loss of retained introns.
- Extraction and labelling of the specific splice junctions associated with each event, adding metadata columns such as the donor–acceptor dinucleotide sequence (e.g. AG-GT) and a logical indicating whether the junction is canonical, for downstream interpretation.
- Facilitating RNA-binding protein (RBP) motif detection.
- Facilitating interpretation of downstream structural consequences.

# Feedback

We would love to hear your feedback. Please post to an 
[Issue on GitHub](https://github.com/thelovelab/splicelogic/issues/new).

# Funding

_splicelogic_ is supported by NHGRI R01-HG009937, and
the Wellcome Trust as part of the EOSS program.
