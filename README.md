# splicelogic: DTU to splice events

_splicelogic_ is an R/Bioconductor package for detecting alternative splicing events from exon-level data stored as _GRanges_ objects. Given a set of exons annotated with a coefficient column indicating differential transcript usage (DTU), _splicelogic_ can be used to identify a variety of splicing events. See the [vignette](articles/splicelogic.html) for more details.

# How to install

`splicelogic` will be submitted to Bioconductor. For now you can test it by
installing from GitHub:

```
devtools::install_github("thelovelab/splicelogic")
```

# Quick start

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
