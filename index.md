# splicelogic: differential transcripts to splice events

# How to install

`splicelogic` will be submitted to Bioconductor. For now you can test it
by installing from GitHub:

    devtools::install_github("thelovelab/splicelogic")

# Quick start

    exons <- prepare_exons(
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
      dtu_table = <DTU_TABLE>,
      coef_col = "estimate"
      )

    processed_exons <- preprocess(exons, coef_col = "estimate")
    events <- processed_exons |> find_all_events()

# Future directions

- Support detection of alternative UTR events (alternative 5’ and 3’
  UTRs), when the reference annotation includes UTR coordinates
  (e.g. GENCODE).
- Support detection of additional event types, such as consecutive
  skipped exons or loss of retained introns.
- Extraction and labelling of the specific splice junctions associated
  with each event, to easily establish whether the junction is canonical
  or not for downstream interpretation.

# Feedback

We would love to hear your feedback. Please post to an [Issue on
GitHub](https://github.com/thelovelab/splicelogic/issues/new).

# Funding

`splicelogic` was supported by NHGRI R01-HG009937, and the Wellcome
Trust as part of the EOSS program.
