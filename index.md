# splicelogic: differential transcripts to splice events

# How to install

*splicelogic* will be submitted to Bioconductor. For now you can test it
by installing from GitHub:

    devtools::install_github("thelovelab/splicelogic")

# Quick start

    exons <- prepare_exons(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
      dtu_table = <DTU_TABLE>,
      coef_col = "estimate")

    processed_exons <- preprocess_input(exons, coef_col = "estimate")
    events <- processed_exons |> calc_all_events()

# Feedback

We would love to hear your feedback. Please post to an [Issue on
GitHub](https://github.com/thelovelab/splicelogic/issues/new).

# Funding

splicelogic was supported by NHGRI R01-HG009937, and the Wellcome Trust
as part of the EOSS program.
