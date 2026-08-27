# to make the mock_exons.bed.gz and mock_exons_mcols.tsv.gz

# This is a small, fully synthetic exon set used to demonstrate importing
# exons from BED and running the splicelogic pipeline, without needing a
# TxDb or a real DTU results table.

# All transcripts belong to one gene, gene_1, and every one of them is
# derived from the same 5-exon baseline, tx_0. tx_0 is the anchor, the
# transcript to pass as `down`: the event finders look for neg exons that
# are missing from, or truncated in, the pos transcripts, so putting the
# baseline in the neg group is what lets each variant be read against it.
 

library(plyranges)
library(dplyr)

# the anchor: 5 evenly spaced exons that every transcript is derived from
baseline <- data.frame(
  exon = 1:5,
  start = c(1001, 1201, 1401, 1601, 1801),
  end = c(1100, 1300, 1500, 1700, 1900)
)

# derive a transcript from the baseline
#   skip: baseline exon numbers to drop entirely
#   alt3: named numeric, names are baseline exon numbers and values the new
#         start — moving the acceptor while keeping the donor is an
#         alternative 3' splice site on the + strand
#   alt5: named numeric, same but values are the new end — moving the donor
#         while keeping the acceptor is an alternative 5' splice site
# exon_rank is renumbered after dropping, so ranks stay consecutive from 1
# exon_id is built from the gene and the final coordinates
make_tx <- function(
  tx_id,
  skip = integer(0),
  alt3 = numeric(0),
  alt5 = numeric(0)
) {
  baseline |>
    dplyr::filter(!exon %in% skip) |>
    dplyr::mutate(
      start = dplyr::coalesce(unname(alt3[as.character(exon)]), start),
      end = dplyr::coalesce(unname(alt5[as.character(exon)]), end),
      seqnames = "chr1",
      strand = "+",
      gene_id = "gene_1",
      tx_id = tx_id,
      exon_id = paste0(gene_id, ":", start, "-", end),
      exon_rank = dplyr::row_number()
    ) |>
    dplyr::select(-exon)
}

exons <- dplyr::bind_rows(
  # anchor: the unmodified baseline, the transcript to pass as `down`
  make_tx("tx_0"),
  # skipped exon
  make_tx("tx_1", skip = 3),
  # alternative 3' splice site
  make_tx("tx_4", alt3 = c("4" = 1651)),
  # alternative 5' splice site on exon 3
  make_tx("tx_2", alt5 = c("3" = 1450)),
  # both at once
  make_tx("tx_5", skip = 3, alt3 = c("4" = 1651))
) |>
  plyranges::as_granges()


# write the ranges to BED and the metadata columns alongside it. The two
# files are matched back up by row order on read, so the BED needs no name
# column of its own.
dir <- "inst/extdata"
plyranges::write_bed(exons, file.path(dir, "mock_exons.bed.gz"))
readr::write_tsv(
  as.data.frame(GenomicRanges::mcols(exons)),
  file.path(dir, "mock_exons_mcols.tsv.gz")
)
