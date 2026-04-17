# to make the dtu_table.tsv

# first we downloaded a DTU results table from this zenodo URL
# https://zenodo.org/records/10381745

# the zenodo entry for the data from this paper lists an MIT license

# after unzipping `mouse_brain_iso_div_data.tar.gz`
# we have a number of directories including `/switchlist_objects`

# from authors on `/switchlist_objects/raw`:
# "This directory contains the initial isoformSwitchAnalyzeR objects, without additional information added." 

library(dplyr)
library(tibble)
x <- readRDS("switchlist_objects/raw/region_sex_switchlist_list.Rds")
dat <- x$cortex$isoformFeatures |> as_tibble()
dat |> filter(gene_switch_q_value < .05)
# this is a 49 x 29 tibble of the transcripts in genes with DTU
# this was then written to the TSV file `dtu_table.tsv`

# to make the exons_M31.bed and exons_mcols.tsv.gz
# we downloaded GENCODE M31 GTF file
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz
# then made into a TxDb using txdbmaker

# get exons grouped by transcript
ebt <- exonsBy(txdb, by="tx")
txp <- transcripts(txdb)
names(ebt) <- txp$tx_name
exons <- unlist(ebt)
exons$tx_id <- names(exons)

# this was then subset to the transcripts in `dtu_table.tsv`
# the ranges were written using rtracklayer::export(), and
# the metadata columns written out using write.delim(mcols(exons), ...)
