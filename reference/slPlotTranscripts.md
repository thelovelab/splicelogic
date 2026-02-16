# Quickly plot transcript structure without read coverage tracks

Quickly plot transcript structure without read coverage tracks

## Usage

``` r
slPlotTranscripts(
  exons,
  cdss = NULL,
  transcript_annotations = NULL,
  rescale_introns = TRUE,
  new_intron_length = 50,
  flanking_length = c(50, 50),
  connect_exons = TRUE,
  transcript_label = TRUE,
  region_coords = NULL
)
```

## Arguments

- exons:

  list of GRanges objects, each object containing exons for one
  transcript. The list must have names that correspond to transcript_id
  column in transript_annotations data.frame.

- cdss:

  list of GRanges objects, each object containing the coding regions
  (CDS) of a single transcript. The list must have names that correspond
  to transcript_id column in transript_annotations data.frame. If cdss
  is not specified then exons list will be used for both arguments.
  (default: NULL)

- transcript_annotations:

  Data frame with at least three columns: transcript_id, gene_name,
  strand. Used to construct transcript labels. (default: NULL)

- rescale_introns:

  Specifies if the introns should be scaled to fixed length or not.
  (default: TRUE)

- new_intron_length:

  length (bp) of introns after scaling. (default: 50)

- flanking_length:

  Lengths of the flanking regions upstream and downstream of the gene.
  (default: c(50,50))

- connect_exons:

  Print lines that connect exons together. Set to FALSE when plotting
  peaks (default: TRUE).

- transcript_label:

  If TRUE then transcript labels are printed above each transcript.
  (default: TRUE).

- region_coords:

  Start and end coordinates of the region to plot, overrides
  flanking_length parameter.

## Value

ggplot2 object

## Examples

``` r
plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE)
#> Error in plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE): could not find function "plotTranscripts"
```
