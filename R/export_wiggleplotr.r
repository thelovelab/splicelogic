library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(rlang)
library(purrr)
library(assertthat)

# from wiggleplotr.R
#-------------------------------------------------------------------------------
#' Quickly plot transcript structure without read coverage tracks
#'
#' @param exons list of GRanges objects, each object containing exons for one transcript.
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame.
#' @param cdss list of GRanges objects, each object containing the coding regions (CDS) of a single transcript. 
#' The list must have names that correspond to transcript_id column in transript_annotations data.frame. 
#' If cdss is not specified then exons list will be used for both arguments. (default: NULL)
#' @param transcript_annotations Data frame with at least three columns: transcript_id, gene_name, strand.
#' Used to construct transcript labels. (default: NULL)
#' @param rescale_introns Specifies if the introns should be scaled to fixed length or not. (default: TRUE)
#' @param new_intron_length length (bp) of introns after scaling. (default: 50)
#' @param flanking_length Lengths of the flanking regions upstream and downstream of the gene. (default: c(50,50))
#' @param connect_exons Print lines that connect exons together. Set to FALSE when plotting peaks (default: TRUE).
#' @param transcript_label If TRUE then transcript labels are printed above each transcript. (default: TRUE). 
#' @param region_coords Start and end coordinates of the region to plot, overrides flanking_length parameter.
#'
#' @return ggplot2 object
#' @examples
#' plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE)
#' 
#' @export
myPlotTranscripts <- function(exons, cdss = NULL, transcript_annotations = NULL, 
                            rescale_introns = TRUE, new_intron_length = 50, 
                            flanking_length = c(50,50), connect_exons = TRUE, 
                            transcript_label = TRUE, region_coords = NULL){
  
  #IF cdss is not specified then use exons instead on cdss
  if(is.null(cdss) || length(ccds) == 0){
    cdss = exons
  }
  
  #Check exons and cdss
  assertthat::assert_that(is.list(exons)|| is(exons, "GRangesList")) #Check that exons and cdss objects are lists
  assertthat::assert_that(is.list(cdss) || is(exons, "GRangesList"))
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(joint_exons)[1])
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    min_start = min(GenomicRanges::start(gene_range))
    max_end = max(GenomicRanges::end(gene_range))
    flanking_length = c(min_start - region_coords[1], region_coords[2] - max_end)
  }
  #Make sure that flanking_length is a vector of two elements
  assertthat::assert_that(length(flanking_length) == 2) 
  
  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length, flanking_length)
    xlabel = "Distance from region start (bp)"
  } else {
    old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges),
                          old_introns = old_introns, new_introns = old_introns)
    
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  
  #If transcript annotations are not supplied then construct them manually from the GRanges list
  if(is.null(transcript_annotations)){
    plotting_annotations = dplyr::tibble(transcript_id = names(exons),
                                         strand = extractStrandsFromGrangesList(exons)) %>%
      prepareTranscriptAnnotations()
  } else{
    plotting_annotations = prepareTranscriptAnnotations(transcript_annotations)
  }
  
  #Plot transcript structures
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  structure = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                    tx_annotations$cds_ranges, plotting_annotations)
  plot = plotTranscriptStructure(structure, limits, connect_exons = connect_exons, xlabel = xlabel, 
                                 transcript_label = transcript_label)
  return(plot)
}


# from makePlots.R with beamimc changes
#-------------------------------------------------------------------------------
plotTranscriptStructure <- function(exons_df, limits = NA, connect_exons = TRUE,  
                                    xlabel = "Distance from gene start (bp)", transcript_label = TRUE){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by_(exons_df, ~transcript_id) %>% 
    dplyr::filter_(~feature_type == "exon") %>%
    dplyr::arrange_('transcript_id', 'start') %>%
    dplyr::filter(row_number() == 1)
  
  # Define default fill column
  color_by <- "feature_type"
  colors <- c("cds" = "#2c7bb6", "exon" = "#abd9e9")
  
  # If 'color_by' column exists, overwrite
  if ("color_by" %in% colnames(exons_df)) {
    color_by <- "color_by"
    unique_colors <- unique(exons_df$color_by)
    colors <- setNames(unique_colors, unique_colors)  # name = value
  }
  
  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes(x = start, y = transcript_rank, 
                                group = transcript_rank, 
                                color = !!rlang::sym(color_by)))
  }
  plot = plot + 
    geom_rect(aes(xmin = start, 
                  xmax = end, 
                  ymax = transcript_rank + 0.25, 
                  ymin = transcript_rank - 0.25, 
                  fill = !!rlang::sym(color_by))) + 
    theme_light() +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) +
    xlab(xlabel) +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0.2,0.15)) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors)
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = limits)
  }
  if(transcript_label){
    plot = plot + geom_text(aes_(x = ~start, 
                                 y = ~transcript_rank + 0.30, 
                                 label = ~transcript_label), 
                            data = transcript_annot, hjust = 0, vjust = 0, size = 4)
    
  }
  return(plot)
}



# from functions.R 
#-------------------------------------------------------------------------------
joinExons <- function(exons) {
  #Join a list of exons into one GRanges object
  
  #Test that all transcripts are on the same chromosome
  chrs = purrr::map_chr(as.list(exons), ~GenomicRanges::seqnames(.)[1] %>% 
                          S4Vectors::as.vector.Rle(mode = "character"))
  if (!all(chrs == chrs[1])){
    stop("Some transcripts are on different chromosomes.")
  }
  
  #Join all exons together
  transcript_ids = names(exons)
  joint_exons = c()
  for(tx_id in transcript_ids){
    tx = exons[[tx_id]]
    if(length(joint_exons) == 0){
      joint_exons = tx
    }
    else{
      joint_exons = c(joint_exons, tx)
    }
  }
  joint_exons = GenomicRanges::reduce(joint_exons)
  return(joint_exons)
}

extractStrandsFromGrangesList <- function(granges_list){
  strands = purrr::map(as.list(granges_list), ~(GenomicRanges::strand(.) %>%
                                                  S4Vectors::as.vector.Rle(.,"character"))[1])
  return(unlist(strands))
}

prepareTranscriptAnnotations <- function(transcript_annotations){
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
  
  
  #Make sure that the strand information is represented correctly
  transcript_annotations = dplyr::mutate(transcript_annotations,
                                         strand = ifelse(strand %in% c("+","*") | strand == 1, 1, -1))
  
  #Add transcript label
  if(assertthat::has_name(transcript_annotations, "gene_name")){
    transcript_annotations = dplyr::select_(transcript_annotations, "transcript_id", "gene_name", "strand") %>% 
      dplyr::mutate(transcript_label = ifelse(strand == 1, 
                                              paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                                              paste("< ",paste(gene_name, transcript_id, sep = ":"),sep ="")))
  } else{
    transcript_annotations = dplyr::mutate(transcript_annotations, transcript_label = ifelse(strand == 1, 
                                                                                             paste(paste(transcript_id, sep = ":")," >",sep =""), 
                                                                                             paste("< ",paste(transcript_id, sep = ":"),sep =""))) 
  }
  return(transcript_annotations)
}

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = purrr::map_df(exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = purrr::map_df(cds_ranges, data.frame, .id = "transcript_id")
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Add transcript label to transcript structure
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id")
  return(transcript_struct)
}

intronsFromJointExonRanges <- function(joint_exon_ranges, flanking_length){
  #Construct intron ranges from joint exon ranges
  introns = IRanges::gaps(joint_exon_ranges, 
                          start = min(IRanges::start(joint_exon_ranges)) - flanking_length[1], 
                          end = max(IRanges::end(joint_exon_ranges)) + flanking_length[2])
  return(introns)
}

# Find the start and end coordinates of the whole gene form joint exons. 
constructGeneRange <- function(joint_exon_ranges, flanking_length){
  gene_range = GenomicRanges::reduce(c(joint_exon_ranges, GenomicRanges::gaps(joint_exon_ranges, start = NA, end = NA)))
  GenomeInfoDb::seqlevels(gene_range) = S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")[1]
  GenomicRanges::start(gene_range) = GenomicRanges::start(gene_range) - flanking_length[1]
  GenomicRanges::end(gene_range) = GenomicRanges::end(gene_range) + flanking_length[2]
  return(gene_range)
}

#from shortenIntrons.R
#-------------------------------------------------------------------------------
shortenIntrons <- function(introns, intron_length){
  #Shorten introns from a fixed length to a variable length
  
  #Calculate neccesary parameters
  exons = IRanges::gaps(introns)
  n_introns = length(introns)
  n_exons = length(exons)
  
  #Calculate cumulative with of introns
  intron_cum_width = seq(intron_length,(n_introns-1)*intron_length,intron_length)
  #Calculate new exon starts ignoring introns
  new_intron_starts = c(1,IRanges::start(introns)[2:n_introns] - (IRanges::end(introns)[1:n_introns-1] - intron_cum_width))
  #Add exon widths to the introns
  new_intron_starts = new_intron_starts + c(0,cumsum(IRanges::width(exons)) - IRanges::width(exons))
  
  new_introns = IRanges::IRanges(start = new_intron_starts, width = rep(intron_length, n_introns))
  return(new_introns)
}

translateExonCoordinates <- function(exons, old_introns, new_introns){
  #Tranlate exon coordinates by shortening introns
  old_exon_starts = IRanges::start(exons)
  old_intron_ends = IRanges::end(old_introns)
  new_intron_ends = IRanges::end(new_introns)
  
  #Translate old exon coordinates to new exon coordinates
  new_exon_starts = rep(0,length(old_exon_starts))
  for (i in seq_along(old_exon_starts)){
    #Find the nearest upstream intron for the current gene
    nearest_intron_number = max(which(old_exon_starts[i] > old_intron_ends))
    new_exon_starts[i] = old_exon_starts[i] - old_intron_ends[nearest_intron_number] + new_intron_ends[nearest_intron_number]
  }
  
  #Create new exon coordinates
  new_exons = IRanges::IRanges(start = new_exon_starts, width = IRanges::width(exons))
  return(new_exons)
}

rescaleIntrons <- function(exons, cdss, joint_exons, new_intron_length, flanking_length){
  
  #Convert exons and cds objects to ranges
  exon_ranges = lapply(exons, GenomicRanges::ranges)
  cds_ranges = lapply(cdss, GenomicRanges::ranges)
  
  #Shorten introns and translate exons into the new exons
  old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  new_cds_ranges = lapply(cds_ranges, translateExonCoordinates, old_introns, new_introns)
  
  return(list(exon_ranges = new_exon_ranges, cds_ranges = new_cds_ranges, 
              old_introns = old_introns, new_introns = new_introns))
}
