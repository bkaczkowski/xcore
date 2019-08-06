### FUNCTION TO PROCESS CTSS FILES TO COUNT TABLE
### USING data.table package for speed

#    1) read_in_ctss : Read in CTSS files to data.table object
#    2) bed_to_anno  : Read in a BED file to crate "anno" data.table object.
#    3) ctss_data_to_counts : Create count table by overlapping CTSS data with genomic coordinates [from 1) and 2) ]
#    4) ctss_to_counts : convinience wrapper for the above functions
# EXAMPLE
         # ctss_files = list_the_ctss_files("~/ctss_files")
         # promoter_counts = ctss_to_counts( ctss_files , "~/promoters.bed" )

#' Read in CTSS files to data.table object
#' @param ctss_files character vector with the paths to CTSS files, names(ctss_files) should contain the sample names
#' @return data.table object with all the data from CTSS files
#' @importFrom data.table fread rbindlist setnames setkey
#' @import R.utils
#' @export
read_in_ctss = function( ctss_files ) {
  #read in the data using fread from data.table package
  ctss_data <- lapply(ctss_files, data.table::fread)

  #insert the name of the file into the 4th column
  for ( i in 1:length(ctss_data)){
    ctss_data[[i]]$V4 = names(ctss_data)[i]
  }
  # rbind the list into one large data.table
  ctss_data <- data.table::rbindlist(ctss_data)

  # prepare the object for further analysis (set names and setkey)
  data.table::setnames(ctss_data, c("chr", "start", "end", "name" , 'score', "strand"))
  data.table::setkey(ctss_data, chr, start, end)
  # return the data.table object ready for the overlap analysis
  ctss_data
}

#' Read in a BED file to crate "anno" data.table object.
#' @param bed_file  paths to bed file (BED4 or BED6)
#' @return data.table object with genomic coordinates (6 columns)
#' @importFrom data.table fread setnames setkey
#' @import R.utils
#' @export
bed_to_anno  = function( bed_file ){
  anno <- data.table::fread(bed_file)
  if( ncol(anno) >= 6){
    anno = anno[ , 1:6 ]
    data.table::setnames(anno, c("chr", "start", "end", "name" , 'score', "strand"))
    data.table::setkey(anno, chr, start, end)
  } else if ( ncol(anno) >= 4 ) {
    anno = anno[ , 1:4 ]
    data.table::setnames(anno, c("chr", "start", "end", "name"))
    anno$score = 0
    anno$strand = "."
    data.table::setkey(anno, chr, start, end)
  } else {
    stop ( "the provided BED file is not compatible with this function")
  }
  anno
}

#' Convert GRanges object to "anno" data.table object
#' @param granges a GRanges object to be converted to data.table
#' @return data.table object with genomic coordinates (6 columns)
#' @importFrom GenomicRanges seqnames start end score strand
#' @export
granges_to_anno = function(granges){

  chr   = as.character(GenomicRanges::seqnames(granges))
  start = GenomicRanges::start(granges) - 1
  end   = GenomicRanges::end(granges)
  name  = granges$name
  score = ifelse ( is.null( GenomicRanges::score(granges))  ,0 , GenomicRanges::score(granges) )
  strand = GenomicRanges::strand(granges)
  strand = sub("\\*", ".", strand)
  anno_df = data.frame( chr = chr, start = start, end = end,name = name, score = score, strand = strand, stringsAsFactors = FALSE )

  anno = data.table::as.data.table(anno_df)
  data.table::setnames(anno, c("chr", "start", "end", "name" , 'score', "strand"))
  data.table::setkey(anno, chr, start, end)
  anno
}


#' Create count table by overlapping CTSS data with genomic coordinates (data.table)
#' @param ctss_data data.table object with CTSS data, created by "read_in_ctss" function
#' @param anno data.table object with genomic coordinates, created by "bed_to_anno" function
#' @return count table with number of reads mapping to genomic regions
#' @importFrom data.table foverlaps
#' @importFrom data.table dcast.data.table
#' @export
ctss_data_to_counts = function( ctss_data , anno) {
  overlaps <- data.table::foverlaps(ctss_data, anno, nomatch=0)
  overlaps <- overlaps[strand==i.strand | strand=='.']
  dt = data.table::dcast.data.table(overlaps, name~i.name, value.var='i.score', fill=0, fun.aggregate=sum)
  count_table = data.matrix(dt [ , -1])
  rownames(count_table) = dt$name
  count_table
}

#' List the CTSS Files
#' @param path path to the directory containing the CTSS files
#' @export
#' @return a character vector of bam files (full path), with "names" representing shortened file names
list_the_ctss_files = function( path ){
  ctss_files = list.files(path = path , pattern = ".ctss.bed.gz$", full.names = T)
  names(ctss_files) = gsub(".*/|.ctss.bed.gz" , "", ctss_files )
  ctss_files
}

#' CTSS files to count table wrapper
#' @param ctss_files a character vector with paths to CTSS files, names(ctss_files) should represent sample names
#' @param regions_bed path to bed file with genomic regions of interest
#' @return (expression) count table of number of reads within provided genomic regions
#' @export
ctss_to_counts = function( ctss_files , regions_bed){
  ctss_data = read_in_ctss(ctss_files)
  anno =  bed_to_anno(regions_bed)
  ctss_data_to_counts ( ctss_data , anno)
}

