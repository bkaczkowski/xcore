### !!! THIS MODULE SHOULD WORK BUT NOT THOROUGHLY TESTED, TESTING IS REQUIRED !!!

### FUNCTION TO PROCESS BAM FILES TO COUNT TABLE
### USING GenomicAlignments and GenomicRanges packages from Bioconductor


#    1) bam_to_counts_5prime : Create count table by overlapping BAM data (5'starts) with genomic coordinates 
#    2) bam_to_counts : convinience wrapper for the above function, EXAMPLE:
                              # bam_files = list_the_bam_files("~/bam_files")
                              # promoter_counts = bam_files( ctss_files , "~/promoters.bed" )

#' Bam to Counts function (5'starts)
#' A function that counts a number of reads mapping to a set of genomic regions. 
#' 
#' @param bam_file The path to the bam file with mapped reads.
#' @param regions  A GRanges object with genomic regions e.g. promoter regions
#' @param mapq a number specifying the minimum MAPping Quality (MAPQ) score (reads with lower mapq will not be use)
#' @return a numeric vector of counts
#' @export
#  @import IRanges this MAY be needed, for testing
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicRanges granges resize countOverlaps
bam_to_counts_5prime = function(bam_file, regions ,mapqFilter = 20){
  aln <- GenomicAlignments::readGAlignments(file=bam_file,
                                            param = ScanBamParam(mapqFilter=mapqFilter))
  aln <- GenomicRanges::granges(aln)
  
  gr5prime = GenomicRanges::resize(aln, width = 1, fix="start")
  hits     = GenomicRanges::countOverlaps(gr5prime , regions)
  counts   = GenomicRanges::countOverlaps(regions, aln[hits==1]) # THIS may have to be adjusted to allow for multiple mapping
  names(counts) = regions$id
  data.matrix(counts)
}

#' List the BAM Files
#' @param path path to the directory containing the BAM files
#' @export
#' @return a character vector of bam files (full path), with "names" representing shortened file names
list_the_bam_files = function( path ){
  bam_files = list.files(path = path, pattern = ".bam", full.names = T)
  names(bam_files) = gsub(".*/|\\.bam", "", bam_files)
  bam_files
}

#' Bam to Counts (5'starts) wrapper
#' A convenience wrapper to apply the bam_to_counts_5prime to list (character vector) of bam files 
#' @param bam_files a character vector with paths to BAM files, names(bam_files) should represent sample names
#' @param regions_bed path to bed file with genomic regions of interest
#' @param mapqFilter a optional MAPQ, read quality threslod
#' @return (expression) count table of number of reads within provided genomic regions
#' @importFrom rtracklayer import.bed
#' @export
bam_to_counts = function( bam_files , regions_bed, mapqFilter = 20){
  regions = rtracklayer::import.bed(regions_bed)
  counts = sapply(bam_files , bam_to_counts_5prime, regions = regions, mapqFilter = mapqFilter )
  colnames(counts) = names(bam_files )
  counts
}




