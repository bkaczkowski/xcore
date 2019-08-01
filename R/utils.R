# FUNCTIONS:
#    1) convert_ctss_to_bigwig : Convers CTSS bed file to  BigWig


#' Convers CTSS bed file to  BigWig
#'
#' @param ctss_file path to ctss.bed file
#' @param genomeInfo genome info file created e.g. by "rtracklayer::SeqinfoForUCSCGenome("hg38")"
#' @return The function writes two files, one for each strand, (".plus.bw" and ".minus.bw"), in the same directory as the original file
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import export
#' @export
convert_ctss_to_bigwig = function( ctss_file, genomeInfo ) {
  tpm_bed = rtracklayer::import(ctss_file)
  tpm_bed = GenomicRanges::GRanges(tpm_bed , seqinfo =  genomeInfo)
  tpm_bed_plus  = tpm_bed[ tpm_bed@strand == "+", ]
  tpm_bed_minus = tpm_bed[ tpm_bed@strand == "-", ]
  rtracklayer::export(object = tpm_bed_plus , paste(ctss_file ,   ".plus.bw" , sep = "")  , format = "BigWig" )
  rtracklayer::export(object = tpm_bed_minus, paste(ctss_file ,  ".minus.bw" , sep = "")  , format = "BigWig" )
}

#' Normalize the count table with edgeR
#'
#' @param counts expression table with counts
#' @param method passed to to edgeR::calcNormFactors see ?edgeR::calcNormFactors
#' @param log passed on to edgeR::cpm , should the normalized counts be log2 transformed (default FALSE)
#' @param prior.count average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
#' @param ... other parameter passed on to edgeR::calcNormFactors or edgeR::cpm
#' @return normalized count table
#' @import edgeR
#' @export
normalize_counts = function( counts, method = c("TMM","TMMwsp","RLE","upperquartile","none") ,
                             log = FALSE , prior.count = 1 , ...){
  dge = DGEList(counts)
  dge = calcNormFactors(dge , method= method , ... )
  cpm <- cpm( dge , log = log, prior.count = prior.count, ... )
  cpm
}
