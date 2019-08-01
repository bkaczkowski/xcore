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
