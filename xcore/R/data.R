#' Enahncers GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 enhancers overlapped with EP300 peaks
#' from ReMap2020 and dyadic regions from Dfam database.
#'
#' @format A GenomicRanges object of length 63285, with 4 metadata columns:
#' \describe{
#'   \item{name}{feature name}
#'   \item{hg38_enhancer}{FANTOM5's hg38 enhancer name}
#'   \item{ep300}{variable indicating overlapp with EP300 ReMap2020 track}
#'   \item{repeat_dfam}{variable indicating overlapp with Dfam dyadic track}
#' }
#'
"enhancers"
