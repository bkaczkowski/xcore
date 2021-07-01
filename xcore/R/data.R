#' Enahncers GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 enhancers overlapped with EP300 peaks
#' from ReMap2020 and dyadic regions from Dfam database.
#'
#' @format A GenomicRanges object of length 63285, with 4 metadata columns:
#' \describe{
#'   \item{name}{Feature name.}
#'   \item{hg38_enhancer}{FANTOM5's hg38 enhancer name.}
#'   \item{ep300}{Variable indicating overlapp with EP300 ReMap2020 track.}
#'   \item{repeat_dfam}{variable indicating overlapp with Dfam dyadic track.}
#' }
#'
"enhancers"

#' Promoters GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 promoters with FANTOM5's annotations.
#' Additionally this dataset was also annotated to nearest features in GENCODE 
#' ver. 38 annotation and UCSC hg38 knownGene annotation ver. 3.13.0.
#'
#' @format A GenomicRanges object of length 209911, with 11 metadata columns:
#' \describe{
#'   \item{name}{Promotor name.}
#'   \item{score}{Numeric vector.}
#'   \item{SYMBOL_F5_annot}{Gene symbol of associated gene as defined by FANTOM5 annotation.}
#'   \item{ENTREZID_F5_annot}{ENTREZ ID of associated gene as defined by FANTOM5 annotation.}
#'   \item{SYMBOL_gencode}{Gene symbol of associated gene as defined by GENCODE annotation.}
#'   \item{ENTREZID_gencode}{ENTREZ ID of associated gene as defined by GENCODE annotation.}
#'   \item{gene_type_gencode}{Gene type of associated gene as defined by GENCODE annotation.}
#'   \item{SYMBOL_ucsc}{Gene symbol of associated gene as defined by UCSC annotation.}
#'   \item{ENTREZID_ucsc}{ENTREZ ID of associated gene as defined by UCSC annotation.}
#'   \item{ENTREZID}{ENTREZ ID of associated gene as defined by all three annotation sources.
#'                   Prevalence of annotations: UCSC > GENCODE > FANTOM5.}
#'   \item{SYMBOL}{Gene symbol of associated gene as defined by all three annotation sources.
#'                 Prevalence of annotations: UCSC > GENCODE > FANTOM5.}
#' }
#'
"promoters"

#' Detailed promoters GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 promoters with FANTOM5's annotations.
#' Additionally this dataset was also annotated to nearest features in GENCODE 
#' ver. 38 annotation and UCSC hg38 knownGene annotation ver. 3.13.0.
#' Further object was overlapped with EP300 peaks from ReMap2020, \code{\link{enhancers}} object 
#' and dyadic regions from Dfam database. Moreover the dataset is overlapped with GENCODE annotation
#' on both strands.
#'
#' @format A GenomicRanges object of length 209911, with 24 metadata columns:
#' \describe{
#'   \item{name}{Promotor name.}
#'   \item{score}{Numeric vector.}
#'   \item{SYMBOL_F5_annot}{Gene symbol of associated gene as defined by FANTOM5 annotation.}
#'   \item{ENTREZID_F5_annot}{ENTREZ ID of associated gene as defined by FANTOM5 annotation.}
#'   \item{distance_F5_annot}{Distance to associated gene as defined by FANTOM5 annotation.}
#'   \item{SYMBOL_gencode}{Gene symbol of associated gene as defined by GENCODE annotation.}
#'   \item{ENTREZID_gencode}{ENTREZ ID of associated gene as defined by GENCODE annotation.}
#'   \item{gene_type_gencode}{Gene type of associated gene as defined by GENCODE annotation.}
#'   \item{distance_gencode}{Distance to associated gene as defined by GENCODE annotation.}
#'   \item{SYMBOL_ucsc}{Gene symbol of associated gene as defined by UCSC annotation.}
#'   \item{ENTREZID_ucsc}{ENTREZ ID of associated gene as defined by UCSC annotation.}
#'   \item{distance_ucsc}{Distance to associated gene as defined by UCSC annotation.}
#'   \item{ENTREZID}{ENTREZ ID of associated gene as defined by all three annotation sources.
#'                   Prevalence of annotations: UCSC > GENCODE > FANTOM5.}
#'   \item{SYMBOL}{Gene symbol of associated gene as defined by all three annotation sources.
#'                 Prevalence of annotations: UCSC > GENCODE > FANTOM5.}
#'   \item{ep300}{variable indicating overlapp with EP300 ReMap2020 track.}
#'   \item{enhancer}{variable indicating overlapp with \code{\link{enhancers}}.}
#'   \item{repeat_dfam}{variable indicating overlapp with Dfam dyadic track.}
#'   \item{annotation}{Feature annotation type as defined by overlapp with GENCODE annotation.}
#'   \item{symbol}{Gene symbol of associated gene as defined by overlapp with GENCODE annotation.}
#'   \item{gene_type}{Gene type of associated gene as defined by overlapp with GENCODE annotation.}
#'   \item{level}{Character vector.}
#'   \item{opposite_strand_annotation}{Feature annotation type as defined by opposite strand
#'                                     overlapp with GENCODE annotation.}
#'   \item{opposite_strand_symbol}{Gene symbol of associated gene as defined by opposite strand
#'                                 overlapp with GENCODE annotation.}
#'   \item{opposite_gene_type}{Gene type of associated gene as defined by opposite strand overlapp
#'                             with GENCODE annotation.}
#' }
#'
"promoters_detailed"

#' ReMap2020 and FANTOM5 promoters intersection matrix
#'
#' An intersection matrix describing overlaps between ReMap2020's ChIP-seq tracks
#' and \code{\link{promoters}}. To find overlapping regions promoters were extended
#' by 500bp in both directions.
#'
#' @format A Matrix with 209911 rows and 5728 columns. Row names corresponds to promoters
#'         names, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_promoters"

#' ReMap2020 and FANTOM5 enhancers intersection matrix
#'
#' An intersection matrix describing overlaps between ReMap2020's ChIP-seq tracks
#' and \code{\link{enhancers}}. To find overlapping regions enhancers were extended
#' by 500bp in both directions.
#'
#' @format A Matrix with 63285 rows and 5728 columns. Row names corresponds to enhancers
#'         names, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_enhancers"

#' ReMap2020 gene level interaction matrix
#'
#' An matrix describing interactions between ReMap2020's ChIP-seq tracks
#' and human genes (hg38) as defined by \code{\link{promoters}}. The interaction score
#' for each gene and transcription factor is a sum of transcription factor occurences
#' in the gene promoter. Promoters were assigned to thier target genes based on ENTREZ
#' IDs.
#'
#' @format A Matrix with 23274 rows and 5728 columns. Row names corresponds to gene's
#'         ENTREZ IDs, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_entrez"

#' ReMap2020 gene level interaction matrix
#'
#' An matrix describing interactions between ReMap2020's ChIP-seq tracks
#' and human genes (hg38) as defined by \code{\link{promoters}}. The interaction score
#' for each gene and transcription factor is a sum of transcription factor occurences
#' in the gene promoter. Promoters were assigned to thier target genes based on gene symbols.
#'
#' @format A Matrix with 25925 rows and 5728 columns. Row names corresponds to gene's
#'         symbols, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_symbol"

