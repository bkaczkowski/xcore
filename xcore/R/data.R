#' Enahncers GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 enhancers overlapped with EP300 peaks
#' from ReMap2020 and dyadic regions from Dfam database.
#'
#' @format A \code{GenomicRanges} object of length 63285, with 4 metadata columns:
#' \describe{
#'   \item{name}{Feature name.}
#'   \item{hg38_enhancer}{FANTOM5's hg38 enhancer name.}
#'   \item{ep300}{Variable indicating overlapp with EP300 ReMap2020 track.}
#'   \item{repeat_dfam}{Variable indicating overlap with Dfam dyadic track.}
#'   \item{tau}{Numeric Tau tissue specificity score calculated based on FANTOM5
#'              expression data.}
#'   \item{encode_blacklist}{Logical indicating overlap with ENCODE black
#'                           listed regions.}
#'   \item{ensembl_sm}{Logical indicating overlap with ENSEMBL soft masked
#'                     regions.}
#' }
#'
"enhancers_f5"

#' Promoters GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 promoters with FANTOM5's annotations.
#' Additionally this dataset was also annotated to nearest features in GENCODE
#' ver. 38 annotation and UCSC hg38 knownGene annotation ver. 3.13.0.
#'
#' @format A \code{GenomicRanges} object of length 209911, with 11 metadata columns:
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
#'   \item{tau}{Numeric Tau tissue specificity score calculated based on FANTOM5
#'              expression data.}
#'   \item{protein_atlas_tissue_specificity}{Character giving tissue specificity
#'                                           class as defined in Protein Atlas.}
#'   \item{protein_atlas_cell_specificity}{Character giving cell type specificity
#'                                           class as defined in Protein Atlas.}
#'   \item{encode_blacklist}{Logical indicating overlap with ENCODE black
#'                           listed regions.}
#'   \item{ensembl_sm}{Logical indicating overlap with ENSEMBL soft masked
#'                     regions.}
#' }
#'
"promoters_f5"

#' Detailed promoters GenomicRanges object
#'
#' A dataset containing FANTOM5's hg38 promoters with FANTOM5's annotations.
#' Additionally this dataset was also annotated to nearest features in GENCODE
#' ver. 38 annotation and UCSC hg38 knownGene annotation ver. 3.13.0.
#' Further object was overlapped with EP300 peaks from ReMap2020, \code{\link{enhancers_f5}} object
#' and dyadic regions from Dfam database. Moreover the dataset is overlapped with GENCODE annotation
#' on both strands.
#'
#' @format A \code{GenomicRanges} object of length 209911, with 24 metadata columns:
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
#'   \item{enhancer}{variable indicating overlapp with \code{\link{enhancers_f5}}.}
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
#' and \code{\link{promoters_f5}}. To find overlapping regions promoters were extended
#' by 500bp in both directions.
#'
#' @format A Matrix with 209911 rows and 5728 columns. Row names corresponds to promoters
#'         names, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_promoters"

#' ReMap2020 and FANTOM5 enhancers intersection matrix
#'
#' An intersection matrix describing overlaps between ReMap2020's ChIP-seq tracks
#' and \code{\link{enhancers_f5}}. To find overlapping regions enhancers were extended
#' by 500bp in both directions.
#'
#' @format A Matrix with 63285 rows and 5728 columns. Row names corresponds to enhancers
#'         names, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_enhancers"

#' ReMap2020 gene level interaction matrix
#'
#' An matrix describing interactions between ReMap2020's ChIP-seq tracks
#' and human genes (hg38) as defined by \code{\link{promoters_f5}}. The interaction score
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
#' and human genes (hg38) as defined by \code{\link{promoters_f5}}. The interaction score
#' for each gene and transcription factor is a sum of transcription factor occurences
#' in the gene promoter. Promoters were assigned to thier target genes based on gene symbols.
#'
#' @format A Matrix with 25925 rows and 5728 columns. Row names corresponds to gene's
#'         symbols, column names are formatted as ExperimentID.TranscriptionFactor.Biotype.
"remap_symbol"

#' ChIP-Atlas FANTOM5 promoters intersection matrix
#'
#' An intersection matrix describing overlaps between ChIP-Atlas's ChIP-seq tracks
#' and \code{\link{promoters_f5}}. To find overlapping regions promoters were extended
#' by 500bp in both directions.
#'
#' @format A \code{Matrix} with 209911 rows and 15133 columns. Row names corresponds to promoters
#'         names, column names are formatted as TranscriptionFactor_Origin_Cell_ExperimentID
#'	   (eg. PARK7_Neural_SH-SY5Y_DRX000550, MLL-AF6_Blood_ML-2_DRX001460).
#'
"chip_atlas_promoters"

#' ChIP-Atlas FANTOM5 enhancers intersection matrix
#'
#' An intersection matrix describing overlaps between ChIP-Atlas's ChIP-seq tracks
#' and \code{\link{enhancers_f5}}. To find overlapping regions promoters were extended
#' by 500bp in both directions.
#'
#' @format A \code{Matrix} with 63285 rows and 15133 columns. Row names corresponds to promoters
#'         names, column names are formatted as TranscriptionFactor_Origin_Cell_ExperimentID
#'	   (eg. PARK7_Neural_SH-SY5Y_DRX000550, MLL-AF6_Blood_ML-2_DRX001460).
#'
"chip_atlas_enhancers"

#' ChIP-Atlas gene level interaction matrix
#'
#' A matrix describing interactions between ChIP-Atlas's ChIP-seq tracks
#' and human genes (hg38) as defined by \code{\link{promoters_f5}}. The interaction score
#' for each gene and transcription factor is a *maximum* of transcription factor occurences
#' among the promoters associated with a given gene. Promoters were assigned to thier target
#' genes based on ENTREZ IDs.
#'
#' @format A \code{Matrix} with 23274 rows and 15133 columns. Row names corresponds to gene's
#'         ENTREZ IDs, column names are formatted as TranscriptionFactor_Origin_Cell_ExperimentID
#'	   (eg. PARK7_Neural_SH-SY5Y_DRX000550, MLL-AF6_Blood_ML-2_DRX001460).
#'
"chip_atlas_entrez"

#' ChIP-Atlas gene level interaction matrix
#'
#' A matrix describing interactions between ChIP-Atlas's ChIP-seq tracks
#' and human genes (hg38) as defined by \code{\link{promoters_f5}}. The interaction score
#' for each gene and transcription factor is a *maximum* of transcription factor occurences
#' among the promoters associated with a given gene. Promoters were assigned to thier target
#' genes based on gene's symbols.
#'
#' @format A \code{Matrix} with 25925 rows and 15133 columns. Row names corresponds to gene's
#'         symbols, column names are formatted as TranscriptionFactor_Origin_Cell_ExperimentID
#'	   (eg. PARK7_Neural_SH-SY5Y_DRX000550, MLL-AF6_Blood_ML-2_DRX001460).
#'
"chip_atlas_symbol"

#' ReMap2020 metadata
#'
#' Metadata associated with \code{remap_promoters} and \code{remap_enhancers}.
#'
#' @format A \code{data.table} with 5798 rows and 5 columns.
#' \describe{
#'   \item{id}{Character giving internal experiment ID.}
#'   \item{tf}{Character giving transcription factor name.}
#'   \item{biotype}{Character giving experiment biological origin.}
#'   \item{study}{Character giving study ID.}
#'   \item{condition}{Character specifiying experiment conditions or treatment.}
#' }
#'
"remap_meta"

#' ChIP-Atlas metadata
#'
#' Metadata associated with \code{chip_atlas_promoters} and
#' \code{chip_atlas_enhancers}.
#'
#' @format A \code{data.frame} with 13923 rows and 5 columns.
#' \describe{
#'   \item{id}{Character giving internal experiment ID.}
#'   \item{tf}{Character giving transcription factor name.}
#'   \item{biotype}{Character giving experiment biological origin.}
#'   \item{study}{Character giving study ID.}
#'   \item{background}{Character specifying biological background as defined in ChIP-Atlas.}
#' }
#'
"chip_atlas_meta"
