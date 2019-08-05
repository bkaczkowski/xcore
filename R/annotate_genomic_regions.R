#' Annotation Worker by Same Strand Overlap
#'
#' Usually used by the gencode_one_direction_annotator function.
#' @param regions query genomic regions as GRanges object
#' @param gencode Gencode annotation as GRanges object, usually of one type like gencode$type == "exon"
#' @param annotatation_label label (character) to mark what annotation was used, e.g. "gene", "exon" or "promoters"
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges findOverlaps
#' @export
gencode_one_direction_findOverlaps = function( regions , gencode, annotatation_label = " "){
  # ordering so higher level gene annotation will overwrie lower level annotation
  # if multiple genes are annotated to same region
  gencode = gencode[ order( gencode$level , decreasing = F),]

  if( is.null(regions$annotatation) ) { regions$annotatation = ""}
  if( is.null(regions$symbol) )       { regions$symbol       = ""}
  if( is.null(regions$gene_type) )    { regions$gene_type    = ""}
  if( is.null(regions$level) )        { regions$level        = ""}

  hits = GenomicRanges::findOverlaps( regions, gencode , ignore.strand=FALSE)
  regions$symbol[hits@from]       = gencode$gene_name[hits@to]
  regions$gene_type[hits@from]    = gencode$gene_type[hits@to]
  regions$level[hits@from]        = gencode$level[hits@to]
  regions$annotatation[hits@from] = annotatation_label

  regions
}

#' Annotating Same Strand by Overlap
#' A wrapper that calls gencode_one_direction_findOverlaps on all types of Gencode Annotation
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @export
gencode_one_direction_annotator = function( regions , gencode ){

  gencode_gene        = gencode[ gencode$type == "gene"]
  gencode_exon        = gencode[ gencode$type == "exon"]
  gencode_3_prime_UTR = gencode[ gencode$type == "three_prime_UTR"]
  gencode_5_prime_UTR = gencode[ gencode$type == "five_prime_UTR"]
  gencode_promoters   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 500, downstream = 500, use.names=TRUE)

  regions = gencode_one_direction_findOverlaps( regions , gencode_gene        , annotatation_label = "intron" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_exon        , annotatation_label = "exon" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_3_prime_UTR , annotatation_label = "three_prime_UTR" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_5_prime_UTR , annotatation_label = "five_prime_UTR" )
  regions = gencode_one_direction_findOverlaps( regions , gencode_promoters   , annotatation_label = "promoter" )

  regions
}

#' Annotating Both Strands by Overlap
#'
#' A wrapper that calls gencode_one_direction_annotator on same and inverted strand
#' to annotate overlap on both strands
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @export
gencode_two_direction_annotator = function(regions , gencode  ){

  regions_inverted_strand = regions
  strand(regions_inverted_strand) = invertStrand(strand(regions))

  regions  = gencode_one_direction_annotator( regions = regions, gencode = gencode)
  regions_inverted_strand  = gencode_one_direction_annotator( regions = regions_inverted_strand, gencode = gencode)

  regions$opposite_strand_annotatation = regions_inverted_strand$annotatation
  regions$opposite_strand_symbol = regions_inverted_strand$symbol
  regions$opposite_gene_type = regions_inverted_strand$gene_type

  regions
}

#' Annotating Strandless regions by Overlap
#'
#' This is mostly for bidirectional/strandless enhancer regions
#' The wrappers call the gencode_one_direction_annotator twice after setting the strand to + and -.
#' Finally the output is nicely summarized
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges strand mcols
#' @export
gencode_strandless_annotator = function( regions, gencode , output_option = 2) {

  regions_tmp = regions
  strand(regions_tmp) = "+"
  mcols(regions_tmp)= NULL
  regions_annotated_plus  = gencode_one_direction_annotator( regions = regions_tmp, gencode = gencode)
  regions_annotated_plus$level = NULL
  regions_annotated_plus  = mcols(regions_annotated_plus)

  strand(regions_tmp) = "-"
  regions_annotated_minus = gencode_one_direction_annotator( regions = regions_tmp, gencode = gencode)
  regions_annotated_minus$level = NULL
  regions_annotated_minus = mcols(regions_annotated_minus)

  if( output_option == 1 ){
    regions$strandless_annotatation = paste( regions_annotated_plus$annotatation, "(+);", regions_annotated_minus$annotatation , "(-)", sep = "" )
    regions$strandless_symbol = paste( regions_annotated_plus$symbol, "(+);", regions_annotated_minus$symbol , "(-)", sep = "" )
    regions$strandless_gene_type = paste( regions_annotated_plus$gene_type, "(+);", regions_annotated_minus$gene_type , "(-)", sep = "" )

  }else if( output_option == 2) {
    colnames(regions_annotated_minus) = paste( colnames(regions_annotated_minus), "_minus",sep = "")
    colnames(regions_annotated_plus)  = paste( colnames(regions_annotated_plus ), "_plus" ,sep = "")
    mcols(regions) = DataFrame( cbind(mcols(regions) ,  regions_annotated_plus, regions_annotated_minus))
  }

  regions
}

#' Annotation Worker by Same Strand DistanceToNearest
#'
#' @param regions query genomic regions as GRanges object
#' @param gencode Gencode annotation as GRanges object, usually of one type like gencode[ gencode$type == "exon"]
#' @param annotatation_label label (character) to mark what annotation was used, e.g. "gene", "exon" or "promoters"
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges distanceToNearest
#' @export
gencode_one_direction_distanceToNearest = function( regions , gencode){

  # ordering so higher level gene annotation will overwrie lower level annotation
  # if multiple genes are annotated to same region
  gencode = gencode[ order( gencode$level , decreasing = F),]

  if( is.null(regions$nearest_symbol) )       { regions$nearest_symbol         = ""}
  if( is.null(regions$nearest_gene_type) )    { regions$nearest_gene_type      = ""}
  if( is.null(regions$nearest_distance) )     { regions$nearest_distance       = -1}

  hits = GenomicRanges::distanceToNearest( regions, gencode , ignore.strand=FALSE)
  regions$nearest_symbol[hits@from]       = gencode$gene_name[hits@to]
  regions$nearest_gene_type[hits@from]    = gencode$gene_type[hits@to]
  regions$nearest_distance                = hits@elementMetadata$distance

  regions
}

#' Annotate with Nearest Promoter for Both Strands
#'
#' The wrappers call the gencode_one_direction_distanceToNearest twice after setting the strand to + and -.
#' @param regions query genomic regions as GRanges object
#' @param gencode complete Gencode annotation as GRanges object, e.g. obrained by rtracklayer::import.gff(con = "gencode.v31.annotation.gff3.gz") ,  (gencode[ gencode$type == "transcript"] is also OK)
#' @return GRanges with new columns added to the feature metadata: mcols(regions)
#' @importFrom GenomicRanges promoters mcols
#' @importFrom S4Vectors DataFrame
#' @export
gencode_nearest_promoter_both_strands = function( regions, gencode ) {

  gencode_promoters   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 0, downstream = 0, use.names=TRUE)

  regions_tmp = regions
  mcols(regions_tmp)= NULL

  strand(regions_tmp) = "+"
  nearest_promoter_plus  = gencode_one_direction_distanceToNearest( regions = regions_tmp, gencode = gencode_promoters)
  nearest_promoter_plus  = mcols(nearest_promoter_plus)
  colnames(nearest_promoter_plus) = paste( colnames(nearest_promoter_plus), "_plus",sep = "")

  strand(regions_tmp) = "-"
  nearest_promoter_minus = gencode_one_direction_distanceToNearest( regions = regions_tmp, gencode = gencode_promoters)
  nearest_promoter_minus = mcols(nearest_promoter_minus)
  colnames(nearest_promoter_minus) = paste( colnames(nearest_promoter_minus), "_minus",sep = "")

  mcols(regions) = DataFrame( cbind( mcols(regions) , nearest_promoter_plus, nearest_promoter_minus))

  regions
}
