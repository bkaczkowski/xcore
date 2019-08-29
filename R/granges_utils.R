
#' Name the regions based on genomic coordinates
#' @param granges a GRanges object
#' @return a character vector of coordinate-based names
#' @export
name_from_granges = function( granges ){
  name =  paste( seqnames(granges), ":",
                 start(granges)-1, "..",
                 end(granges), "," ,
                 strand( granges ),
                 sep = "")
  name
}

#' Merge list of GRanges
#' @param list_of_granges a list, each element should be a GRanges object
#' @return a GRanges object with all the regions merged
#' @importFrom S4Vectors append
#' @importFrom GenomicRanges GRanges
#' @export
merge_list_of_granges = function( list_of_granges ){
  granges_merged <- GenomicRanges::GRanges()
  for (gr in list_of_granges ) {
    granges_merged <- S4Vectors::append(granges_merged, gr)
  }
  granges_merged
}
#' Create Intersection Matrix
#' this functions creates a "set intersections in a matrix layout", the format used by UpSetR package to visualize intersections/overlaps bewteen different sets.
#' @param list_of_granges a list of GRanges containing regions of intereset, names(x) will be used to name the columns of the output matrix
#' @param pre_defined_regions a list of predefined regions to be used at rows, if NULL (default), the union of all genomic regions in the list_of_granges will be used.
#' @return intersection matrix
#' @export
#' @importFrom GenomicRanges reduce findOverlaps
#'
create_intersection_matrix   = function( list_of_granges , pre_defined_regions = NULL ) {

  if ( is.null( pre_defined_regions) ) {
    granges_merged   = merge_list_of_granges(list_of_granges)
    granges_reduced  = GenomicRanges::reduce( granges_merged)
    granges_reduced$name  = name_from_granges ( granges_reduced )
    regions = granges_reduced
  } else {
    regions = pre_defined_regions
  }

  intersection_matrix = matrix (0L, nrow = length(regions), ncol = length( list_of_granges ) )
  colnames(intersection_matrix) = names(list_of_granges)
  rownames(intersection_matrix) = regions$name

  for (i in 1:length( list_of_granges ) ) {
    peaks = list_of_granges[[i]]
    hits = GenomicRanges::findOverlaps(regions ,peaks)
    intersection_matrix[ hits@from  , i] =  1L
  }
  data.frame(intersection_matrix)
}


#' Create Intersection Vector
#'
#' @param query GRanges containing query regions e.g. Transcription Factor Chip-Seq data
#' @param pre_defined_regions a set set of reference regions like known promoters or enhancers
#'
#' @return a vector of the same lenght as pre_defined_regions, 0 if no overlap, 1 if there is overlap with query regions
#' @export

create_intersection_vector   = function( query , pre_defined_regions ) {

  intersection_vector = rep (0L,times =  length( pre_defined_regions ) )

  hits = GenomicRanges::findOverlaps(pre_defined_regions ,query)
  intersection_vector[ hits@from ] =  1L

  intersection_vector
}
