#' Process count matrix for modeling # TODO
#' 
#' @param counts matrix of read counts.
#' @param groups factor giving the experimental group for each sample.
#' @param base_lvl string indicating level in \code{groups} corresponding to
#'   basal expression level. The reference samples to which expression change
#'   will be compared.
#' @param log2 logical flag indicating if counts should be log2(counts per 
#'   million) should be returned.
#' @param pseudo_count integer count to be added before taking log2.
#' @param drop_base_lvl logical flag indicating if \code{base_lvl} samples
#'   should be dropped from resulting MultiAssayExperiment object.
#'   
#' @return MultiAssayExperiment object with two experiments:
#'   \describe{
#'     \item{U}{matrix giving expression values averaged over basal level samples}
#'     \item{Y}{matrix of expression values}
#'   }
#'
#' @importFrom edgeR DGEList calcNormFactors cpm filterByExpr
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#'
#' @export
prepareCountsForRegression <- function(counts,
                                       base_lvl,
                                       groups, # TODO change to design, should be stored in MAE for further use
                                       log2 = TRUE,
                                       pseudo_count = 1,
                                       drop_base_lvl = TRUE) {
  deglist <- edgeR::DGEList(counts = mat, group = groups)
  keep <- edgeR::filterByExpr(deglist)
  deglist <- deglist[keep, , keep.lib.sizes=FALSE] # if not keep some genes have 0 sd
  deglist <- edgeR::calcNormFactors(deglist)
  cpm <- edgeR::cpm(deglist)
  
  if (log2) {
    cpm <- log2(cpm + pseudo_count)
  }
  # TODO mean log2 is disputable
  U <- matrix(
    data = rowMeans(cpm[, groups == base_lvl, drop = FALSE]),
    ncol = 1L,
    dimnames = list(rownames(cpm), "u"))
  
  ymask <- if (drop_base_lvl) { groups != base_lvl } else { rep(TRUE, length(groups)) }
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(
      U = U, 
      Y = cpm[, ymask, drop = FALSE]))
}

#' Add signatures matrices to MAE
#'
#'
addSignatures <- function(mae, ..., intersect_rows = TRUE) {
  ex <- list(...)
  ex_names <- names(ex)
  if (is.null(ex_names) || any(ex_names == "")) {
    stop("All experiments must be named")
  }
  
  mae_map <- MultiAssayExperiment::sampleMap(mae)
  ex_map <-
    lapply(ex, function(x)
      S4Vectors::DataFrame(primary = colnames(x), colname = colnames(x)))
  names(ex_map) <- ex_names
  newListMap <- c(MultiAssayExperiment::mapToList(mae_map),
                  IRanges::SplitDataFrameList(ex_map))
  newSampleMap <- MultiAssayExperiment::listToMap(newListMap)
  
  newExperimentList <- c(MultiAssayExperiment::experiments(mae), ex)
  
  # TODO preserve mae colData
  newColData <- S4Vectors::DataFrame(row.names = newSampleMap[["primary"]])
  
  new_mae <- BiocGenerics:::replaceSlots(object = mae,
                                         ExperimentList = newExperimentList,
                                         sampleMap = newSampleMap,
                                         colData = newColData)
  
  if (intersect_rows) new_mae <- MultiAssayExperiment::intersectRows(new_mae)
  
  return(new_mae)
}

#' Compute interaction matrix
#'
#' \code{getInteractionMatrix} construct interaction matrix between two Granges
#' objects. Names of object \code{a} became row names and names of \code{b}
#' column names.
#'
#' @param a GRanges object.
#' @param b GRanges object.
#' @param ext Integer specifying number of base pairs the \code{a} coordinates
#'   should be extended in upstream and downstream directions.
#' @param count Logical indicating if matrix should hold number of overlaps
#'   between \code{a} and \code{b} or if FALSE presence / absence indicators.
#'
#' @return Sparse matrix of class dgCMatrix, with rows corresponding to
#'   \code{a} and columns to \code{b}. Each cell holds a number indicating
#'   how many times \code{a} and \code{b} overlapped.
#'
#' @importFrom GenomicRanges start end findOverlaps
#' @importFrom Matrix sparseMatrix
#' 
#' @examples
#' a <- GenomicRanges::GRanges(
#'   seqnames = c("chr20", "chr4"),
#'   ranges = IRanges::IRanges(
#'     start = c(62475984L, 173530220L), 
#'     end = c(62476001L, 173530236L)), 
#'   strand = c("-", "-"),
#'   name = c("hg19::chr20:61051039..61051057,-;hg_188273.1", 
#'            "hg19::chr4:174451370..174451387,-;hg_54881.1"))  
#' b <- GenomicRanges::GRanges( 
#'   seqnames = c("chr4", "chr20"),
#'   ranges = IRanges::IRanges(
#'     start = c(173530229L, 63864270L), 
#'     end = c(173530236L, 63864273L)), 
#'   strand = c("-", "-"),
#'   name = c("HAND2", "GATA5"))
#' getInteractionMatrix(a, b)
#'
#' @export
getInteractionMatrix <- function(a, b, ext = 500, count = FALSE) {
  stopifnot("a must be an instance of class 'GRanges'" = is(a, "GRanges"))
  stopifnot("a must have 'name' attribute in it's metadata" = "name" %in% colnames(GenomicRanges::elementMetadata(a)))
  stopifnot("b must be an instance of class 'GRanges'" = is(b, "GRanges"))
  stopifnot("b must have 'name' attribute in it's metadata" = "name" %in% colnames(GenomicRanges::elementMetadata(b)))
  stopifnot("ext must be numeric larger or equal to 0" = is.numeric(ext) & ext >= 0)
  stopifnot("count must be TRUE or FALSE" = isTRUEorFALSE(count))
  
  GenomicRanges::start(a) <- GenomicRanges::start(a) - ext
  GenomicRanges::end(a) <- GenomicRanges::end(a) + ext
  
  hits <- GenomicRanges::findOverlaps(a, b)
  i <- factor(x = a$name[hits@from], levels = unique(a$name))
  j <- factor(x = b$name[hits@to], levels = unique(b$name))
  overlap_mat <- Matrix::sparseMatrix(
    i = as.integer(i),
    j = as.integer(j),
    x = 1,
    dims = c(nlevels(i), nlevels(j)),
    dimnames = list(levels(i), levels(j))
  )
  
  if (! count) {
    attr(overlap_mat, "x") <- rep(1, length(attr(overlap_mat, "x")))
  }
  
  return(overlap_mat)
}
