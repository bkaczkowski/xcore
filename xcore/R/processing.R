#' Process count matrix for expression modeling
#'
#' Expression counts are processed using \link[edgeR]{edgeR} following
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf}{User's Guide}.
#' Shortly, counts for each sample are filtered for lowly expressed promoters,
#' normalized for the library size and transformed into counts per million (CPM).
#' Optionally, CPM are log2 transformed with addition of pseudo count. Basal
#' level expression is calculated by averaging \code{base_lvl} samples
#' expression values.
#'
#' @param counts matrix of read counts.
#' @param design matrix giving the design matrix for the samples. Columns
#'   corresponds to samples groups and rows to samples names.
#' @param base_lvl string indicating group in \code{design} corresponding to
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
#'   design with \code{base_lvl} dropped is stored in metadata and directly
#'   available for \code{modelGeneExpression}.
#'
#' @examples
#' base_lvl <- "00hr"
#' design <- matrix(
#'   data = c(1, 0, 0,
#'            1, 0, 0,
#'            1, 0, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 0, 1,
#'            0, 0, 1,
#'            0, 0, 1),
#'   ncol = 3,
#'   nrow = 9,
#'   byrow = TRUE,
#'   dimnames = list(colnames(rinderpest_mini), c("00hr", "12hr", "24hr")))
#' mae <- prepareCountsForRegression(
#'   counts = rinderpest_mini,
#'   design = design,
#'   base_lvl = base_lvl)
#'
#' @importFrom edgeR DGEList calcNormFactors cpm filterByExpr
#' @importFrom MultiAssayExperiment ExperimentList MultiAssayExperiment
#'
#' @export
prepareCountsForRegression <- function(counts,
                                       design,
                                       base_lvl,
                                       log2 = TRUE,
                                       pseudo_count = 1L,
                                       drop_base_lvl = TRUE) {
  stopifnot("counts must be an integer matrix" = is.matrix(counts) && (typeof(counts) == "integer"))
  stopifnot("design must be a matrix" = is.matrix(design)) # && (typeof(design) == "integer"))
  stopifnot("number of rows in design must equal to number of columns in counts" = nrow(design) == ncol(counts))
  stopifnot("design rownames must be the same as counts colnames" = all(rownames(design) == colnames(counts)))
  stopifnot("base_lvl must match to one of design colnames" = sum(colnames(design) == base_lvl) == 1L)
  stopifnot("log2 must be TRUE or FALSE" = isTRUEorFALSE(log2))
  stopifnot("pseudo_count must be an positive integer or zero" = is.integer(pseudo_count) && (length(pseudo_count) == 1L) && (pseudo_count >= 0L))
  stopifnot("drop_base_lvl must be TRUE or FALSE" = isTRUEorFALSE(drop_base_lvl))

  groups <- design2factor(design)
  deglist <- edgeR::DGEList(counts = counts, group = groups)
  keep <- edgeR::filterByExpr(deglist)
  deglist <- deglist[keep, , keep.lib.sizes=FALSE] # if not keep some genes have 0 sd
  deglist <- edgeR::calcNormFactors(deglist)
  cpm <- edgeR::cpm(deglist)

  if (log2) {
    cpm <- log2(cpm + pseudo_count)
  }

  U <- matrix(
    data = rowMeans(cpm[, groups == base_lvl, drop = FALSE]),
    ncol = 1L,
    dimnames = list(rownames(cpm), "u"))

  ymask <- if (drop_base_lvl) { groups != base_lvl } else { rep(TRUE, length(groups)) }
  design[, base_lvl] <- 0 # design without base_lvl
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(
      U = U,
      Y = cpm[, ymask, drop = FALSE]),
    metadata = list(design = design))
}

#' Add molecular signatures to MultiAssayExperiment
#'
#' \code{addSignatures} extends \code{mae} by adding to it new experiments.
#' Rows consistency is ensured by taking an intersection of rows after new
#' experiments are added.
#'
#' @param mae MultiAssayExperiment object.
#' @param ... named experiments to be added to \code{mae}.
#' @param intersect_rows logical flag indicating if only common rows across
#'   experiments should be included. Only set to \code{FALSE} if you know what
#'   you are doing.
#'
#' @return MultiAssayExperiment object with new experiments added.
#'
#' @examples
#' base_lvl <- "00hr"
#' design <- matrix(
#'   data = c(1, 0, 0,
#'            1, 0, 0,
#'            1, 0, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 0, 1,
#'            0, 0, 1,
#'            0, 0, 1),
#'   ncol = 3,
#'   nrow = 9,
#'   byrow = TRUE,
#'   dimnames = list(colnames(rinderpest_mini), c("00hr", "12hr", "24hr")))
#' mae <- prepareCountsForRegression(
#'   counts = rinderpest_mini,
#'   design = design,
#'   base_lvl = base_lvl)
#' mae <- addSignatures(mae, remap = remap_mini)
#'
#' @importFrom IRanges SplitDataFrameList
#' @importFrom methods is
#' @importFrom MultiAssayExperiment experiments intersectRows listToMap mapToList sampleMap
#' @importFrom S4Vectors DataFrame
#'
#' @export
addSignatures <- function(mae, ..., intersect_rows = TRUE) {
  stopifnot("mae must be an instance of class 'MultiAssayExperiment'" = is(mae, "MultiAssayExperiment"))
  stopifnot("intersect_rows must be TRUE or FALSE" = isTRUEorFALSE(intersect_rows))
  ex <- list(...)
  ex_names <- names(ex)
  stopifnot("experiments must be named" = ! (is.null(ex_names) || any(ex_names == "")))
  stopifnot("experiments names must be unique" = ! anyDuplicated(c(names(mae), ex_names)))

  mae_map <- MultiAssayExperiment::sampleMap(mae)
  ex_map <-
    lapply(ex, function(x)
      S4Vectors::DataFrame(primary = colnames(x), colname = colnames(x)))
  names(ex_map) <- ex_names
  newListMap <- c(MultiAssayExperiment::mapToList(mae_map),
                  IRanges::SplitDataFrameList(ex_map))
  newSampleMap <- MultiAssayExperiment::listToMap(newListMap)

  newExperimentList <- c(MultiAssayExperiment::experiments(mae), ex)

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

#' Filter signatures by coverage
#'
#' Filter signatures overlapping low or high number of promoters. Useful to get
#' rid of signatures that have very low variance.
#'
#' @param mae MultiAssayExperiment object.
#' @param min length one numeric between 0 and 1 defining minimum promoter
#'   coverage for the signature to pass filtering.
#' @param max length one numeric between 0 and 1 defining maximum promoter
#'   coverage for the signature to pass filtering.
#' @param ref_experiment string giving name of experiment to use for inferring
#'   total number of promoters.
#' @param omit_experiments character giving names of experiments to exclude from
#'   filtering.
#'
#' @return MultiAssayExperiment object with selected experiments filtered.
#'
#' @importFrom Matrix colSums
#'
#' @examples
#' base_lvl <- "00hr"
#' design <- matrix(
#'   data = c(1, 0, 0,
#'            1, 0, 0,
#'            1, 0, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 1, 0,
#'            0, 0, 1,
#'            0, 0, 1,
#'            0, 0, 1),
#'   ncol = 3,
#'   nrow = 9,
#'   byrow = TRUE,
#'   dimnames = list(colnames(rinderpest_mini), c("00hr", "12hr", "24hr")))
#' mae <- prepareCountsForRegression(
#'   counts = rinderpest_mini,
#'   design = design,
#'   base_lvl = base_lvl)
#' mae <- addSignatures(mae, remap = remap_mini)
#' mae <- filterSignatures(mae)
#'
#' @export
filterSignatures <- function(mae,
                             min = 0.05,
                             max = 0.95,
                             ref_experiment = "Y",
                             omit_experiments = c("Y", "U")) {
  stopifnot("mae must be an instance of class 'MultiAssayExperiment'" = is(mae, "MultiAssayExperiment"))
  stopifnot("min must be a length one numeric between 0 and 1" = is.numeric(min) && length(min) == 1 && min >= 0 && min <= 1)
  stopifnot("max must be a length one numeric between 0 and 1" = is.numeric(max) && length(max) == 1 && max >= 0 && max <= 1)
  stopifnot("ref_experiment must be a length one character" = is.character(ref_experiment) && length(ref_experiment) == 1)
  stopifnot("ref_experiment must match one of mae names" = ref_experiment %in% names(mae))
  stopifnot("omit_experiments must be a character" = is.character(omit_experiments))
  stopifnot("omit_experiments must match mae names" = all(omit_experiments %in% names(mae)))

  mae_nrow <- nrow(mae[[ref_experiment]])
  signatures <- setdiff(names(mae), omit_experiments)
  for (sig in signatures) {
    if (! is(mae, "sparseMatrix")) {
      warning(sprintf("Only 'sparseMatrix' experiments filtering is supported. Omitting %s", sig))
      next()
    }
    mask <- Matrix::colSums(mae[[sig]] != 0)
    mask <- mask >= (min * mae_nrow) & mask <= (max * mae_nrow)
    suppressMessages(mae[[sig]] <- mae[[sig]][, mask])
  }

  return(mae)
}
