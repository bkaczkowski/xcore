#' re-export magrittr pipe operator
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

#' Intersect GRanges objects
#'
#' \code{intersectGR} intersect two GRanges objects and return ranges from
#' \code{a} that were overlapped by ranges in \code{b}.
#'
#' @param a GRanges object.
#' @param b GRanges object.
#' @param ... other arguments internally passed to
#'   \code{\link[GenomicRanges]{findOverlaps}}.
#'
#' @return GRanges object which is a subset of \code{a}.
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#'
#' @export
intersectGR <- function(a, b, ...) {
  stopifnot(is(a, "GRanges"))
  stopifnot(is(b, "GRanges"))

  hits <- GenomicRanges::findOverlaps(query = b, subject = a, ...)
  i <- S4Vectors::subjectHits(hits)

  return(a[i, ])
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
#'   between \code{a} and \code{b} or only presence / absence indicators.
#'
#' @return Sparse matrix of class dgCMatrix, with rows corresponding to
#'   \code{a} and columns to \code{b}. Each cell holds a number indicating
#'   how many times \code{a} and \code{b} overlapped.
#'
#' @importFrom GenomicRanges start end findOverlaps
#' @importFrom Matrix sparseMatrix
#'
#' @export
getInteractionMatrix <- function(a, b, ext = 500, count = FALSE) {
  stopifnot(is(a, "GRanges"))
  stopifnot(is(b, "GRanges"))
  stopifnot(is.numeric(ext) & ext >= 0)
  stopifnot(is.logical(count))

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

#' Calculate Tau
#'
#' \code{tau} calculate Tau tissue-specificity metric.
#' \deqn{\tau = \frac{\sum_{i=1}^{n} (1-\hat{x_i})}{n-1}; \hat{x_i} = \frac{x_i}{max(x)}}
#' \link[https://academic.oup.com/bib/article/18/2/205/2562739]{Reference}.
#'
#' @param x Numeric vector giving expression of one gene in different tissues.
#'
#' @return Tau tissue-specificity metric.
#'
tau <- function(x) {
  stopifnot(is.numeric(x) && length(x) > 1)

  if (all(!is.na(x)))
  {
    if (min(x, na.rm = TRUE) >= 0)
    {
      if (max(x) != 0)
      {
        x <- (1 - (x / max(x)))
        res <- sum(x, na.rm = TRUE)
        res <- res / (length(x) - 1)
      } else {
        res <- 0
      }
    } else {
      stop("Expression values have to be positive!")
    }
  } else {
    stop("No data for this gene avalable.")
  }

  return(res)
}

#' Calculate regions coverage
#'
#' \code{getCoverage} calculates coverage of regions (rows in interaction
#' matrix) by features (columns). It is possible to specify features grouping
#' variable \code{gr} then coverage tells how many distinct groups the region
#' overlap with.
#'
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}.
#' @param gr Character vector specifying features groups. Must have length equal
#'   to number of columns in \code{mat}.
#'
#' @return Numeric vector.
#'
#' @importFrom DelayedArray colsum
#'
#' @export
getCoverage <- function(mat, gr) {
  stopifnot(is(mat, "dgCMatrix"))
  stopifnot(is.character(gr))
  stopifnot(length(gr) == ncol(mat))

  group <- as.factor(gr)
  DelayedArray::colsum(x = mat, group = group) %>%
    `>`(1) %>%
    rowSums()
}

#'
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}.
#' @param gr Character vector specifying features groups. Must have length equal
#'   to number of columns in \code{mat}.
#' @param method String. Currently only majority voting strategy is available.
#'
#' @return dgCMatrix
#'
simplifyInteractionMatrix <- function(mat, gr, method = "majority") {
  stopifnot(is(mat, "dgCMatrix"))
  stopifnot(is.character(rownames(mat)))
  stopifnot(is.character(gr))
  stopifnot(length(gr) == ncol(mat))

  group <- as.factor(gr)
  i.list <- lapply(X = levels(group),
         FUN = function(x) {
           m <- group == x
           crit <- floor(sum(m) / 2)
           mask <- Matrix::rowSums(mat[, m, drop = FALSE]) > crit
           which(mask)
         })
  i_lens <- vapply(i.list, length, numeric(1L))
  j <- rep(seq(1, nlevels(group)), times = i_lens)
  Matrix::sparseMatrix(i = unlist(i.list, use.names = FALSE),
                       j = j,
                       x = 1,
                       dims = c(nrow(mat), nlevels(group)),
                       dimnames = list(rownames(mat), levels(group)))
}

# TODO test
# identical(
#   simplifyInteractionMatrix(Matrix(
#     matrix(
#       c(1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1),
#       ncol = 6,
#       byrow = T,
#       dimnames = list(1:3)
#     ), sparse = TRUE
#   ),
#   c("A", "A", "A", "B", "B", "B")),
#   new(
#     "dgCMatrix",
#     i = c(0L, 2L, 1L),
#     p = c(0L, 2L, 3L),
#     Dim = 3:2,
#     Dimnames = list(c("1", "2", "3"), c("A", "B")),
#     x = c(1,
#           1, 1),
#     factors = list()
#   )
# )
