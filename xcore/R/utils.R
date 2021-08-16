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
#'   \code{\link{getInteractionMatrix}}.
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

#' Simplify Interaction Matrix
#'
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}}.
#' @param alpha Number between 0 and 1 specifying voting threshold. Eg. for 3
#'   column matrix alpha 0.5 will give voting criteria >= 2.
#'
#' @return dgCMatrix
#'
#' @importFrom Matrix rowSums sparseMatrix
#'
#' @export
simplifyInteractionMatrix <- function(mat, alpha = 0.5, colname = NA) {
  stopifnot(is(mat, "dgCMatrix"))
  stopifnot(is.numeric(alpha) && (length(alpha) == 1))
  stopifnot((alpha > 0) && (alpha <= 1))

  crit <- ceiling(ncol(mat) * alpha)
  i <- which(Matrix::rowSums(mat) >= crit)
  Matrix::sparseMatrix(i = i,
                       j = rep(1L, length(i)),
                       x = 1,
                       dims = c(nrow(mat), 1L),
                       dimnames = list(rownames(mat), colname))
}

#' Calculate vector purity
#'
#' @export
vectorPurity <- function(x) {
  xtab <- table(x, useNA = "always")
  max(xtab) / length(x)
}

#' Get most frequent value
#'
#'
mostFreqValue <- function(vec) {
  freq <- table(vec, useNA = "always")
  freq <- sort(freq, decreasing = TRUE)
  mostfreq <- names(freq)[1]

  return(mostfreq)
}

#' Prune metadata
#'
#' Prunes metadata based on selected column such that only entries with most
#' frequent items remains.
#'
#' @param col pruning column
#'
#' @export
pruneClusterMeta <- function(meta, col) {
  idx <- meta[[col]]
  mostfreq <- mostFreqValue(idx)
  mask <- idx == mostfreq

  return(meta[mask, ])
}

#' Make cluster name based on metadata
#'
#'
makeClusterName <- function(meta,
                            columns = c("tf", "biotype", "study", "tf_dbd", "cov_type"),
                            mixed_th = 0.5) {
  nms <- vapply(X = columns,
                FUN = function(x) {
                  vec <- meta[[x]]
                  nm <-
                    ifelse(vectorPurity(vec) > mixed_th, mostFreqValue(vec), "mixed")
                  paste(x, nm, sep = ".")
                },
                FUN.VALUE = character(1L))
  nm <- paste(nms, collapse = "_")

  return(nm)
}

#'
#' @param pruning_purity if purity > pruning_purity then prune
#'
collapseInteractionMatrix <- function(mat,
                                      meta,
                                      cl,
                                      alpha = 0.5,
                                      purity_feature = "tf",
                                      min_pruning_purity = 0.5
) {
  stopifnot(is(mat, "dgCMatrix"))
  stopifnot(is(meta, "data.table"))
  stopifnot("id" %in% colnames(meta))
  stopifnot(all(colnames(mat) == meta$id))
  stopifnot(length(cl) == ncol(mat))

  cl <- as.factor(cl)
  out <- list()
  for (clu in levels(cl)) {
    m <- cl == clu
    purity <- vectorPurity(meta[m, ][[purity_feature]])

    if (purity > min_pruning_purity) {
      pruned_meta <- pruneClusterMeta(meta[m, ], purity_feature)
      m <- pruned_meta[["id"]]
      nice_name <- makeClusterName(pruned_meta, mixed_th = min_pruning_purity)
    } else {
      nice_name <- makeClusterName(meta[m, ], mixed_th = min_pruning_purity)
    }
    nice_name <- paste(clu, nice_name, sep = "_")

    # alpha <- 0.3 # calculate alpha based on

    out[[nice_name]] <- simplifyInteractionMatrix(mat = mat[, m],
                                                  alpha = alpha,
                                                  colname = nice_name)
  }

  new_mat <- do.call(cbind, out)

  return(new_mat)
}

#' Select cutHeight and deepSplit parameters for dynamicTreeCut
#'
#'
selectParams4dynamicTreeCut <- function(dendro,
                                        distM,
                                        ref_cl,
                                        minClusterSize = 3,
                                        c_deepSplit = 1:4,
                                        c_cutHeight = seq(from = 0.5, to = 6, by = 0.2)) {
  res <- list(deepSplit = c(), cutHeight = c(), ARI = c(), AMI = c(), N = c())
  for (ds in c_deepSplit) {
    for (ch in c_cutHeight) {
      cl <- dynamicTreeCut::cutreeDynamic(
        dendro = dendro,
        cutHeight = ch,
        minClusterSize = minClusterSize,
        method = "hybrid",
        distM = distM,
        deepSplit = ds,
        pamStage = TRUE)


      res[["deepSplit"]] <- c(res[["deepSplit"]], ds)
      res[["cutHeight"]] <- c(res[["cutHeight"]], ch)
      res[["AMI"]] <- c(res[["AMI"]], aricode::AMI(ref_cl, cl))
      res[["ARI"]] <- c(res[["ARI"]], aricode::ARI(ref_cl, cl))
      res[["N"]] <- c(res[["N"]], cl %>% unique() %>% length())
    }
  }

  res <- do.call(cbind, res) %>% as.data.table()

  return(res)
}

#' Paint vector
#'
#' Get colors along vector
#'
paintVector <- function(vec) {
  colorspace::qualitative_hcl(length(unique(vec)))[as.integer(factor(vec))]
}
