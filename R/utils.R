#' re-export magrittr pipe operator
#'
#' @name %>%
#' @rdname pipe
#'
#' @importFrom magrittr %>%
#'
#' @export
NULL

# #' Intersect GRanges objects
# #'
# #' \code{intersectGR} intersect two GRanges objects and return ranges from
# #' \code{a} that were overlapped by ranges in \code{b}.
# #'
# #' @param a GRanges object.
# #' @param b GRanges object.
# #' @param ... other arguments internally passed to
# #'   \code{\link[GenomicRanges]{findOverlaps}}.
# #'
# #' @return GRanges object which is a subset of \code{a}.
# #'
# intersectGR <- function(a, b, ...) {
#   warning("deprecated! use IRanges::subsetByOverlaps instead") # TODO drop from the package
#   stopifnot(is(a, "GRanges"))
#   stopifnot(is(b, "GRanges"))
#
#   IRanges::subsetByOverlaps(x = a, ranges = b, ...)
# }

#' Calculate regions coverage
#'
#' \code{getCoverage} calculates coverage of regions (rows in interaction
#' matrix) by features (columns). It is possible to specify features grouping
#' variable \code{gr} then coverage tells how many distinct groups the region
#' overlap with.
#'
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}}.
#' @param gr factor specifying features groups. Must have length equal
#'   to number of columns in \code{mat}.
#'
#' @return Numeric vector.
#'
#' @examples
#' data("remap_mini")
#' y <- colnames(remap_mini)
#'
#' # simple coverage
#' gr <- seq_along(y) %>% as.factor()
#' getCoverage(remap_mini, gr)
#'
#' # per cell type coverage
#' gr <- sub(".*\\.", "", y) %>% as.factor()
#' getCoverage(remap_mini, gr)
#'
#' @importFrom DelayedArray colsum
#'
#' @export
getCoverage <- function(mat, gr) {
  stopifnot("mat must be an instance of class 'dgCMatrix'" = is(mat, "dgCMatrix"))
  stopifnot("gr must be a factor" = is.factor(gr))
  stopifnot("gr length must equal to number of columns in mat" = length(gr) == ncol(mat))

  DelayedArray::colsum(x = mat, group = gr) %>%
    `>=`(1) %>%
    rowSums()
}

#' Simplify Interaction Matrix
#'
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}}.
#' @param alpha Number between 0 and 1 specifying voting threshold. Eg. for 3
#'   column matrix alpha 0.5 will give voting criteria >= 2.
#' @param colname character giving new column name.
#'
#' @return dgCMatrix
#'
#' @importFrom Matrix rowSums sparseMatrix
#'
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

#' Apply function over groups of columns
#'
#' Returns a array obtained by applying a function to rows of submatrices of the
#' input matrix, where the submatrices are divided into specified groups of
#' columns.
#'
#' @param mat a matrix.
#' @param groups a vector giving columns grouping.
#' @param f function to be applied.
#' @param ... optional arguments to \code{f}.
#'
#' @return a matrix of dimensions \code{nrow(mat) x nlevels(groups)}.
#'
applyOverColumnGroups <- function(mat, groups, f, ...) {
  stopifnot(is.matrix(mat))
  stopifnot(is.vector(groups) || is.factor(groups))
  stopifnot(length(groups) == ncol(mat))
  stopifnot(is.function(f))

  if (is.vector(groups)) {
    groups <- factor(groups)
  }

  group_f <- lapply(X = levels(groups),
                    FUN = function(g) {
                      apply(X = mat,
                            MARGIN = 1,
                            FUN = function(row) f(row[groups == g], ...))
                    })
  new_mat <- do.call(cbind, group_f)
  colnames(new_mat) <- levels(groups)

  return(new_mat)
}

#' Estimate linear models goodness of fit statistic
#'
#' Estimate goodness of fit statistic of penalized linear regression models.
#' Works with different goodness of fit statistic functions.
#'
#' @inheritParams glmnet::glmnet
#' @param u offset vector as in \code{\link{glmnet}}. \code{"U"} experiment in
#'   mae.
#' @param s user supplied lambda.
#' @param method currently only cross-validation is implemented.
#' @param nfold number of fold to use in cross-validation.
#' @param statistic function computing goodness of fit statistic. Should accept
#'   \code{y}, \code{x}, \code{offset} arguments and return
#'   a numeric vector of the same length. See \code{rsq}, \code{mse} for examples.
#'
#' @return numeric vector of \code{statistic} estimates.
#'
#' @importFrom foreach foreach
#'
estimateStat <- function(x, y, u, s, method = "cv", nfold = 10, statistic = rsq, alpha = 0) {
  if (method == "cv") {
    out <- c()
    part <- sample(seq_len(nfold), size = length(y), replace = TRUE)

    out <- foreach::foreach(p_ = seq_len(nfold), .combine = c) %do%
      {
        py <- y[part != p_]
        px <- x[part != p_, ]
        poffset <- u[part != p_, ]
        mod <- glmnet::glmnet(x = px, y = py, offset = poffset, lambda = s, alpha = alpha)

        # evaluate on held-out fold
        py <- y[part == p_]
        px <- x[part == p_, ]
        poffset <- u[part == p_, ]
        yhat <- stats::predict(mod, newx = px, newoffset = poffset, s = s)
        stat <- statistic(py, yhat, px, poffset)
        out <- c(out, stat)
      }
  }

  return(out)
}

# declare estimateStat foreach variables
utils::globalVariables("p_")

#' Calculate $R^2$
#'
#' @param y numeric vector of observed expression values.
#' @param yhat numeric vector of predicted expression values.
#' @param offset numeric vector giving basal expression level.
#'
#' @return numeric vector
#'
rsq <- function(y, yhat, offset) {
  y <- y - offset
  yhat <- yhat - offset
  1 - (sum((y - yhat)^2) / (var(y) * (length(y) - 1)))
}

#' Calculate Mean Squared Error
#'
#' @param y numeric vector of observed expression values.
#' @param yhat numeric vector of predicted expression values.
#' @param ... not used.
#'
#' @return numeric vector
#'
mse <- function(y, yhat, ...) mean((y - yhat)^2)

#' Calculate Mean Absolute Error
#'
#' @param y numeric vector of observed expression values.
#' @param yhat numeric vector of predicted expression values.
#' @param ... not used.
#'
#' @return numeric vector
#'
mae <- function(y, yhat, ...) mean(abs(y - yhat))

#' Transform design matrix to factor
#'
#' @param design design matrix
#'
#' @return factor
#'
#' @examples
#' \dontrun{
#' design <- matrix(data = c(1, 1, 0, 0, 0, 0, 1, 1),
#'                  nrow = 4,
#'                  ncol = 2,
#'                  dimnames = list(c(paste("sample", 1:4)), c("gr1", "gr2")))
#' design2factor(design)
#' }
#'
design2factor <- function(design) {
  # based on edgeR::designAsFactor, but jokes aside
  groups <- factor(rowMeans(design * col(design) * ncol(design)))
  samples_to_keep <- groups != 0
  groups <- groups[samples_to_keep] # omit empty groups
  groups <- droplevels(groups)
  levels(groups) <- colnames(design)[colSums(design) > 0]

  return(groups)
}

#' Check if argument is a binary flag
#'
#' @param x object to test
#'
#' @return binary flag
#'
isTRUEorFALSE <- function(x) {
  (length(x) == 1) && is.logical(x) && (! is.na(x))
}

#' Apply function over selected column in list of data frames
#'
#' \code{applyOverDFList} operates on a list of data frames where all data frames
#' has the same size and columns. Column of interest is extracted from each data
#' frame and column binded in \code{groups}, next \code{fun} is applied over
#' rows. Final result is a matrix with result for each group on a separate column.
#' Function is parallelized over groups.
#'
#' @param list_of_df list of \code{data.frame}s.
#' @param col_name string specifying column in \code{data.frame}s to apply
#'   \code{fun} on.
#' @param fun function to apply, should take a single vector as a argument.
#' @param groups factor defining how elements of \code{list_of_df} should be
#'   grouped.
#'
#' @return matrix with \code{nrow(list_of_df[[1]])} rows and
#'   \code{nlevels(groups)} columns.
#'
applyOverDFList <- function(list_of_df, col_name, fun, groups) {
  stopifnot("all list_of_df names must be included in groups" = setequal(names(list_of_df), names(groups)))
  stopifnot("groups must not have unused levels" = setdiff(levels(groups), groups) == character(0))

  col_fun_mat <- foreach::foreach(gr_ = levels(groups), .combine = cbind) %dopar% # parallel version is bit faster
    {
      i <- groups == gr_
      attr <- lapply(list_of_df[i], function(df) df[[col_name]])
      attr <- do.call(cbind, attr)
      matrix(
        data = apply(X = attr, MARGIN = 1, FUN = fun),
        ncol = 1L,
        dimnames = list(rownames(attr), gr_)
      )
    }

  return(col_fun_mat)
}

# declare applyOverDFList foreach variables
utils::globalVariables("gr_")

#' Subset keeping missing
#'
#' Subset matrix keeping unmatched rows as NA.
#'
#' @param mat matrix
#' @param rows character
#'
#' @return a matrix
#'
subsetWithMissing <- function(mat, rows) {
  i <- match(x = rows, table = rownames(mat), nomatch = 0)
  i[i == 0] <- NA
  smat <- mat[i, ]
  rownames(smat) <- rows
  smat
}

#' Helper summarizing MAE object
#'
#' @param mae MultiAssayExperiment object.
#'
#' @return named list giving number of rows and columns, overall mean and
#'   standard deviation in \code{mae}'s experiments.
#'
#' @importFrom MultiAssayExperiment experiments
#' @importFrom Matrix mean
#' @importFrom stats sd
#'
#'
maeSummary <- function(mae) {
  lapply(MultiAssayExperiment::experiments(mae), function(x) c(dim(x), Matrix::mean(x), sd(x)))
}


#' Translate counts matrix rownames
#'
#' \code{translateCounts} renames counts matrix rownames according to supplied
#' \code{dict}ionary. Function can handle many to one assignments by summing
#' over \code{counts} rows. Other types of ambiguous assignments are not
#' supported.
#'
#' @param counts matrix of expression values.
#' @param dict named character vector mapping \code{counts} rownames to new
#'   values.
#'
#' @return matrix of expression values with new rownames.
#'
#' @examples
#' #TODO also add tests
#'
#' @export
translateCounts <- function(counts, dict) {
  stopifnot("counts must be an instance of class 'matrix'" = is(counts, "matrix"))
  stopifnot("counts must have its rownames defined" = ! is.null(rownames(counts)))
  stopifnot("counts rownames must be unique" = ! anyDuplicated(rownames(counts)))
  stopifnot("dict must be an instance of class 'character'" = is(dict, "character"))
  stopifnot("dict must be a named character vector" = ! is.null(names(dict)))
  stopifnot("dict names must be unique" = ! anyDuplicated(names(dict)))
  stopifnot("counts rownames must correspond to dict names" = any(rownames(counts) %in% names(dict)))

  new_counts <- counts[rownames(counts) %in% names(dict), ]
  new_counts <- rowsum(x = new_counts, group = dict[rownames(new_counts)]) # TODO a simple and useful extension would be to implement option to take an average instead of sum, this could be easly done by dividing by the groups size

  return(new_counts)
}
