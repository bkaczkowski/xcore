#' re-export magrittr pipe operator
#'
#' @name %>%
#' @rdname pipe
#'
#' @importFrom magrittr %>%
#'
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
intersectGR <- function(a, b, ...) {
  warning("deprecated! use IRanges::subsetByOverlaps instead") # TODO drop from the package
  stopifnot(is(a, "GRanges"))
  stopifnot(is(b, "GRanges"))

  IRanges::subsetByOverlaps(x = a, ranges = b, ...)
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
tau <- function(x) { # TODO likely to be removed from the package
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
#' @param gr factor specifying features groups. Must have length equal
#'   to number of columns in \code{mat}.
#'
#' @return Numeric vector.
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
                            columns = c("tf", "biotype", "study", "tf_dbd"),
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

#' Make cluster name based on metadata ver. 2
#'
#'
makeClusterName2 <- function(meta,
                             columns = c("tf", "biotype", "study", "tf_dbd"),
                             inclusion_th = 0.25,
                             mixed_th = 0.5) {
  nms <- vapply(X = columns,
                FUN = function(x) {
                  vec <- meta[[x]]
                  ifreq <- table(vec) / length(vec)
                  ifreq <- ifreq[ifreq >= inclusion_th]
                  if (length(ifreq)) {
                    ifreq <- sort(ifreq, decreasing = TRUE)
                    nms <- names(ifreq)
                    nms <- nms[nms != "" | is.na(nms)] # drop empty names
                    if (all(ifreq < mixed_th)) nms <- c("mixed", nms) # if none dominates add mixed prefix
                    if (! length(nms)) nms <- "mixed" # if all names were empty
                    nm <- paste(nms, collapse = "_")
                  } else {
                    nm <- "mixed"
                  }
                  paste(x, nm, sep = ".")
                },
                FUN.VALUE = character(1L))
  nms <- nms[! grepl(".*\\.mixed$", nms)]
  nm <- if (length(nms)) { paste(nms, collapse = ".") } else { "mixed" }

  return(nm)
}

#' collapseInteractionMatrix
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
      nice_name <- makeClusterName2(pruned_meta, mixed_th = min_pruning_purity)
    } else {
      nice_name <- makeClusterName2(meta[m, ], mixed_th = min_pruning_purity)
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

#' Estimate a goodness of fit stat
estimateStat <- function(x, y, u, s, method = "cv", nfold = 10, statistic = R2, alpha = 0) {
  if (method == "cv") {
    out <- c()
    part <- sample(1:nfold, size = length(y), replace = TRUE)

    out <- foreach::foreach(p = seq_len(nfold), .combine = c) %do% # TODO parallelization should be optional maybe?
      {
        py <- y[part != p]
        px <- x[part != p, ]
        poffset <- u[part != p, ]
        mod <- glmnet::glmnet(x = px, y = py, offset = poffset, lambda = s, alpha = alpha)

        # evaluate on held-out fold
        py <- y[part == p]
        px <- x[part == p, ]
        poffset <- u[part == p, ]
        yhat <- predict(mod, newx = px, newoffset = poffset, s = s)
        stat <- statistic(py, yhat, px, poffset)
        out <- c(out, stat)
      }
  }

  return(out)
}

# https://stats.stackexchange.com/questions/186396/appropriate-way-to-calculate-cross-validated-r-square
# 1-(sum((data[,1] - predictions)^2) / ((n-1) * var(test[,1]))
rsq <- function(y, yhat, x, offset) {
  y <- y - offset
  yhat <- yhat - offset
  1 - (sum((y - yhat)^2) / (var(y) * (length(y) - 1)))
}
pearson.sq <- function(y, yhat, x, offset) cor(y - offset, yhat - offset)^2
mse <- function(y, yhat, ...) mean((y - yhat)^2)
mae <- function(y, yhat, ...) mean(abs(y - yhat))
spearman.sq <- function(y, yhat, x, offset) cor(y - offset, yhat - offset, method = "spearman")^2

#'
getAvgCoeff <- function(models, lambda = "lambda.min", drop_intercept = TRUE) {
  coefs <- lapply(models, function(m) coef(m, s = m[[lambda]]))
  if (drop_intercept) {
    coefs <- lapply(coefs, function(m) {
      keep <- grep(pattern = "(Intercept)", x = rownames(m), invert = TRUE)
      m[keep, ]
    })
  }
  coefs <- do.call(cbind, coefs)
  coefs_avg <- rowMeans(coefs)
  coefs_sd <- apply(coefs, 1, sd)
  res <-
    cbind(
      estimate = coefs_avg,
      sd = coefs_sd,
      z = (coefs_avg - mean(coefs_avg)) / coefs_sd)
  res
}

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
  # TODO consider adding argument check
  # stopifnot("groups must not have unused levels" = setdiff(levels(groups), groups) == character(0))

  col_fun_mat <- foreach::foreach(gr = levels(groups), .combine = cbind) %dopar% # parallel version is bit faster
    {
      i <- groups == gr
      attr <- lapply(list_of_df[i], function(df) df[[col_name]])
      attr <- do.call(cbind, attr)
      matrix(
        data = apply(X = attr, MARGIN = 1, FUN = fun),
        ncol = 1L,
        dimnames = list(rownames(attr), gr)
      )
    }

  return(col_fun_mat)
}
