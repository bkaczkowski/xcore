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
  warning("deprecated! use IRanges::subsetByOverlaps instead")
  stopifnot(is(a, "GRanges"))
  stopifnot(is(b, "GRanges"))

  IRanges::subsetByOverlaps(x = a, ranges = b, ...)
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

#' Normalize the count table with edgeR
#'
#' @param counts expression table with counts
#' @param method passed to to edgeR::calcNormFactors see ?edgeR::calcNormFactors
#' @param log passed on to edgeR::cpm , should the normalized counts be log2 transformed (default FALSE)
#' @param prior.count average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
#' @param ... other parameter passed on to edgeR::calcNormFactors or edgeR::cpm
#' @return normalized count table
#' @import edgeR
#' @export
normalize_counts = function(counts,
                            method = "RLE",
                            log = FALSE ,
                            prior.count = 1 ,
                            ...) {
  dge = edgeR::DGEList(counts)
  dge = edgeR::calcNormFactors(dge , method = method , ...)
  cpm = edgeR::cpm(dge , log = log, prior.count = prior.count, ...)
  cpm
}

#' Select top variable rows
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_var = function(matrix, n = 30) {
  act_var <- apply(matrix, 1, sd)
  matrix <- matrix[order(act_var, decreasing = TRUE), ]
  matrix[1:n, ]
}

#' Select rows with highest mean
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_mean_up = function(matrix,  n = 30) {
  act_mean <- apply(matrix, 1, mean)
  matrix <- matrix[order(act_mean, decreasing = TRUE), ]
  matrix[1:n, ]
}

#' Select rows with lowest mean
#'
#' @param matrix a matrix or data.frame.
#' @param n a number of rows to return.
#'
#' @return a matrix or data.frame with \code{n} rows.
#'
#' @export
select_top_mean_down = function(matrix,  n = 30) {
  act_mean <- apply(matrix, 1, mean)
  matrix <- matrix[order(act_mean, decreasing = FALSE), ]
  matrix[1:n, ]
}

#' Process count matrix for modeling
#'
#'
prepareCountsForRegression <- function(mat,
                                       base_idx,
                                       group = NULL,
                                       log2 = TRUE,
                                       pseudo_count = 1) {
  deglist <- edgeR::DGEList(counts = mat, group = group)
  keep <- edgeR::filterByExpr(deglist)
  deglist <- deglist[keep, , keep.lib.sizes=FALSE] # if not keep some genes have 0 sd
  deglist <- edgeR::calcNormFactors(deglist)
  cpm <- edgeR::cpm(deglist)

  if (is.logical(base_idx)) base_idx <- which(base_idx)

  if (log2) {
    cpm <- log2(cpm + pseudo_count)
  }

  U <- matrix(
    data = rowMeans(cpm[, base_idx]),
    ncol = 1L,
    dimnames = list(rownames(cpm), "u"))

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(U = U, Y = cpm))
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

#' Estimate a goodness of fit stat
estimateStat <- function(x, y, u, s, method = "cv", nfold = 10, statistic = R2, alpha = 0) {
  if (method == "cv") {
    out <- c()
    part <- sample(1:nfold, size = length(y), replace = TRUE)

    for (p in 1:nfold) {
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

  out
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
#' design <- matrix(data = c(1, 1, 0, 0, 0, 0, 1, 1),
#'                  nrow = 4,
#'                  ncol = 2,
#'                  dimnames = list(c(paste("sample", 1:4)), c("gr1", "gr2")))
#' design2factor(design)
#' 
design2factor <- function(design) {
  # based on edgeR::designAsFactor, but jokes aside
  groups <- factor(rowMeans(design * col(design) * ncol(design)))
  groups <- groups[groups != 0] # omit empty groups
  levels(groups) <- colnames(design)
  
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
