#' Combine p-values using Fisher method
#'
#' Fisher's method is a meta-analysis technique used to combine the results from
#' independent statistical tests with the same hypothesis
#' (\href{https://en.wikipedia.org/wiki/Fisher%27s_method}{Wikipedia article}).
#'
#' @inheritParams stats::pchisq
#' @param p.value a numeric vector of p-values to combine.
#'
#' @return a number giving combined p-value.
#'
#' @importFrom stats pchisq
#'
fisherMethod <- function(p.value, lower.tail = FALSE, log.p = TRUE) {
  stopifnot("p.value must be numeric" = is.numeric(p.value))
  stopifnot("p.value must be longer than 1" = length(p.value) > 1L)
  stopifnot("lower.tail must be TRUE or FALSE" = isTRUEorFALSE(lower.tail))
  stopifnot("log.p must be TRUE or FALSE" = isTRUEorFALSE(log.p))

  K <- 2 * length(p.value)
  X <- -2 * sum(log(p.value))
  cmbp <- stats::pchisq(X, df = K, lower.tail = lower.tail, log.p = log.p)

  return(cmbp)
}

#' Linear ridge regression
#'
#' Wrapper around \code{\link[glmnet]{cv.glmnet}} to run linear ridge regression
#' with lambda selection using cross-validation.
#'
#' @inheritParams glmnet::cv.glmnet
#'
#' @return an object of class "cv.glmnet" is returned. See
#'   \code{\link[glmnet]{cv.glmnet}} for more details.
#'
#' @importFrom glmnet cv.glmnet
#'
runLinearRidge <-
  function(x,
           y,
           offset,
           alpha = 0,
           standardize = TRUE,
           ...) {
    cv <- glmnet::cv.glmnet(
      x = x,
      y = y,
      offset = offset,
      alpha = alpha,
      standardize = standardize,
      ...)

    return(cv)
  }

#' Significance testing in linear ridge regression
#'
#' Standard error estimation and significance testing for coefficients
#' estimated in linear ridge regression. \code{ridgePvals} re-implement
#' original method by (Cule et al. BMC Bioinformatics 2011.) found in
#' \link[ridge]{ridge-package}. This function is intended to use with
#' \code{\link{runLinearRidge}} output.
#'
#' @param x input matrix, same as used in \code{\link{runLinearRidge}}.
#' @param y response variable, same as used in \code{\link{runLinearRidge}}.
#' @param beta matrix of coefficients, estimated using
#'   \code{\link{runLinearRidge}}.
#' @param lambda lambda value for which \code{beta} was estimated.
#' @param standardizex logical flag for x variable standardization, should be
#'   set to same value as \code{standarize} flag in \code{\link{runLinearRidge}}.
#' @param svdX optional singular-value decomposition of \code{x} matrix. One can
#'   be obtained using \code{link[base]{svd}}. Passing this argument omits
#'   internal call to \code{link[base]{svd}}, this is useful when calling
#'   \code{ridgePvals} repeatedly using same \code{x}.
#'
#' @return a data.frame with columns
#'   \describe{
#'     \item{coef}{\code{beta}'s names}
#'     \item{se}{\code{beta}'s standard errors}
#'     \item{tstat}{\code{beta}'s test statistic}
#'     \item{pval}{\code{beta}'s p-values}
#'   }
#'
ridgePvals <- function (x, y, beta, lambda, standardizex = TRUE, svdX = NULL) {
  n <- length(y)
  if (standardizex) x <- scale(x)
  if (is.null(svdX)) svdX <- svd(x)
  U <- svdX$u
  D <- svdX$d
  D2 <- svdX$d ^ 2
  V <- svdX$v
  div <- D2 + lambda
  # estimate sig^2
  sig2hat <-                      # adding ( here helps a lot!    # t(U) %*% y; crossprod is way faster
    as.numeric(crossprod(y - U %*% (diag((D2) / (div)) %*% crossprod(U, y)))) /
    (n - sum(D2 * (D2 + 2 * lambda) / (div ^ 2)))
  varmat <- tcrossprod(V %*% diag(D2 / (div ^ 2)), V)
  varmat <- sig2hat * varmat
  se <- sqrt(diag(varmat))
  tstat <- abs(beta / se)
  pval <- 2 * (1 - pnorm(tstat))

  res <-
    list(
      coef = beta,
      se = se,
      tstat = tstat,
      pval = pval)
  class(res) <- "data.frame"
  attr(res, "row.names") <- colnames(x)

  return(res)
}

#' Gene expression modeling pipeline
#'
#' \code{modelGeneExpression} uses parallelization if parallel backend is
#' registered. For that reason we advise agains passing \code{parallel} argument
#' to internally called \code{\link[glmnet]{cv.glmnet}} routine.
#'
#' For speeding up the calculations consider lowering number of folds used in
#' internally run \code{\link[glmnet]{cv.glmnet}} by specifying \code{nfolds}
#' argument. By default 10 fold cross validation is used.
#'
#' The relationship between the expression (Y) and molecular signatures (X) is
#' described using linear model formulation. The pipeline attempts to model the
#' change in expression between basal expression level (u) and each sample, with
#' the goal of finding the unknown molecular signatures activities. Linear
#' models are fit using popular ridge regression implementation
#' \link[glmnet]{glmnet} (Friedman, Hastie, and Tibshirani 2010).
#'
#' If \code{pvalues} is set to \code{TRUE} the significance of the estimated
#' molecular signatures activities is tested using methodology introduced by
#' (Cule, Vineis, and De Iorio 2011) which original implementation can be found
#' in \link[ridge]{ridge-package}.
#'
#' If replicates are available the signatures activities estimates and
#' their standard error estimates can be combined. This is done by averaging
#' signatures activities estimates and pooling standard errors following
#' (Cohen 1977) formulation for pooled standard deviation.
#'
#' For detailed pipeline description we refer interested user to paper
#' accompanying this package.
#'
#' @param mae MultiAssayExperiment object such as produced by
#'   \code{\link{prepareCountsForRegression}}.# TODO X has to be filtered at least for zeros and ones !!!
#' @param yname string indicating experiment in \code{mae} to use as the
#'   expression input.
#' @param uname string indicating experiment in \code{mae} to use as the basal
#'   expression level.
#' @param xnames character indicating experiments in \code{mae} to use as
#'   molecular signatures.
#' @param design matrix giving the design matrix for the samples. Columns
#'   corresponds to samples groups and rows to samples names. Only samples
#'   included in the design will be processed.
#' @param standardize logical flag indicating if the molecular signatures should
#'   be scaled. Advised to be set to \code{TRUE}.
#' @param parallel parallel argument to internally used
#'   \code{\link[glmnet]{cv.glmnet}} function. Advised to be set to \code{FALSE}
#'   as it might interfere with parallelization used in \code{modelGeneExpression}.
#' @param pvalues logical flag indicating if significance testing for the
#'   estimated molecular signatures activities should be performed.
#' @param precalcmodels optional list of precomputed \code{'cv.glmnet'} objects
#'   for each molecular signature and sample. The elements of this list should
#'   be matching the \code{xnames} vector. Each of those elements should be a
#'   named list holding \code{'cv.glmnet'} objects for each sample. If provided
#'   those models will be used instead of running regression from scratch.
#' @param ... arguments passed to glmnet::cv.glmnet.
#'
#' @return Nested list with following elements
#'   \describe{
#'     \item{regression_models}{Named list with elements corresponding to
#'       signatures specified in \code{xnames}. Each of these is a list holding
#'       \code{'cv.glmnet'} objects corresponding to each sample.}
#'     \item{pvalues}{Named list with elements corresponding to
#'       signatures specified in \code{xnames}. Each of these is a list holding
#'       \code{data.frame} of signature's p-values and test statistics
#'       estimated for each sample.}
#'     \item{replicate_avg}{Named list with elements corresponding to
#'       signatures specified in \code{xnames}. Each of these is a \code{matrix}
#'       holding replicate average Z-scores with columns corresponding to groups
#'       in the design. TODO this output need better description!}
#'   }
#'
#' @examples
#' TODO
#'
#' @importFrom foreach foreach
#' @importFrom iterators iter
#'
#' @export
modelGeneExpression <- function(mae,
                                yname,
                                uname,
                                xnames,
                                design,
                                standardize = TRUE,
                                parallel = FALSE,
                                pvalues = TRUE,
                                precalcmodels = NULL,
                                ...) {
  stopifnot("mae must be an instance of class 'MultiAssayExperiment'" = is(mae, "MultiAssayExperiment"))
  stopifnot("yname must be a length one character" = is.character(yname) && length(yname) == 1)
  stopifnot("uname must be a length one character" = is.character(uname) && length(uname) == 1)
  stopifnot("yname must be distinct from uname" = yname != uname)
  stopifnot("xnames must be a character vector" = is.character(xnames))
  stopifnot("xnames must be distinct from yname and uname" = base::intersect(c(yname, uname), xnames) == 0)
  stopifnot("yname, uname and xnames must match mae names" = all(c(yname, uname, xnames) %in% names(mae)))
  stopifnot("design must be a matrix" = is.matrix(design))
  stopifnot("design rownames must correspond to mae[[yname]] columns" = (! is.null(rownames(design))) && all(rownames(design) %in% colnames(mae[[yname]])))
  stopifnot("each sample in design can be assigned only to one group" = all(rowSums(design) == 1 | rowSums(design) == 0))
  stopifnot("at least one sample in design must be assigned to a group" = sum(rowSums(design)) > 0)
  stopifnot("standardize must be TRUE or FALSE" = isTRUEorFALSE(standardize))
  stopifnot("parallel must be TRUE or FALSE" = isTRUEorFALSE(parallel))
  stopifnot("pvalues must be TRUE or FALSE" = isTRUEorFALSE(pvalues))
  if (! is.null(precalcmodels)) {
    stopifnot("precalcmodels must be a list" = is.list(precalcmodels))
    stopifnot("precalcmodels elements names must be included in xnames" = (! is.null(names(precalcmodels))) && all(names(precalcmodels) %in% xnames))
    stopifnot("precalcmodels must contain only objects of class 'cv.glmnet'" = all(vapply(X = unlist(precalcmodels, recursive = FALSE), FUN = is, FUN.VALUE = logical(1L), class2 = "cv.glmnet")))
    stopifnot("precalcmodels is not compatible with mae" =
                all(vapply(
                  X = names(precalcmodels),
                  FUN = function(xnm) {
                    all(vapply(
                      X = precalcmodels[[xnm]],
                      FUN = function(x) {
                        identical(rownames(x[["glmnet.fit"]][["beta"]]), colnames(mae[[xnm]]))
                      },
                      FUN.VALUE = logical(1L)
                    ))
                  },
                  FUN.VALUE = logical(1L)
                )))
  }

  groups <- design2factor(design)

  iter_to_pass <- lapply(precalcmodels, names)

  args <- list(...)
  args[["offset"]] <- mae[[uname]]
  args[["standardize"]] <- standardize
  args[["parallel"]] <- parallel

  # extract X from MAE structure -- easier to work with
  X <- foreach::foreach(
    xnm = xnames,
    .inorder = TRUE,
    .final = function(x) setNames(x, xnames)
  ) %do% mae[[xnm]]

  message("##------ modelGeneExpression: started ridge regression ", timestamp(prefix = "", quiet = TRUE))
  regression_models <- foreach::foreach(
    x = iterators::iter(X),
    xn = xnames,
    .inorder = TRUE,
    .final = function(x) setNames(x, xnames),
    .packages = "xcore"
  ) %:%
    foreach::foreach(
      y = iterators::iter(mae[[yname]][, names(groups), drop = FALSE]),
      yn = names(groups),
      .inorder = TRUE,
      .final = function(x) setNames(x, names(groups)),
      .packages = "xcore"
    ) %dopar% {
      if (yn %in% iter_to_pass[[xn]]) {
        res <- "precalc"
      } else {
        args[["x"]] <- x
        args[["y"]] <- y
        res <- do.call(runLinearRidge, args) # TODO sometimes I am getting NULL here, not sure why though maybe I am running out of memory? should protect somehow against it; but then it is not like the most numerous context are getting those nulls...; eg. "ARPE-19_EMT_induced_with_TGF-beta_and_TNF-alpha" takes very long to complete even for small number of samples; it would appear that some particular samples take longer to run than others?
      }

      res
    }
  # add precalculated models
  for (xn in names(iter_to_pass)) {
    for (yn in iter_to_pass[[xn]])
      print(iter_to_pass[[xn]][[yn]])
    # regression_models[[xn]][[yn]] <- precalcmodels[[xn]][[yn]]
  }
  message("##------ modelGeneExpression: finished begining ridge  ", timestamp(prefix = "", quiet = TRUE))

  if (pvalues) {
    message("##------ modelGeneExpression: started significance testing  ", timestamp(prefix = "", quiet = TRUE))
    if (standardize) { # scale here to avoid scaling multiple times in ridgePvals
      X <- foreach::foreach(
        x = iterators::iter(X),
        .inorder = TRUE,
        .final = function(x) setNames(x, xnames)
      ) %dopar% scale(x)
    }
    svdX <- foreach::foreach(
      x = iterators::iter(X),
      .inorder = TRUE,
      .final = function(x) setNames(x, xnames)
    ) %dopar% svd(x)

    pvalues <- foreach::foreach(
      x = iterators::iter(X),
      xnm = xnames,
      .inorder = TRUE,
      .final = function(x) setNames(x, xnames)
      # .packages = "xcore" TODO set after development
    ) %:%
      foreach::foreach(
        y = iterators::iter(mae[[yname]][, names(groups), drop = FALSE]),
        id = names(groups),
        .inorder = TRUE,
        .final = function(x) setNames(x, names(groups))
        # .packages = "xcore"
      ) %dopar% {
        y <- y - mae[[uname]]
        lambda <- regression_models[[xnm]][[id]]$lambda.min
        beta <- coef(regression_models[[xnm]][[id]], s = lambda)
        beta <- beta[-1, ] # drop intercept TODO set after development
        ridgePvals(
          x = x,
          y = y,
          beta = beta,
          lambda = lambda,
          standardizex = FALSE,
          svdX = svdX[[xnm]]
        )
      }

    # TODO test if replicates are available
    replicate_avg <- lapply(pvalues, repAvgZscore, groups = groups)
    message("##------ modelGeneExpression: finished significance testing  ", timestamp(prefix = "", quiet = TRUE))
  } else {
    pvalues <- NULL
    replicate_avg = NULL
  }

  return(list(
    regression_models = regression_models,
    pvalues = pvalues,
    replicate_avg = replicate_avg
  ))
}

#' Calculate replicate averaged Z-scores
#'
#' Replicate averaged Z-scores for is calculated by dividing replicate average
#' coefficient by replicate pooled standard error.
#'
#' @param pvalues Data frame with \code{'se'} (standard error) and \code{'coef'}
#'   (coefficient) columns. Such as in \code{pvalues} output of
#'   \code{modelGeneExpression} .
#' @param groups Factor giving group membership for samples in  \code{pvalues}.
#'
#' @return Numeric matrix of averaged Z-scores. Columns correspond to groups and
#'   rows to predictors.
#'
repAvgZscore <- function(pvalues, groups) {
  stopifnot("pvalues elements must be instances of class data.frame" = all(vapply(pvalues, function(x) is(x, "data.frame"), logical(1L))))
  stopifnot("pavalues must have 'se' and 'coef' columns" = all(c("se", "coef") %in% colnames(pvalues)))
  stopifnot("groups must be a factor" = is.factor(groups))
  stopifnot("groups length must equal pvalues length" = length(groups) == length(pvalues))
  stopifnot("groups must not have unused levels" = setdiff(levels(groups), groups) == character(0))

  se <- applyOverDFList(list_of_df = pvalues, col_name = "se", fun = poolSE, groups = groups)
  estimate <- applyOverDFList(list_of_df = pvalues, col_name = "coef", fun = mean, groups = groups)
  zscore <- estimate / se

  return(zscore)
}

#' Pool Standard Error / Standard Deviation
#'
#' Pooled standard erro is calculated following (Cohen 1977) formulation for
#' pooled standard deviation.
#' TODO check out https://www.statisticshowto.com/find-pooled-sample-standard-error/, https://www.statisticshowto.com/pooled-standard-deviation/
#'
#' @param x Numeric vector of standard errors to pool.
#'
#' @return Number giving pooled standard error.
#'
poolSE <- function(x) {
  sqrt(sum(x * x) / length(x))
}
