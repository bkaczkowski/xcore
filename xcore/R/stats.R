#' Combine p-values using Fisher method
#'
#' @inheritParams stats::pchisq
#' @param p.value a numeric vector of p-values to combine.
#'
#' @return a number giving comined p-value.
#'
fisherMethod <- function(p.value, lower.tail = FALSE, log.p = TRUE) {
  stopifnot(is.numeric(p.value))
  stopifnot(length(p.value) > 1L)

  K <- 2 * length(p.value)
  X <- -2 * sum(log(p.value))
  stats::pchisq(X, df = K, lower.tail = lower.tail, log.p = log.p)
}

#' Fisher's test signatures
#'
#' Perform Fisher test between binary matrices columns. The function can be work
#' using parallelization  if  a parallel backend is registered (eg. \code{doMC}).
#'
#' @inheritParams stats::fisher.test
#' @param exprss_mat a binary expression matrix.
#' @param sign_mat a binary signature matrix.
#'
#' @return a list with \code{"p.value"} matrix of p-values and \code{"odds"}
#'   matrix of odds ratios from the Fisher test.
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom stats stats
#'
#' @export
fisherTestSignature <- function(exprss_mat, sign_mat, alternative = "greater") {
  stopifnot(is.matrix(exprss_mat)) # || is(exprss_mat, "dgCMatrix")) #TODO if we support sparseMatrices all methods eg. colSums must also support it
  stopifnot(is.matrix(sign_mat)) # || is(sign_mat, "dgCMatrix"))
  stopifnot(nrow(exprss_mat) == nrow(sign_mat))

  ecs <- colSums(exprss_mat)
  if (any(drop <- ecs == nrow(exprss_mat) | ecs == 0)) {
    warning(sprintf("Droping %i constant columns from expression", sum(drop)))
    exprss_mat <- exprss_mat[, ! drop]
  }
  esi <- colSums(sign_mat)
  if (any(drop <- esi == nrow(sign_mat) | esi == 0)) {
    warning(sprintf("Droping %i constant columns from signatures", sum(drop)))
    sign_mat <- sign_mat[, ! drop]
  }

  fishers <- foreach(j = seq_len(ncol(exprss_mat))) %dopar%
    apply(X = sign_mat,
          MARGIN = 2,
          FUN = function(sig) fisher.test(x = exprss_mat[, j, drop = TRUE],
                                          y = sig,
                                          alternative = alternative)
    )

  .getVal <- function(fishers, val, colnm) {
    mat <- do.call(what = cbind,
                   args = lapply(X = fishers,
                                 FUN = function(x) {
                                   do.call(what = rbind,
                                           args = lapply(x, `[[`, val))
                                 }))
    colnames(mat) <- colnm
    mat
  }

  list(
    p.value = .getVal(fishers, "p.value", colnames(exprss_mat)),
    odds = .getVal(fishers, "estimate", colnames(exprss_mat))
  )
}

#' Run linear ridge regression
#'
#' Linear ridge regression with lambda selection via CV.
#'
#' @inheritParams glmnet::glmnet
#'
#' @return An object of class "ridgeLinear".
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom ridge linearRidge
#'
#' @export
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

#' Calculate se and p-values for ridge regression betas
#'
ridgePValues <- function (x, y, beta, lambda, standardizex = TRUE, svdX = NULL) {
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

#'
#' @param mae MultiAssayExperiment object.
#' @param yname String response matrix name.
#' @param uname String response offset name.
#' @param xnames Character signatures names.
#' @param standardize Logical should signatures be scaled? scale(X).
#' @param parallel  parallel argument to glmnet::cv.glmnet.
#' @param ... arguments passed to glmnet::cv.glmnet.
#'
linearRidgePipeline <- function(mae,
                                yname,
                                uname,
                                xnames,
                                design,
                                standardize = TRUE,
                                parallel = TRUE,
                                simple_replicates_avg = TRUE,
                                elaborate_replicates_avg = TRUE,
                                precalcmodels = NULL,
                                ...) {
  stopifnot(is(mae, "MultiAssayExperiment"))
  stopifnot(base::intersect(yname, uname) == 0)
  stopifnot(base::intersect(yname, xnames) == 0)
  stopifnot(base::intersect(yname, uname) == 0)
  stopifnot(is.matrix(design))
  stopifnot(all(rownames(design) %in% colnames(mae[[yname]])))
  stopifnot(all(rowSums(design) == 1))
  if (! is.null(precalcmodels)) {
    stopifnot(all(names(precalcmodels) %in% xnames))
  }

  print("started regression")
  regression_models <- list()
  for (xnm in xnames) {
    regression_models[[xnm]] <- list()
    if (xnm %in% names(precalcmodels)) {
      regression_models[[xnm]] <- precalcmodels[[xnm]]
      next()
    }
    for (j in seq_len(ncol(mae[[yname]]))) {
      regression_models[[xnm]][[j]] <- runLinearRidge(
        x = mae[[xnm]],
        y = mae[[yname]][, j, drop = TRUE],
        offset = mae[[uname]],
        standardize = standardize,
        parallel = parallel,
        ...)
    }
    names(regression_models[[xnm]]) <- colnames(mae[[yname]])
  }

  # following edgeR::designAsFactor, but jokes aside
  # function (design)
  # {
  #   design <- as.matrix(design)
  #   z <- (exp(1) + pi)/5
  #   g <- factor(rowMeans(design * z^(col(design) - 1)))
  #   levels(g) <- seq_len(length(levels(g)))
  #   g
  # }
  groups <- factor(rowMeans(design * col(design) * ncol(design)))
  levels(groups) <- colnames(design)
  # can be tested like this table((names(groups) %>% sub(pattern = "_R[0-9][0-9]*$", replacement = "")) == groups)

  # simple replicates handling
  print("simple replicates handling")
  regression_simple_avg <- list()
  if (simple_replicates_avg) {
    for (xnm in xnames) {
      regression_simple_avg[[xnm]] <- foreach(gr = levels(groups)) %dopar%
          getAvgCoeff(regression_models[[xnm]][groups == gr])

      regression_simple_avg[[xnm]] <-
        do.call(what = cbind,
                args = lapply(regression_simple_avg[[xnm]], function (x) x[, "z"]))
      colnames(regression_simple_avg[[xnm]]) <- levels(groups)
    }
  }

  # elaborate replicates handling
  print("elaborate replicates handling")
  regression_elaborate_pvalues <- list()
  regression_elaborate_avg <- list()
  if (elaborate_replicates_avg) {
    for (xnm in xnames) {
      # calculate p-values (Cule E. 2011)
      X <- mae[[xnm]]
      if (standardize) X <- scale(X)
      svdX <- svd(X)
      regression_elaborate_pvalues[[xnm]] <-
        foreach(id = colnames(mae[[yname]])) %dopar%
        {
          y <- mae[[yname]][, id] - mae[[uname]]
          lambda <- regression_models[[xnm]][[id]]$lambda.min
          beta <- coef(regression_models[[xnm]][[id]], s = lambda)
          beta <- beta[-1, ] # drop intercept
          pval_test <- ridgePValues(
            x = X,
            y = y,
            beta = beta,
            lambda = lambda,
            standardizex = FALSE,
            svdX = svdX
          )
        }
      names(regression_elaborate_pvalues[[xnm]]) <- colnames(mae[[yname]])

      # calculate averaged z-scores
      regression_elaborate_avg[[xnm]] <- foreach(gr = levels(groups)) %dopar%
        {
          i <- groups == gr
          se <- lapply(regression_elaborate_pvalues[[xnm]][i], function (x) x$se)
          se <- do.call(cbind, se)
          se <- apply(se, 1, function (x) sqrt(sum(x * x) / length(x))) # Cohen, J. (1988) TODO check out https://www.statisticshowto.com/find-pooled-sample-standard-error/, https://www.statisticshowto.com/pooled-standard-deviation/
          estimate <- lapply(regression_elaborate_pvalues[[xnm]][i], function (x) x$coef)
          estimate <- do.call(cbind, estimate)
          estimate <- apply(estimate, 1, function (x) mean(x))
          estimate / se
        }
      regression_elaborate_avg[[xnm]] <-
        do.call(cbind, regression_elaborate_avg[[xnm]])
      rownames(regression_elaborate_avg[[xnm]]) <- colnames(X)
      colnames(regression_elaborate_avg[[xnm]]) <- levels(groups)
    }
  }

  return(list(
    regression_models = regression_models,
    regression_simple_avg = regression_simple_avg,
    regression_elaborate_pvalues = regression_elaborate_pvalues,
    regression_elaborate_avg = regression_elaborate_avg
  ))
}
