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
#' @param lambda a string specifying if "lambda.min" or "lambda.1se" should be 
#'   used.
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
           standardize = TRUE,
           ...) {
  cv <- glmnet::cv.glmnet(
    x = x,
    y = y,
    offset = offset,
    alpha = 0,
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