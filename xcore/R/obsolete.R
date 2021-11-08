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