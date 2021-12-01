#' xcore example expression data
#'
#' Expression data intended for use in xcore vignette and examples. It is build
#' from FANTOM5's 293SLAM rinderpest infection time course dataset. Here the
#' data is only a subset limited to core promoters (\code{promoters_f5_core}).
#'
#' @format A \code{matrix} with 14191 rows and 6 columns holding expression
#'   counts from CAGE-seq experiment. Rows corresponds to FANTOM5 promoters and
#'   columns to time points at which expression was measured 0 and 24 hours post
#'   infection.
#'
"rinderpest_mini"

#' xcore example molecular signatures
#'
#' Molecular signatures data intended for use in xcore vignette and examples. 
#' It is build ReMap2020 molecular signatures constructed against FANTOM5 
#' annotation, which can be found in \code{\link[xcoredata]{xcoredata}} package. 
#' Here the data is only a subset limited to core promoters 
#' (\code{promoters_f5_core}) and randomly selected 600 signatures.
#'
#' @format A \code{dgCMatrix} with 14191 rows and 600 columns holding 
#'   interaction matrix for subset of ReMap2020 molecular signatures against 
#'   FANTOM5 annotation. Rows corresponds to FANTOM5 promoters and columns to 
#'   signatures.
#'
"remap_mini"
