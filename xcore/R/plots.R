#' Object Summary as a string
#'
#' Return object summary as a string that can be placed on a plot.
#'
#' @param ... arguments passed to \code{summary}.
#'
#' @return String.
#'
summaryString <- function(...) {
  s <- summary(...)
  snm <- names(s)
  sval <- unclass(s) %>% round(digits = 2L) %>% as.character()
  width <- c(snm, sval) %>% nchar() %>% max()
  snm <- format(x = snm, width = width, justify = "right") %>%
    paste(collapse = " ")
  sval <- format(x = sval, width = width, justify = "right") %>%
    paste(collapse = " ")
  paste0(snm, "\n", sval)
}

#' Histogram with summary
#'
#' \code{plotHistogram} create a histogram with summary on the upper margin.
#' Refer to \code{\link[graphics]{hist}} documentation for more informations.
#'
#' @inheritParams graphics::hist
#'
#' @return an object of class "histogram".
#'
plotHistogram <- function(x, breaks = 40, cex = 0.7, xlab = "", ...) {
  par(mar = c(5, 5, 6, 2))
  graphics::hist(x = x,
                 breaks = breaks,
                 cex = 0.7,
                 xlab = xlab,
                 ...)
  mtext(text = summaryString(x),
        side = 3,
        cex=0.6,
        family = "mono")
}

#' Experiment's heatmap
#'
#' \code{plotExperimentsHeatmap} crate a heatmap for a group o experiments
#' (columns) based on interaction matrix. As function does not use sparsity
#' aware methods for distance calculation it is advised to not use too large
#' matrices. \code{\link[colorspace]{qualitative_hcl}} color palette is used
#' to allow higer number of categories in the annotation. No scaling is applied
#' to \code{mat}!
#'
#' @inheritParams pheatmap::pheatmap
#' @param mat dgCMatrix interaction matrix such as produced by
#'   \code{\link{getInteractionMatrix}}.
#' @param meta data.table metadata associated with \code{mat}. Should contain
#'   \code{id} column, which is used to subset \code{mat}. Other character
#'   columns are used as heatmap annotations.
#' @param method a character string indicating which method should be used to
#'   calculate distances between \code{mat} columns. Supported values: "pearson",
#'   "kendall", "spearman" (the actual distance is calculated as
#'   (1 - cor) /2), "binary".
#' @param ... arguments passed to \code{\link[pheatmap]{pheatmap}} function.
#'
#' @return Invisibly a pheatmap object. See \code{\link[pheatmap]{pheatmap}} for
#'   more informations.
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
plotExperimentsHeatmap <- function(mat,
                                   meta,
                                   method = "pearson",
                                   cluster_rows = TRUE,
                                   cluster_cols = TRUE,
                                   angle_col = 45,
                                   treeheight_col = 0,
                                   ...) {
  stopifnot(is(mat, "dgCMatrix"))
  stopifnot(is(meta, "data.table"))
  stopifnot("id" %in% colnames(meta))
  stopifnot(sum(meta$id %in% colnames(mat)) >= 2)
  stopifnot(method %in% c("pearson", "kendall", "spearman", "binary"))

  mat <- mat[, meta$id]

  # experiments distance
  if (method == "binary") {
    tf_cor_d <- dist(Matrix::t(mat), method = "binary") # binary
    tf_cor <- as.matrix(tf_cor_d)
  } else {
    tf_cor <- cor(mat %>% as.matrix(), method = method)
    tf_cor_d <- as.dist(abs((1 - tf_cor) / 2))
  }

  # heatmap
  cols <- vapply(
    X = meta,
    FUN = function(x)
      is.character(x) || is.factor(x),
    FUN.VALUE = logical(1L)
  ) &
    (!colnames(meta) %in% "id")
  annotaion <- as.data.frame(meta[, cols])
  rownames(annotaion) <- meta$id
  annotaion <- rapply(annotaion, as.factor, how = "replace")
  annotation_colors <- lapply(
    X = annotaion,
    FUN = function(x) {
      stats::setNames(colorspace::qualitative_hcl(nlevels(x)), levels(x))
    }
  )

  pheatmap::pheatmap(
    mat = tf_cor,
    scale = "none",
    clustering_distance_rows = tf_cor_d,
    clustering_distance_cols = tf_cor_d,
    annotation_row = annotaion,
    annotation_col = annotaion,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    angle_col = angle_col,
    treeheight_col = treeheight_col,
    ...
  )
}
