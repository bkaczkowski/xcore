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
#' @param mat Numeric matrix of pairwise correlations to be plotted.
#' @param meta data.table metadata associated with \code{mat}. Should contain
#'   \code{id} column, which is used to subset \code{mat}. Other character
#'   columns are used as heatmap annotations.
#' @param ... arguments passed to \code{\link[pheatmap]{pheatmap}} function.
#'
#' @return Invisibly a pheatmap object. See \code{\link[pheatmap]{pheatmap}} for
#'   more informations.
#'
#' @importFrom pheatmap pheatmap
#' @importFrom proxy dist
#'
#' @export
plotExperimentsHeatmap <- function(mat,
                                   meta,
                                   cluster_rows = TRUE,
                                   cluster_cols = TRUE,
                                   angle_col = 45,
                                   treeheight_col = 0,
                                   fontsize = 8,
                                   fontsize_row = 7,
                                   fontsize_col = 7,
                                   ...) {
  stopifnot(is(mat, "matrix"))
  stopifnot(is(meta, "data.table"))
  stopifnot("id" %in% colnames(meta))
  stopifnot(sum(meta$id %in% colnames(mat)) >= 2)

  mat <- mat[meta$id, meta$id]

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
    mat = mat,
    scale = "none",
    annotation_row = annotaion,
    annotation_col = annotaion,
    annotation_colors = annotation_colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    angle_col = angle_col,
    treeheight_col = treeheight_col,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    ...
  )
}

#' Plot subtree
#'
#'
plotSubtree <- function(hc, meta, j) {
  tree <- as.dendrogram(hc)
  subtree <- dendextend::find_dendrogram(tree, meta[["id"]])
  if (is.null(subtree)) {stop("Could not find subtree. It happens sometimes.")}
  data.table::setkeyv(meta, "id") 
  cols <- meta[stats:::labels.dendrogram(subtree)][[j]]
  dendextend::labels_colors(subtree) <- paintVector(cols)
  par(cex=0.8, mar=c(15, 4, 1, 1))
  plot(subtree)
}
