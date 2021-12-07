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
  annotaion <- as.data.frame(meta[, cols, with = FALSE])
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

#' Customized color pallette
#' 
#'  
rdBuPalette <- function(n) {
  colorspace::diverge_hsv(n = n, power = 0.75, h = c(245, 0), s = 0.9)
}

#' Scale bar to use with heatmaps
#' 
#' 
plotHeatmapScaleBar <-
  function(breaks,
           col = rdBuPalette(length(breaks) - 1),
           side = 2,
           nticks = 3,
           mar = c(2.1, 3.1, 4.1, 1.1),
           cex = 0.75,
           cex.main = 0.9,
           ...) {
    omar <- par()$mar
    on.exit(par(mar = omar))
    par(mar = mar, cex = cex, cex.main = cex.main)
    
    image(
      matrix(data = breaks, nrow = 1),
      breaks = breaks,
      col = col,
      axes = FALSE,
      ...
    )
    axis(
      side = 2,
      at = seq(
        from = 0,
        to = 1,
        length.out = nticks
      ),
      labels = seq(
        from = min(breaks),
        to = max(breaks),
        length.out = nticks
      )
    )
    
    invisible()
  }

#' Heatmap
#' 
#' 
plotHeatmap <-
  function(x,
           breaks,
           col = rdBuPalette(length(breaks) - 1),
           rownames = TRUE,
           colnames = TRUE,
           clust_rows = TRUE,
           na.color = "grey",
           mar = c(6.1, 0.1, 4.1, 22.1),
           cex = 1,
           cex.main = 0.9,
           cex.axis = 0.8,
           ...) {
    omar <- par()$mar
    on.exit(par(mar = omar))
    graphics::par(mar = mar, cex = cex, cex.main = cex.main, cex.axis = cex.axis)

    # order rows
    if (clust_rows) {
      # tmp substitute NA with 0
      i_na <- is.na(x)
      x[i_na] <- 0
      
      if (nrow(x) > 1) {
        hc <- hclust(dist(x))
        ord <- order.dendrogram(as.dendrogram(hc))
        x <- x[rev(ord), , drop = FALSE] # ordering in image is reversed
      }
      
      # restore NA
      x[i_na[rev(ord)]] <- NA
    }
    
    # cap data to fit breaks
    x <- ifelse(test = x < min(breaks, na.rm = TRUE), 
                yes = min(breaks, na.rm = TRUE), 
                no = x)
    x <- ifelse(test = x > max(breaks, na.rm = TRUE), 
                yes = max(breaks, na.rm = TRUE), 
                no = x)
    
    graphics::image(
      x = t(x),
      breaks = breaks,
      col = col,
      axes = FALSE,
      ...
    )
    
    # plot NAs if needed -- taken from gplots::heatmap.2
    if (any(is.na(x))) {
      namat <- ifelse(is.na(x), 1, NA)
      image(
        x = t(namat),
        col = na.color,
        axes = FALSE,
        xlab = "",
        ylab = "",
        add = TRUE
      )
    }
    
    if (colnames) {
      graphics::axis(
        side = 1,
        at = seq(
          from = 0,
          to = 1,
          length.out = ncol(x)
        ),
        labels = colnames(x),
        las = 2,
        line = -0.5, 
        tick = 0,
        family = "mono"
      )
    }
    
    if (rownames) {
      graphics::axis(
        side = 4,
        at = seq(
          from = 0,
          to = 1,
          length.out = nrow(x)
        ),
        labels = rownames(x),
        las = 2,
        line = -0.5, 
        tick = 0,
        family = "mono"
      )
    }
    
    invisible()
  }

#' plot title
#' 
#' 
plotTitle <- function(labels, mar = c(0.1, 1.1, 0.1, 1.1)) {
  omar <- par()$mar
  on.exit(par(mar = omar))
  par(mar = mar, oma = c(0, 0, 0, 0))
  
  plot.new()
  graphics::text(
    x = 0.5,
    y = 0.5,
    adj = c(0.5, 0.5),
    labels = labels,
    cex = 1.75,
    font = 2
  )
  
  invisible()
}
