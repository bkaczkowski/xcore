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
