test_that("ridgePValues", {
  ridge_mod <- ridge::linearRidge(formula = mpg ~ ., 
                                  data = mtcars, 
                                  lambda = 2.747, 
                                  scaling = "scale")
  ridge_pvals <- ridge::pvals(ridge_mod)
  xcore_pvals <- ridgePValues(x = scale(mtcars[, -1]),
                              y = mtcars[, 1] - mean(mtcars[, 1]), # ridge subtract avg from y
                              beta = ridge_mod$coef[, 1], 
                              lambda = 2.747, standardizex = FALSE)
  testthat::expect_equal(
    round(as.numeric(ridge_pvals$se), 10),
    round(xcore_pvals$se, 10))
  testthat::expect_equal(
    round(as.numeric(ridge_pvals$tstat), 10),
    round(xcore_pvals$tstat, 10))
  testthat::expect_equal(
    round(as.numeric(ridge_pvals$pval), 10),
    round(xcore_pvals$pval, 10))
})
