test_that("simplifyInteractionMatrix", {
  testthat::expect_equal(
    simplifyInteractionMatrix(Matrix::Matrix(
      matrix(
        c(1, 1, 0, 0, 0, 1, 1, 1, 0),
        ncol = 3,
        byrow = T,
        dimnames = list(1:3, c("A", "A", "A"))
      ), sparse = TRUE),
      0.5,
      "AAA"),
    new(
      "dgCMatrix",
      i = c(0L, 2L),
      p = c(0L, 2L),
      Dim = c(3L, 1L),
      Dimnames = list(c("1", "2", "3"), "AAA"),
      x = c(1, 1),
      factors = list())
  )
})

test_that("design2factor", {
  design <- matrix(data = c(1, 1, 0, 0, 0, 0, 1, 1),
                   nrow = 4,
                   ncol = 2,
                   dimnames = list(c(paste("sample", 1:4)), c("gr1", "gr2")))
  testthat::expect_equal(
    design2factor(design), 
    factor(x = setNames(c("gr1", "gr1", "gr2", "gr2"), rownames(design)))
  )
})
