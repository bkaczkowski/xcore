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

test_that("isTRUEorFALSE", {
  testthat::expect_equal(isTRUEorFALSE(TRUE), TRUE)
  testthat::expect_equal(isTRUEorFALSE(FALSE), TRUE)
  testthat::expect_equal(isTRUEorFALSE(NA), FALSE)
  testthat::expect_equal(isTRUEorFALSE(LETTERS), FALSE)
})

test_that("applyOverDFList", {
  # works on one element list
  list_of_df <- list(cars = cars)
  col_name <- "speed"
  fun <- mean
  groups <- factor("group1")
  names(groups) <- "cars"
  testthat::expect_equal(
    applyOverDFList(list_of_df, col_name, fun, groups),
    matrix(data = cars$speed, ncol = 1, dimnames = list(NULL, "group1"))
  )

  # works on multiple element list
  list_of_df <- list(a = cars, b = cars + 1, c = cars + 2)
  col_name <- "speed"
  fun <- mean
  groups <- factor(c("group1", "group1", "group2"))
  names(groups) <- c("a", "b", "c")
  testthat::expect_equal(
    applyOverDFList(list_of_df, col_name, fun, groups),
    matrix(data = c(cars$speed + 0.5, cars$speed + 2), ncol = 2, dimnames = list(NULL, c("group1", "group2")))
  )
})

test_that("translateCounts", {
  counts <- matrix(data = 1:8,
                   nrow = 4,
                   ncol = 2,
                   dimnames = list(c(paste0("gene", 1:4)), c("gr1", "gr2")))
  dict <- c(gene1 = "foo", gene2 = "bar", gene3 = "foo", gene4 = "zzz")

  testthat::expect_error(
    translateCounts("foo"),
    "counts must be an instance of class 'matrix'"
  )

  testthat::expect_error(
    translateCounts(matrix()),
    "counts must have its rownames defined"
  )

  bad_counts <- counts
  rownames(bad_counts)[2] <- rownames(bad_counts)[1]
  testthat::expect_error(
    translateCounts(bad_counts),
    "counts rownames must be unique"
  )

  testthat::expect_error(
    translateCounts(counts, 1:6),
    "dict must be an instance of class 'character'"
  )

  testthat::expect_error(
    translateCounts(counts, letters),
    "dict must be a named character vector"
  )

  bad_dict <- dict
  names(bad_dict)[2] <- names(bad_dict)[1]
  testthat::expect_error(
    translateCounts(counts, bad_dict),
    "dict names must be unique"
  )

  bad_dict <- dict
  names(bad_dict) <- letters[1:4]
  testthat::expect_error(
    translateCounts(counts, bad_dict),
    "counts rownames must correspond to dict names"
  )

  testthat::expect_equal(
    translateCounts(counts, dict),
    structure(c(2L, 4L, 4L, 6L, 12L, 8L),
              dim = 3:2,
              dimnames = list(c("bar", "foo", "zzz"), c("gr1", "gr2")))
  )
})
