test_that("prepareCountsForRegression", {
  data("rinderpest_mini")
  design <- matrix(
    data = c(1, 0,
             0, 1),
    ncol = 2,
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("00hr_rep1", "24hr_rep3"), c("00hr", "24hr")))

  testthat::expect_error(
    prepareCountsForRegression(counts = LETTERS),
    "counts must be an integer matrix"
  )
  testthat::expect_error(
    prepareCountsForRegression(counts = rinderpest_mini, design = LETTERS),
    "design must be a matrix"
  )
  testthat::expect_error(
    prepareCountsForRegression(counts = rinderpest_mini, design = design),
    "number of rows in design must equal to number of columns in counts"
  )
  rinderpest_mini <- rinderpest_mini[, rownames(design)]
  fake_design <- design
  rownames(fake_design) <- paste0(rownames(fake_design), "FOO")
  testthat::expect_error(
    prepareCountsForRegression(counts = rinderpest_mini, design = fake_design),
    "design rownames must be the same as counts colnames"
  )
  testthat::expect_error(
    prepareCountsForRegression(
      counts = rinderpest_mini,
      design = design,
      base_lvl = "foo"),
    "base_lvl must match to one of design colnames"
  )
  testthat::expect_error(
    prepareCountsForRegression(
      counts = rinderpest_mini,
      design = design,
      base_lvl = "00hr",
      log2 = "foo"),
    "log2 must be TRUE or FALSE"
  )
  testthat::expect_error(
    prepareCountsForRegression(
      counts = rinderpest_mini,
      design = design,
      base_lvl = "00hr",
      pseudo_count = -1),
    "pseudo_count must be an positive integer or zero"
  )
  testthat::expect_error(
    prepareCountsForRegression(
      counts = rinderpest_mini,
      design = design,
      base_lvl = "00hr",
      drop_base_lvl = "foo"),
    "drop_base_lvl must be TRUE or FALSE"
  )

  mae <- prepareCountsForRegression(
    counts = rinderpest_mini,
    design = design,
    base_lvl = "00hr")
  testthat::expect_equal(digest::digest(mae), "a4fff024e3fa50599bf59de21684234a")
})

test_that("addSignatures", {
  data("rinderpest_mini", "remap_mini")
  base_lvl <- "00hr"
  design <- matrix(
    data = c(1, 0,
             1, 0,
             1, 0,
             0, 1,
             0, 1,
             0, 1),
    ncol = 2,
    nrow = 6,
    byrow = TRUE,
    dimnames = list(c("00hr_rep1", "00hr_rep2", "00hr_rep3", "24hr_rep1", "24hr_rep2", "24hr_rep3"), c("00hr", "24hr"))
  )
  mae <- prepareCountsForRegression(
    counts = rinderpest_mini[, c("00hr_rep1", "00hr_rep2", "00hr_rep3", "24hr_rep1", "24hr_rep2", "24hr_rep3")],
    design = design,
    base_lvl = base_lvl)

  testthat::expect_error(
    addSignatures(mae = LETTERS),
    "mae must be an instance of class 'MultiAssayExperiment'"
  )
  testthat::expect_error(
    addSignatures(mae = mae, intersect_rows = "foo"),
    "intersect_rows must be TRUE or FALSE"
  )
  testthat::expect_error(
    addSignatures(mae = mae, remap_mini),
    "experiments must be named"
  )
  testthat::expect_error(
    addSignatures(mae = mae, Y = remap_mini),
    "experiments names must be unique"
  )

  mae <- addSignatures(mae, remap = remap_mini)
  testthat::expect_equal(digest::digest(mae), "b0287dfe29a52559bab00b0ed8ff8a56")
})

test_that("getInteractionMatrix", {
  a <- GenomicRanges::GRanges(
    seqnames = c("chr20", "chr4"),
    ranges = IRanges::IRanges(start = c(62475984L, 173530220L), end = c(62476001L, 173530236L)),
    strand = c("-", "-"),
    name = c("hg19::chr20:61051039..61051057,-;hg_188273.1", "hg19::chr4:174451370..174451387,-;hg_54881.1")
  )
  b <- GenomicRanges::GRanges(
    seqnames = c("chr4", "chr20"),
    ranges = IRanges::IRanges(start = c(173530229L, 63864270L), end = c(173530236L, 63864273L)),
    strand = c("-", "-"),
    name = c("HAND2", "GATA5")
  )

  testthat::expect_error(
    getInteractionMatrix(a = 1L),
    "a must be an instance of class 'GRanges'"
  )
  testthat::expect_error(
    getInteractionMatrix(a = GenomicRanges::GRanges()),
    "a must have 'name' attribute in it's metadata"
  )
  testthat::expect_error(
    getInteractionMatrix(a = a, b = 1L),
    "b must be an instance of class 'GRanges'"
  )
  testthat::expect_error(
    getInteractionMatrix(a = a, b = GenomicRanges::GRanges()),
    "b must have 'name' attribute in it's metadata"
  )
  testthat::expect_error(
    getInteractionMatrix(a = a, b = b, ext = letters),
    "ext must be numeric larger or equal to 0"
  )
  testthat::expect_error(
    getInteractionMatrix(a = a, b = b, ext = -1),
    "ext must be numeric larger or equal to 0"
  )
  testthat::expect_error(
    getInteractionMatrix(a = a, b = b, ext = 0, count = NA),
    "count must be TRUE or FALSE"
  )

  testthat::expect_equal(
    getInteractionMatrix(a = a, b = b, ext = 0, count = FALSE),
    as(matrix(
      data = c(0, 1, 0, 0),
      nrow = 2,
      ncol = 2,
      dimnames = list(
        c(
          "hg19::chr20:61051039..61051057,-;hg_188273.1",
          "hg19::chr4:174451370..174451387,-;hg_54881.1"
        ),
        c("HAND2", "GATA5")
      )
    ), "dgCMatrix")
  )

  testthat::expect_equal(
    getInteractionMatrix(a = a, b = b, ext = 1500000, count = FALSE),
    as(matrix(
      data = c(0, 1, 1, 0),
      nrow = 2,
      ncol = 2,
      dimnames = list(
        c(
          "hg19::chr20:61051039..61051057,-;hg_188273.1",
          "hg19::chr4:174451370..174451387,-;hg_54881.1"
        ),
        c("HAND2", "GATA5")
      )
    ), "dgCMatrix")
  )
})

test_that("filterSignatures", {
  data("rinderpest_mini", "remap_mini")
  base_lvl <- "00hr"
  design <- matrix(
    data = c(1, 0,
             1, 0,
             1, 0,
             0, 1,
             0, 1,
             0, 1),
    ncol = 2,
    nrow = 6,
    byrow = TRUE,
    dimnames = list(c("00hr_rep1", "00hr_rep2", "00hr_rep3", "24hr_rep1", "24hr_rep2", "24hr_rep3"), c("00hr", "24hr"))
  )
  mae <- prepareCountsForRegression(
    counts = rinderpest_mini[, c("00hr_rep1", "00hr_rep2", "00hr_rep3", "24hr_rep1", "24hr_rep2", "24hr_rep3")],
    design = design,
    base_lvl = base_lvl)
  mae <- addSignatures(mae, remap = remap_mini)

  testthat::expect_error(
    filterSignatures(
      mae = 1L
    ),
    "mae must be an instance of class 'MultiAssayExperiment'"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      min = "foo"
    ),
    "min must be a length one numeric between 0 and 1"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      max = "foo"
    ),
    "max must be a length one numeric between 0 and 1"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      ref_experiment = 1
    ),
    "ref_experiment must be a length one character"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      ref_experiment = "foo"
    ),
    "ref_experiment must match one of mae names"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      omit_experiments = 1
    ),
    "omit_experiments must be a character"
  )
  testthat::expect_error(
    filterSignatures(
      mae = mae,
      omit_experiments = "foo"
    ),
    "omit_experiments must match mae names"
  )

  to_drop <-
    c(
      "GSE73492.EWSR1.HEK293",
      "GSE23436.GATA6.CACO2_DIFF",
      "GSE53809.DPPA3.HEK293T",
      "GSE54592.ESR1.MCF-7_LETR",
      "GSE47987.AR.DU145_ARC562S",
      "GSE90454.FOXA2.BJ1-hTERT_CDT1_Mimo",
      "GSE37345.AR.LNCaP_SHCTR_R1881_HD",
      "GSE63209.NR3C1.MM1-S_DEX"
    )
  test_mae <- mae
  suppressMessages(
    test_mae[["remap"]] <- test_mae[["remap"]][, ! colnames(test_mae[["remap"]]) %in% to_drop]
  )
  testthat::expect_equal(filterSignatures(mae, min = 0.0001, max = 1), test_mae)
})
