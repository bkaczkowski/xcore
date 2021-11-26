#!/usr/bin/env R
# create rinderpest_mini dataset to use as an example dataset in xcore
devtools::load_all()

data("promoters_f5_core", package = "xcoredata")
F5_counts_file <-
  system.file(
    "inst",
    "extdata",
    "hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz",
    package = "xcore"
  )
F5_counts <-
  read.table(
    file = gzfile(F5_counts_file),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
F5_counts <- F5_counts[-1, ]
rinderpest <- F5_counts[, grep("293SLAM.*rinderpest", colnames(F5_counts))]
colnm <- sub(".*(..hr).*(rep[0-9]).*", "\\1_\\2", colnames(rinderpest))
rinderpest <-
  matrix(
    data = unlist(rinderpest),
    ncol = 12,
    byrow = FALSE,
    dimnames = list(F5_counts[["00Annotation"]], colnm)
  )

# subset core promoters only
rinderpest_mini <- rinderpest[promoters_f5_core$name, ]

usethis::use_data(rinderpest_mini, overwrite = TRUE)
