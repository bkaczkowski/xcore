#!/usr/bin/env R
# create remap_mini interaction matrix to use as an example in xcore
devtools::load_all()

set.seed(123453)
data("remap_promoters", "promoters_f5_core", package = "xcoredata")

remap_mini <- remap_promoters[promoters_f5_core$name, ]
remap_mini <- remap_mini[, sample(1:ncol(remap_mini), 600)]

usethis::use_data(remap_mini, overwrite = TRUE)