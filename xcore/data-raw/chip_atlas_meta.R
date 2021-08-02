#!/usr/bin/env R
# ChIP-Atlas metada
devtools::load_all()

# Studies ids fix
srx2srastudy <- data.table::fread(
  file = system.file(
    "inst",
    "extdata",
    "srx2study.csv",
    package = "xcore"
  )
)

translate <- function(x, key, value) {
  newx <- value[match(x, key)]
}

chip_atlas_id <- c(colnames(chip_atlas_promoters), colnames(chip_atlas_enhancers)) %>%
  unique()

chip_atlas_meta <- chip_atlas_id %>%
  stringr::str_split(pattern = "_", n = 4) %>%
  do.call(what = rbind) %>%
  data.table::as.data.table()
colnames(chip_atlas_meta) <- c("tf", "background", "biotype", "id")

# fix study ids
chip_atlas_meta$study <- 
  srx2srastudy[match(chip_atlas_meta$id, srx2srastudy$Experiment), ]$SRAStudy
chip_atlas_meta$study <-
  ifelse(is.na(chip_atlas_meta$study),
         chip_atlas_meta$id,
         chip_atlas_meta$stud)

# restore ids
chip_atlas_meta$id <- chip_atlas_id

data.table::setcolorder(chip_atlas_meta, c("id", "tf", "biotype", "study", "background"))

usethis::use_data(chip_atlas_meta, overwrite = TRUE)
