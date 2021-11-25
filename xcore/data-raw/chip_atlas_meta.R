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

chip_atlas_id <- colnames(chip_atlas_promoters) %>% unique()

chip_atlas_meta <- chip_atlas_id %>%
  stringr::str_split(pattern = "\\.", n = 3) %>%
  do.call(what = rbind) %>%
  data.table::as.data.table()
colnames(chip_atlas_meta) <- c("tf", "biotype", "id")

# fix study ids
chip_atlas_meta$study <-
  srx2srastudy[match(chip_atlas_meta$id, srx2srastudy$Experiment), ]$SRAStudy
chip_atlas_meta$study <-
  ifelse(is.na(chip_atlas_meta$study),
         chip_atlas_meta$id,
         chip_atlas_meta$study)

# restore ids
chip_atlas_meta$id <- chip_atlas_id

# CIS-BP TF classification
cis_bp <-
  data.table::fread(system.file("inst", "extdata", "cis_bp_tf_class.txt",
                                package = "xcore"))
cis_bp <-
  cis_bp[, .(tf_dbd = unique(DBDs)),
         by = TF_Name]
data.table::setnames(cis_bp, "TF_Name", "tf")
chip_atlas_meta <- cis_bp[chip_atlas_meta, on = c("tf" = "tf")]

data.table::setcolorder(chip_atlas_meta, c("id", "tf", "tf_dbd", "biotype", "study"))
usethis::use_data(chip_atlas_meta, overwrite = TRUE)
