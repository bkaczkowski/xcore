#!/usr/bin/env R
# ReMap2020 metada
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

remap_id <- c(colnames(remap_promoters), colnames(remap_enhancers)) %>%
  unique()

remap_meta <- remap_id %>%
  stringr::str_split(pattern = "\\.", n = 3) %>%
  do.call(what = rbind) %>%
  data.table::as.data.table()
colnames(remap_meta) <- c("id", "tf", "background")
remap_meta$biotype <- remap_meta$background %>% 
  gsub(pattern = "_.*", replacement = "")
remap_meta$condition <- remap_meta$background %>% 
  sub(pattern = "(.*?)_", replacement = "")
remap_meta$condition[! grepl("_", remap_meta$background)] <- ""
remap_meta[, background := NULL]

# fix study ids
remap_meta$study <- remap_meta$id
remap_meta$study[grepl(pattern = "ENCSR", x = remap_meta$id)] <- 
  translate(grep(pattern = "ENCSR", x = remap_meta$id, value = TRUE),
            srx2srastudy$Experiment,
            srx2srastudy$SRAStudy)
remap_meta$study <- ifelse(is.na(remap_meta$study), remap_meta$id, remap_meta$study)

# restore ids
remap_meta$id <- remap_id

data.table::setcolorder(remap_meta, c("id", "tf", "biotype", "study", "condition"))

usethis::use_data(remap_meta, overwrite = TRUE)
