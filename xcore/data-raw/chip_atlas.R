#!/usr/bin/env R
# Bogumil Kaczkowski
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions
# with ChipAtlas, Chip-Seq peak data-base
#
# Following steps are implemented
# 1) Chip-seq peaks
#   1a) Download ChipAtlas (all) peaks from http://chip-atlas.org (inst/extdata/download.sh)
#   1b) SHELL commands to process the downloaded BED file (inst/extdata/download.sh)
#   1c) read in the peak into R
# 2) get F5 promoter region
# 3) create intersection matrix
# 4) Fix the meta data for the chipseq experiments
devtools::load_all()

# 1c) read in the peak into R
chip_atlas_file <-
  system.file("inst", "extdata", "chip_atlas_hg38.Oth.ALL.05.AllAg.AllCell_SRX_only.bed.gz", package = "xcore")
chip_atlas <- rtracklayer::import(chip_atlas_file)
chip_atlas$score <- as.integer(chip_atlas$score)

# remove ambiguous records
ambiguous_tf <- c("5-hmC", "5-mC", "8-Hydroxydeoxyguanosine", "Biotin", "BMI",
                  "BrdU", "Cas9", "Cyclobutane", "MethylCap", "O-GlcNAc",
                  "Pan-acetyllysine", "pFM2", "RTA", "SVS-1", "EBNA1", "EBNA2",
                  "EBNA3", "EBV-ZEBRA", "AML1-ETO", "VSV-G", "MLL-AF4",
                  "MLL-AF6", "Hepatitis", "HIV", "KSHV", "MCPV", "Epitope tags",
                  "GFP", "HIV Tat", "MCPV ST", "KSHV LANA", 
                  "Cyclobutane pyrimidine dimers", "Hepatitis B Virus X antigen")
ambiguous_srx <-
  read.delim(
    file = system.file("inst", "extdata", "experimentList_TF_hg38.txt", package = "xcore"),
    sep = "\t",
    header = FALSE,
    quote = "\"") %>%
  dplyr::filter(V4 %in% ambiguous_tf) %>% # V4 == TF column
  dplyr::pull(V1) # V1 == experiment ID
chip_atlas <- chip_atlas[! chip_atlas$name %in% ambiguous_srx, ]

unique_srxIDs <- unique(chip_atlas$name)
unique_srxIDs <- unique_srxIDs[order(unique_srxIDs)]

# 2) get F5 promoter region
data("promoters_f5", package = "xcore")

# 3) create intersection matrises
chip_atlas_promoters <-
  xcore::getInteractionMatrix(a = promoters_f5,
                              b = chip_atlas,
                              ext = 500,
                              count = FALSE)
chip_atlas_promoters <- chip_atlas_promoters[, unique_srxIDs]

# 4) Meta data for the chipseq experiments
sxr_meta_file <-
  system.file("inst", "extdata", "experimentList_TF_hg38.txt", package = "xcore")
sxr_meta <- read.delim(file = sxr_meta_file,
		      sep ="\t",
		      header = FALSE,
		      quote = "\"")
sxr_meta <- sxr_meta[sxr_meta$V1 %in% unique_srxIDs, ]
sxr_meta <- sxr_meta[, c(1, 4:7)]
colnames(sxr_meta) <- c("id", "TF", "origin", "cell", "description")
sxr_meta <- data.frame(apply(sxr_meta, 2, as.character), stringsAsFactors = FALSE)
sxr_meta <- sxr_meta[order(sxr_meta$id), ]
sxr_meta$cell <- gsub(" ", "_", sxr_meta$cell)
sxr_meta$cell <- gsub("\\.", "_", sxr_meta$cell) # could not find any other way
sxr_meta$name <- paste(sxr_meta$TF, sxr_meta$cell, sxr_meta$id, sep = ".")

colnames(chip_atlas_promoters) <- sxr_meta$name

sxr_meta$number_of_promoters_with_peak <- colSums(chip_atlas_promoters > 0)

# 5) save
usethis::use_data(chip_atlas_promoters, internal = FALSE, overwrite = TRUE)
