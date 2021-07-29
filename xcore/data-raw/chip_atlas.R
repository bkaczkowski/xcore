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
                  "GFP")
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
promoters_ext500 <- promoters_f5
GenomicRanges::start(promoters_ext500) <- GenomicRanges::start(promoters_f5) - 500
GenomicRanges::end(promoters_ext500) <- GenomicRanges::end(promoters_f5) + 500

data("enhancers", package = "xcore")
enhancers_ext500 <- enhancers
GenomicRanges::start(enhancers_ext500) <- GenomicRanges::start(enhancers) - 500
GenomicRanges::end(enhancers_ext500) <- GenomicRanges::end(enhancers) + 500

rm(enhancers, promoters_f5)
gc()

# 3) create intersection matrises
overlap_mat_promoters <- matrix(data = 0L,
                                nrow = length(promoters_ext500),
                                ncol = length(unique_srxIDs))
colnames(overlap_mat_promoters) <- unique_srxIDs
rownames(overlap_mat_promoters) <- promoters_ext500$name

overlap_mat_enhancers <- matrix(data = 0L,
                                nrow = length(enhancers_ext500),
                                ncol = length(unique_srxIDs))
colnames(overlap_mat_enhancers) <- unique_srxIDs
rownames(overlap_mat_enhancers) <- enhancers_ext500$name

total_number_of_peaks <- rep(NA, length(unique_srxIDs))
names(total_number_of_peaks) <- unique_srxIDs
gc()

for (j in 1:length(unique_srxIDs)){
  query_tf <- chip_atlas[chip_atlas$name == unique_srxIDs[j]]
  total_number_of_peaks[j] <- length(query_tf)

  hits <- GenomicRanges::findOverlaps(promoters_ext500, query_tf)
  hits <- S4Vectors::DataFrame(from = hits@from,
                               to = hits@to,
                               score = query_tf$score[hits@to])
  hits <- hits[order(hits$score, decreasing = TRUE), ]
  hits <- hits[! duplicated(hits$from), ]
  overlap_mat_promoters[hits$from, j] <- query_tf$score[hits$to]

  hits <- GenomicRanges::findOverlaps(enhancers_ext500, query_tf)
  hits <- S4Vectors::DataFrame(from = hits@from,
                               to = hits@to,
                               score = query_tf$score[hits@to])
  hits <- hits[order(hits$score, decreasing = T), ]
  hits <- hits[! duplicated(hits$from), ]
  overlap_mat_enhancers[hits$from, j] <- query_tf$score[hits$to]

  print(j)
}

rm(hits, query_tf, j, chip_atlas) # , export_dir)
gc()

# 4) Meta data for the chipseq experiments
sxr_meta_file <-
  system.file("inst", "extdata", "experimentList_TF_hg38.txt", package = "xcore")
sxr_meta <- read.delim(file = sxr_meta_file,
		      sep ="\t",
		      header = FALSE,
		      quote = "\"")
dim(sxr_meta)
sxr_meta <- sxr_meta[sxr_meta$V1 %in% unique_srxIDs, ]
sxr_meta <- sxr_meta[, c(1, 4:7)]
colnames(sxr_meta) <- c("id", "TF", "origin", "cell", "description")
sxr_meta <- data.frame(apply(sxr_meta, 2, as.character), stringsAsFactors = FALSE)
sxr_meta <- sxr_meta[order(sxr_meta$id), ]
table(sxr_meta$origin)

sxr_meta$name <- paste(sxr_meta$TF, sxr_meta$origin, sxr_meta$cell, sxr_meta$id, sep = "_")
sxr_meta$name <- gsub(" " , ".", sxr_meta$name)

table(sxr_meta$id == colnames(overlap_mat_enhancers))
table(sxr_meta$id == colnames(overlap_mat_promoters))

colnames(overlap_mat_enhancers) <- sxr_meta$name
colnames(overlap_mat_promoters) <- sxr_meta$name

sxr_meta$number_of_promoters_with_peak <- colSums(overlap_mat_promoters > 0)
sxr_meta$number_of_enhancers_with_peak <- colSums(overlap_mat_enhancers > 0)

# 5) Gene level summarization
empty <- promoters_ext500$SYMBOL == ""
promoters_ext500_for_gene_collapsing <- promoters_ext500[! empty]
overlap_mat_promoters_for_gene_collapsing <- overlap_mat_promoters[! empty, ]

empty <- promoters_ext500$ENTREZID == ""
promoters_ext500_for_entrez_collapsing <- promoters_ext500[! empty]
overlap_mat_promoters_for_entrez_collapsing <- overlap_mat_promoters[! empty, ]

max_score_symbol <- function(x) {
  tapply(X = x,
	 INDEX = promoters_ext500_for_gene_collapsing$SYMBOL,
	 FUN = max)
}
max_score_entrez <- function(x) {
  tapply(X = x,
	 INDEX = promoters_ext500_for_entrez_collapsing$ENTREZID,
	 FUN = max)
}

overlap_mat_symbol <- apply(X = overlap_mat_promoters_for_gene_collapsing,
                            MARGIN = 2,
			    FUN = max_score_symbol)
overlap_mat_entrez <- apply(X = overlap_mat_promoters_for_entrez_collapsing,
                            MARGIN = 2,
			    FUN = max_score_entrez)

sxr_meta$number_of_geneSymbols_with_peak <- colSums(overlap_mat_symbol > 0)
sxr_meta$number_of_geneEntrez_with_peak <- colSums(overlap_mat_entrez > 0)
sxr_meta$total_number_of_peaks <- total_number_of_peaks

# Export data
chip_atlas_entrez <- Matrix::drop0(overlap_mat_entrez)
chip_atlas_symbol <- Matrix::drop0(overlap_mat_symbol)
# chip_atlas_meta <- sxr_meta

# promters to absence/presence matrix
chip_atlas_promoters <- Matrix::sparseMatrix(
  dims = c(nrow(overlap_mat_promoters), ncol(overlap_mat_promoters)),
  i = {},
  j = {})
colnames(chip_atlas_promoters) <- colnames(overlap_mat_promoters)
rownames(chip_atlas_promoters) <- rownames(overlap_mat_promoters)
for (j in 1:ncol(chip_atlas_promoters)) {
  chip_atlas_promoters[, j] <- as.logical(overlap_mat_promoters[, j])
}
chip_atlas_promoters <- as(chip_atlas_promoters, "dgCMatrix") # how remap is stored

# enhancers to absence/presence matrix
chip_atlas_enhancers <- Matrix::sparseMatrix(
  dims = c(nrow(overlap_mat_enhancers), ncol(overlap_mat_enhancers)),
  i = {},
  j = {})
colnames(chip_atlas_enhancers) <- colnames(overlap_mat_enhancers)
rownames(chip_atlas_enhancers) <- rownames(overlap_mat_enhancers)
for (j in 1:ncol(chip_atlas_enhancers)) {
  chip_atlas_enhancers[, j] <- as.logical(overlap_mat_enhancers[, j])
}
chip_atlas_enhancers <- as(chip_atlas_enhancers, "dgCMatrix") # how remap is stored

# usethis::use_data
usethis::use_data(chip_atlas_entrez, internal = FALSE, overwrite = TRUE)
usethis::use_data(chip_atlas_symbol, internal = FALSE, overwrite = TRUE)
# usethis::use_data(chip_atlas_meta, internal = FALSE, overwrite = TRUE)
usethis::use_data(chip_atlas_promoters, internal = FALSE, overwrite = TRUE)
usethis::use_data(chip_atlas_enhancers, internal = FALSE, overwrite = TRUE)
