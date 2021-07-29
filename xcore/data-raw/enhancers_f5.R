# Bogumil Kaczkowski August 5, 2019
# The purpose of the script is to combine FANTOM5 hg19 and hg38 enhnacer
# and to create a non redundant set
# Additionally the enhancers are overlapped with ROADMAP enhancer, promoter and dyadic regios
# as well as with EP300 peaks.
devtools::load_all()

# HG38 ENHANCERS (from reprocessed/remapped data)
hg38_enh_file <- system.file("inst", "extdata", "F5.hg38.enhancers.bed.gz", package = "xcore")
hg38_enh <- rtracklayer::import.bed(hg38_enh_file)
hg38_enh$itemRgb <- NULL
hg38_enh$thick <- NULL
hg38_enh$blocks <- NULL
hg38_enh$score <- NULL

# # merging hg38 enhancers with lifted-over and adjusted hg19-only enhancer
# all_enh = c(hg38_enh ,hg19_enh_reduced )
# all_enh <- sortSeqlevels(all_enh)
# all_enh <- sort(all_enh)
all_enh <- hg38_enh
all_enh <- sortSeqlevels(all_enh)
all_enh <- sort(all_enh)

seqinfo_hg38 <- rtracklayer::SeqinfoForUCSCGenome("hg38")
seqinfo_hg38 <- seqinfo_hg38[
  names(seqinfo_hg38)[names(seqinfo_hg38) %in% seqnames(all_enh)]
  ]
seqinfo(all_enh) <- seqinfo_hg38

#.... hg38 enhancers...
all_enh$hg38_enhancer <- ""
hits <- findOverlaps(all_enh, hg38_enh)
all_enh$hg38_enhancer[hits@from] <- hg38_enh$name[hits@to]

#.... EP300 peaks...
all_enh$ep300 <- ""
ep300_file <-
  system.file("inst", "extdata", "remap2020_EP300_nr_macs2_hg38_v1_0.bed.gz", package = "xcore")
ep300 <- rtracklayer::import.bed(ep300_file)
hits <- findOverlaps(all_enh, ep300)
all_enh$ep300[hits@from] <- "EP300"
all_enh$ep300 <- as.factor(all_enh$ep300)

# couldn't find this dataset
# # ROADMAP
# roadmap_hg19 = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_bed6.bed", colnames = c("chrom", "start", "end", "name", "score"))
# roadmap = rtracklayer::liftOver(roadmap_hg19, chain = hg19_hg38_chain )
# roadmap = unlist(roadmap)
#
# all_enh$roadmap = ""
# hits = findOverlaps(all_enh , roadmap )
# all_enh$roadmap[hits@from] = roadmap$name[hits@to]
# all_enh$roadmap_type    =  as.factor( sub("\\|.*" , "" , all_enh$roadmap) )
# rm(roadmap_hg19, roadmap ); gc()

dfam_file <- system.file("inst", "extdata", "hg38_dfam.3.1.nrph.hits.bed.gz", package = "xcore")
dfam <- rtracklayer::import.bed(dfam_file)
all_enh$repeat_dfam <- ""
hits <- findOverlaps(all_enh, dfam)
all_enh$repeat_dfam[hits@from] <- dfam$name[hits@to]

enhancers_f5 <- all_enh

# tissue specificity tau
# Matrix of promoters normalized expression across datasets in FANTOM5,
# column names are not imported!
enhancers_f5_expression <-
  data.table::fread(
    file = system.file("inst",
                       "extdata",
                       "F5.hg38.enhancers.expression.tpm.matrix.gz",
                       package = "xcore"),
    header = TRUE)
enhancers_f5_expression <-
  data.table::melt(data = enhancers_f5_expression,
                   id.vars = "V1")
enhancers_f5_tau <-
  enhancers_f5_expression[,
                          .(t = tau(log(value + 1))),
                          by = "V1"
  ]

enhancers_f5$tau <- GenomicRanges::mcols(enhancers_f5) %>%
  as.data.frame() %>%
  dplyr::left_join(y = enhancers_f5_tau[, c("V1", "t")],
                   by = c("name" = "V1")) %>%
  dplyr::pull(t)

rm(enhancers_f5_expression, enhancers_f5_tau); gc()

# ENCODE blacklisted regions
blacklist <- rtracklayer::import(system.file("inst",
                                             "extdata",
                                             "hg38-blacklist.v2.bed.gz",
                                             package = "xcore")) %>%
  intersectGR(a = enhancers_f5,
              type = "any") %>%
  GenomicRanges::mcols() %>%
  `[[`("name")
enhancers_f5$encode_blacklist <- enhancers_f5$name %in% blacklist

rm(blacklist); gc()

# ENSEMBL soft masked regions
ensembl_sm <- rtracklayer::import(system.file("inst",
                                              "extdata",
                                              "hg38_ensembl_sm.bed.gz",
                                              package = "xcore")) %>%
  intersectGR(a = enhancers_f5,
              type = "any") %>%
  GenomicRanges::mcols() %>%
  `[[`("name")
enhancers_f5$ensembl_sm <- enhancers_f5$name %in% ensembl_sm

rm(ensembl_sm); gc()

usethis::use_data(enhancers_f5, internal = FALSE, overwrite = TRUE)

