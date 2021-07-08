# A script to Annotate the FANTOM5 promoters/TSS regions
# getting BED files from FANTOM5 website one for the coordinates (bed) and one with official annotation
devtools::load_all()

dpi_bed_file  <-
  system.file("inst", "extdata", "hg38_fair+new_CAGE_peaks_phase1and2.bed.gz", package = "xcore")
dpi_annot_file <- 
  system.file("inst", "extdata", "hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz", package = "xcore")

dpi <- rtracklayer::import.bed(dpi_bed_file)
dpi$itemRgb <- NULL
dpi$thick <- NULL
dpi <- dpi[order(dpi$name), ]

seqinfo_hg38 <- rtracklayer::SeqinfoForUCSCGenome("hg38")
seqinfo_hg38 <- seqinfo_hg38[
  names(seqinfo_hg38) [names(seqinfo_hg38) %in% seqnames(dpi)]
  ]
seqinfo_hg38 <- seqinfo_hg38[names(seqinfo_hg38)[c(1:22, 25, 23, 24)]]

seqinfo(dpi) <- seqinfo_hg38

dpi_annot <- data.table::fread(dpi_annot_file, header = T, sep = "\t")
data.table::setDF(dpi_annot)
dpi_annot <- dpi_annot[dpi_annot$`00Annotation` %in% dpi$name, ]
dpi_annot <- dpi_annot[order(dpi_annot$`00Annotation`), ]

table(dpi_annot$`00Annotation` == dpi$name)
# TRUE
# 209911

# This is how example dpi_annot line looks like
# 00Annotation    short_description       description     association_with_transcript     entrezgene_id   hgnc_id uniprot_id
# hg19::chr1:564571..564600,+;hg_1.1      p1@MTND1P23     CAGE_peak_1_at_MTND1P23_5end    -130bp_to_ENST00000416931.1_5end        NA      HGNC:42092      NA
# hg19::chr1:840185..840202,+;hg_20.1     hg_20.1 CAGE_peak_hg_20.1       12bp_to_ENST00000607769.1,uc057axn.1_5end       NA      NA      NA
dpi$SYMBOL_F5_annot <- gsub("\\.*", "", gsub(".*@", "", dpi_annot$short_description)) # sub(" .*", "", dpi_annot$Gene_symbol)
# in $short_description we have a lot of hg_1234... which are not symbols acctually
dpi$SYMBOL_F5_annot <- gsub("hg_[0-9]+", "", dpi$SYMBOL_F5_annot)
dpi$ENTREZID_F5_annot <- sub(" .*", "", dpi_annot$entrezgene_id)
dpi$ENTREZID_F5_annot[is.na(dpi$ENTREZID_F5_annot)] <- ""
dpi$distance_F5_annot <- as.integer(gsub("bp_to_.*", "", dpi_annot$association_with_transcript)) # dpi_annot$Distance
dpi$distance_F5_annot[is.na(dpi$distance_F5_annot)] <- ""

# GENCODE 38 annotation
gencode_file <- system.file("inst", "extdata", "gencode.v38.annotation.gff3.gz", package = "xcore")
gencode <- rtracklayer::import.gff(con = gencode_file)
dpi_gencode <- xcore::gencode_nearest_promoter_same_strand(regions = dpi,
	        					   gencode = gencode,
							   cut_off_distance = 500)
dpi$SYMBOL_gencode <- dpi_gencode$nearest_symbol
select_first <- function(x) { x[1] }
dpi$ENTREZID_gencode <- as.character(lapply(
  X = AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = dpi$SYMBOL_gencode,
    column = 'ENTREZID',
    keytype = 'SYMBOL'),
  FUN = select_first))
dpi$ENTREZID_gencode[dpi$ENTREZID_gencode == "NULL"] <- ""
dpi$ENTREZID_gencode[is.na(dpi$ENTREZID_gencode)] <- ""
dpi$gene_type_gencode <- as.factor(dpi_gencode$nearest_gene_type)
dpi$distance_gencode <- dpi_gencode$nearest_distance

# UCSC.hg38.knownGene annotation
ucsc_knownGene <-
  GenomicFeatures::asGFF(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
ucsc_knownGene <- ucsc_knownGene[ucsc_knownGene$type == "mRNA"]
ucsc_knownGene_promoters <- GenomicRanges::promoters(ucsc_knownGene, upstream = 0, downstream = 0, use.names=TRUE)
ucsc_knownGene_promoters$Parent <- drop(ucsc_knownGene_promoters$Parent)
ucsc_knownGene_promoters <- ucsc_knownGene_promoters [! is.na(ucsc_knownGene_promoters$Parent)]
ucsc_knownGene_promoters$Parent <- sub("GeneID:", "", ucsc_knownGene_promoters$Parent)

dpi$SYMBOL_ucsc <- ""
dpi$ENTREZID_ucsc <- ""
dpi$distance_ucsc <- 10^6

hits <- GenomicRanges::distanceToNearest(dpi, ucsc_knownGene_promoters, ignore.strand=FALSE)
dpi$ENTREZID_ucsc[hits@from] <- ucsc_knownGene_promoters$Parent[hits@to]
dpi$distance_ucsc[hits@from] <- hits@elementMetadata$distance
select_first <- function(x) { x[1] }
dpi$SYMBOL_ucsc  = as.character(lapply(
  X = AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = dpi$ENTREZID_ucsc,
    column = 'SYMBOL',
    keytype = 'ENTREZID'),
  FUN = select_first))
dpi$SYMBOL_ucsc[dpi$SYMBOL_ucsc == "NULL"] <- ""
dpi$SYMBOL_ucsc[is.na(dpi$SYMBOL_ucsc)] <- ""

too_far <- dpi$distance_ucsc > 500
dpi$ENTREZID_ucsc[too_far] <- ""
dpi$SYMBOL_ucsc[too_far] <- ""

# simplifying/collapsing ENTREZID
dpi$ENTREZID <- dpi$ENTREZID_ucsc
empty <- dpi$ENTREZID == ""
dpi$ENTREZID[empty] <- dpi$ENTREZID_gencode[empty]
empty <- dpi$ENTREZID == ""
dpi$ENTREZID[empty] <- dpi$ENTREZID_F5_annot[empty]

# simplifying/collapsing SYMBOLs
dpi$SYMBOL <- dpi$SYMBOL_ucsc
empty <- dpi$SYMBOL == ""
dpi$SYMBOL[empty] <- dpi$SYMBOL_gencode[empty]
empty <- dpi$SYMBOL == ""
dpi$SYMBOL[empty] <- dpi$SYMBOL_F5_annot[empty]

# some check, stats and overview
s <- dpi$SYMBOL_F5_annot != ""
table(dpi$SYMBOL_F5_annot[s] == dpi$SYMBOL[s]) # the FALSE should be about 1.8%, not 15%!!!!
# FALSE  TRUE
#  3305 96576
# FALSE = 3.3089%
s <- dpi$ENTREZID_F5_annot != ""
table(dpi$ENTREZID_F5_annot[s] == dpi$ENTREZID[s]) # the FALSE should be below 1%
# FALSE  TRUE
#   969 96786
# FALSE = 0.99125%

promoters_f5 <- dpi
promoters_f5$distance_F5_annot <- NULL
promoters_f5$distance_gencode <- NULL
promoters_f5$distance_ucsc <- NULL

usethis::use_data(promoters_f5 ,internal = FALSE, overwrite = TRUE)

# # ROADMAP
# roadmap = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_hg38.bed", colnames = c("chrom", "start", "end", "name", "score"))
# dpi$roadmap = ""
# hits = findOverlaps(dpi , roadmap )
# dpi$roadmap[hits@from] = roadmap$name[hits@to]
# rm(roadmap_hg19, roadmap ); gc()

#EP300
dpi$ep300 <- ""
ep300_file <-
  system.file("inst", "extdata", "remap2020_EP300_nr_macs2_hg38_v1_0.bed.gz", package = "xcore")
ep300 <- rtracklayer::import.bed(ep300_file)
hits <- findOverlaps(dpi, ep300)
dpi$ep300[hits@from] <- "EP300"
dpi$ep300 <- as.factor(dpi$ep300)
rm(ep300); gc()

#ENHANCERS
data("enhancers", package = "xcore")
dpi$enhancer <- ""
hits <- findOverlaps(dpi, enhancers)
dpi$enhancer[hits@from] <- enhancers$name[hits@to]

#DFAM
dfam_file <- system.file("inst", "extdata", "hg38_dfam.3.1.nrph.hits.bed.gz", package = "xcore")
dfam <- rtracklayer::import.bed(dfam_file)
dpi$repeat_dfam <- ""
hits <- findOverlaps(dpi, dfam)
dpi$repeat_dfam[hits@from] <- dfam$name[hits@to]

# Detailed two direction annotation from Gencode
promoters_detailed <- xcore::gencode_two_direction_annotator(regions = dpi, gencode = gencode)

usethis::use_data(promoters_detailed, internal = FALSE, overwrite = TRUE)

