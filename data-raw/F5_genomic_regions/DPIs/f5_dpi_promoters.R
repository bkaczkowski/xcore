
# getting BED files from FANTOM5 website one for the coordinates (bed) and one with official annotation
dpi_bed_file  = "~/projects/resources/FANTOM5/hg38/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
dpi_annot_file = "~/projects/resources/FANTOM5/hg38/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz"

dpi   = rtracklayer::import.bed( dpi_bed_file )
dpi$itemRgb  = NULL ; dpi$thick  = NULL
dpi = dpi [  order(dpi$name),]

dpi_annot   = data.frame( data.table::fread(dpi_annot_file, header = T, sep = "\t"), stringsAsFactors = F)
dpi_annot   = dpi_annot [dpi_annot$X.CAGE_Peak_ID %in% dpi$name , ]
dpi_annot   = dpi_annot [ order( dpi_annot$X.CAGE_Peak_ID) , ]

table(dpi_annot$X.CAGE_Peak_ID == dpi$name)

dpi$SYMBOL_F5_annot    = sub( " .*", "", dpi_annot$Gene_symbol )
dpi$ENTREZID_F5_annot  = sub( " .*", "", dpi_annot$GeneID     )
dpi$distance_F5_annot  = dpi_annot$Distance

# ROADMAP
roadmap_hg19 = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_bed6.bed", colnames = c("chrom", "start", "end", "name", "score"))
hg19_hg38_chain = rtracklayer::import.chain("~/projects/resources/liftOver_chains/hg19ToHg38.over.chain")
roadmap = rtracklayer::liftOver(roadmap_hg19, chain = hg19_hg38_chain )
roadmap = unlist(roadmap)
dpi$roadmap = ""
hits = findOverlaps(dpi , roadmap )
dpi$roadmap[hits@from] = roadmap$name[hits@to]
rm(roadmap_hg19, roadmap ); gc()

dpi$ep300 = ""
ep300  = rtracklayer::import.bed("~/projects/resources/remap/remap2018_EP300_nr_macs2_hg38_v1_2.bed.gz" )
hits = findOverlaps(dpi ,ep300)
dpi$ep300[hits@from] = "EP300"
dpi$ep300 = as.factor(dpi$ep300)

#promoters = dpi
#usethis::use_data( promoters , internal = FALSE, overwrite = TRUE)

data("enhancers", package = "xcore")
dpi$enhancer = ""
hits = findOverlaps(dpi ,enhancers)
dpi$enhancer[hits@from] = enhancers$name[hits@to]

dfam  = rtracklayer::import.bed("data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed.gz" )
dpi$repeat_dfam = ""
hits = findOverlaps(dpi , dfam )
dpi$repeat_dfam [hits@from] = dfam$name[hits@to]

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
uscs_knownGene = GenomicFeatures::asGFF(TxDb.Hsapiens.UCSC.hg38.knownGene)
uscs_knownGene = uscs_knownGene[ uscs_knownGene$type == "mRNA"]
uscs_knownGene_promoters = GenomicRanges::promoters( uscs_knownGene , upstream = 0, downstream = 0, use.names=TRUE)
uscs_knownGene_promoters$Parent = drop( uscs_knownGene_promoters$Parent )
uscs_knownGene_promoters = uscs_knownGene_promoters [ !is.na( uscs_knownGene_promoters$Parent )]
uscs_knownGene_promoters$Parent = sub( "GeneID:" , "", uscs_knownGene_promoters$Parent )

dpi$SYMBOL_uscs     = ""
dpi$ENTREZID_uscs   = ""
dpi$distance_uscs   = 10^6

hits = GenomicRanges::distanceToNearest( dpi , uscs_knownGene_promoters , ignore.strand=FALSE)
dpi$ENTREZID_uscs [hits@from]   =  uscs_knownGene_promoters$Parent [hits@to]
dpi$distance_uscs [hits@from]   =  hits@elementMetadata$distance
library(org.Hs.eg.db)
select_first = function(x) {x[1]}
dpi$SYMBOL_uscs  =  lapply ( mapIds(org.Hs.eg.db, dpi$ENTREZID_uscs, 'SYMBOL' , 'ENTREZID' ) , select_first)

too_far = dpi$distance_uscs > 500

dpi$ENTREZID_uscs[ too_far ]  = NA
dpi$distance_uscs[ too_far ]  = NA
dpi$SYMBOL_uscs  [ too_far ]  = NA


promoters = dpi
usethis::use_data( promoters , internal = FALSE, overwrite = TRUE)

