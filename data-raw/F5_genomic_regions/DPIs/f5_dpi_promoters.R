
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

# GENCODE 31 annotation
gencode = rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
dpi_gencode = xcore::gencode_nearest_promoter_same_strand(regions = dpi , gencode = gencode, cut_off_distance = 500)
dpi$SYMBOL_gencode    = dpi_gencode$nearest_symbol
library(org.Hs.eg.db)
select_first = function(x) {x[1]}
dpi$ENTREZID_gencode  = as.character( lapply ( mapIds(org.Hs.eg.db, dpi$SYMBOL_gencode, 'ENTREZID' , 'SYMBOL'  ) , select_first))
dpi$ENTREZID_gencode[ dpi$ENTREZID_gencode =="NULL"] =""
dpi$ENTREZID_gencode[ is.na(dpi$ENTREZID_gencode) ] =""
dpi$gene_type_gencode = as.factor(dpi_gencode$nearest_gene_type)
dpi$distance_gencode  = dpi_gencode$nearest_distance

# UCSC.hg38.knownGene annotation
library(GenomicRanges)
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
dpi$SYMBOL_uscs  = as.character( lapply ( mapIds(org.Hs.eg.db, dpi$ENTREZID_uscs, 'SYMBOL' , 'ENTREZID' ) , select_first))
dpi$SYMBOL_uscs[ dpi$SYMBOL_uscs =="NULL"] =""

too_far = dpi$distance_uscs > 500

dpi$ENTREZID_uscs[ too_far ]  = ""
dpi$SYMBOL_uscs  [ too_far ]  = ""

promoters = dpi
promoters$distance_F5_annot = NULL
promoters$distance_gencode = NULL
promoters$distance_uscs = NULL

usethis::use_data( promoters , internal = FALSE, overwrite = TRUE)

# ROADMAP
roadmap = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_hg38.bed", colnames = c("chrom", "start", "end", "name", "score"))
dpi$roadmap = ""
hits = findOverlaps(dpi , roadmap )
dpi$roadmap[hits@from] = roadmap$name[hits@to]
rm(roadmap_hg19, roadmap ); gc()

#EP300
dpi$ep300 = ""
ep300  = rtracklayer::import.bed("~/projects/resources/remap/remap2018_EP300_nr_macs2_hg38_v1_2.bed.gz" )
hits = findOverlaps(dpi ,ep300)
dpi$ep300[hits@from] = "EP300"
dpi$ep300 = as.factor(dpi$ep300)

#ENHANCERS
data("enhancers", package = "xcore")
dpi$enhancer = ""
hits = findOverlaps(dpi ,enhancers)
dpi$enhancer[hits@from] = enhancers$name[hits@to]

#DFAM
dfam  = rtracklayer::import.bed("data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed.gz" )
dpi$repeat_dfam = ""
hits = findOverlaps(dpi , dfam )
dpi$repeat_dfam [hits@from] = dfam$name[hits@to]

dpi = xcore::gencode_two_direction_annotator(regions = dpi, gencode = gencode)

promoters_detailed = dpi
usethis::use_data( promoters_detailed , internal = FALSE, overwrite = TRUE)

# some checks, stats and overview
table(dpi$SYMBOL_F5_annot == dpi$SYMBOL_gencode)
table(dpi$ENTREZID_F5_annot == dpi$ENTREZID_gencode)
table(dpi$ENTREZID_uscs == dpi$ENTREZID_gencode)

table(dpi$SYMBOL_F5_annot != "")
table(dpi$SYMBOL_gencode != "")
table(dpi$SYMBOL_uscs != "")

table(dpi$ENTREZID_F5_annot != "")
table(dpi$ENTREZID_gencode != "")
table(dpi$ENTREZID_uscs    != "")

length(unique(dpi$SYMBOL_F5_annot))
length(unique(dpi$SYMBOL_gencode))
length(unique(dpi$SYMBOL_uscs))

length(unique(dpi$ENTREZID_F5_annot))
length(unique(dpi$ENTREZID_gencode))
length(unique(dpi$ENTREZID_uscs))

table( dpi$SYMBOL_gencode != "", dpi$annotation == "promoter")

s = dpi$SYMBOL_gencode != ""
table(dpi$SYMBOL_gencode[s]  == dpi$symbol[s])

dpi [ dpi$SYMBOL_gencode  != dpi$symbol & s ]

table(dpi$SYMBOL_F5_annot[s] == dpi$symbol[s])
table(dpi$SYMBOL_F5_annot[s] == dpi$SYMBOL_gencode[s])

table(dpi$SYMBOL_uscs[s]     == dpi$symbol[s])
table(dpi$SYMBOL_uscs[s]     == dpi$SYMBOL_gencode[s])

table(dpi$SYMBOL_uscs[s]     == dpi$SYMBOL_F5_annot[s])

#there are some annotation cases that are confusing
dpi [ dpi$SYMBOL_gencode  != dpi$symbol & s ,c(3,6,10,18)]
