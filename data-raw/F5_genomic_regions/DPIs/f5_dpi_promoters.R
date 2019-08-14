
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
hits = findOverlaps(all_enh ,ep300)
dpi$ep300[hits@from] = "EP300"
dpi$ep300 = as.factor(dpi$ep300)

promoters = dpi
usethis::use_data( promoters , internal = FALSE)

data("enhancers", package = "xcore")
dpi$enhancer = ""
hits = findOverlaps(dpi ,enhancers)
dpi$enhancer[hits@from] = enhancers$name[hits@to]



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
dpi$distance_uscs [hits@from]   =  hits@elementMetadata$distance   [hits@to]

dpi$ENTREZID_uscs = as.character( dpi$ENTREZID_uscs )
dpi$ENTREZID_uscs = sub( "GeneID:" , "" , dpi$ENTREZID_uscs )





# Original FANTOM5 enhancers based on mapping to hg19 ####
f5_enhancers_hg19 = rtracklayer::import.bed("~/projects/resources/FANTOM5/human_permissive_enhancers_phase_1_and_2.bed.gz" )
# lifting over to hg38
hg19_hg38_chain = rtracklayer::import.chain("~/projects/resources/liftOver_chains/hg19ToHg38.over.chain")
f5_enhancers_hg19_liftOver_to_hg38 = rtracklayer::liftOver(f5_enhancers_hg19, chain = hg19_hg38_chain )
f5_enhancers_hg19_liftOver_to_hg38 = unlist(f5_enhancers_hg19_liftOver_to_hg38)
rm(f5_enhancers_hg19)

f5_regions@elementMetadata$hg19_enhancer = ""
library(GenomicRanges)
hits = findOverlaps(f5_regions ,f5_enhancers_hg19_liftOver_to_hg38 )
f5_regions@elementMetadata$hg19_enhancer[hits@from] = f5_enhancers_hg19_liftOver_to_hg38$name[hits@to]

#.... ROADMAP annotations...
roadmap_hg19 = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_bed6.bed", colnames = c("chrom", "start", "end", "name", "score"))
roadmap = rtracklayer::liftOver(roadmap_hg19, chain = hg19_hg38_chain )

roadmap = unlist(roadmap)

f5_regions$roadmap = ""
hits = findOverlaps(f5_regions , roadmap )
f5_regions$roadmap[hits@from] = roadmap$name[hits@to]
f5_regions$roadmap_type    = sub("\\|.*" , "" , f5_regions$roadmap)
rm(roadmap_hg19, roadmap ); gc()

#.... EP300 peaks...
f5_regions$ep300 = ""
ep300  = rtracklayer::import.bed("~/projects/resources/remap/remap2018_EP300_nr_macs2_hg38_v1_2.bed.gz" )
hits = findOverlaps(f5_regions ,ep300)
f5_regions$ep300[hits@from] = "EP300"

#ANNOTATION
library(ChIPseeker)

gencode = rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
gencode_TxDb = GenomicFeatures::makeTxDbFromGRanges(gencode, drop.stop.codons = TRUE)

gencode_anno <- ChIPseeker::annotatePeak( f5_regions, tssRegion=c(-500, 500),
                                          TxDb=gencode_TxDb,
                                          level = "transcript", sameStrand = TRUE)
gencode_anno = gencode_anno@anno
gencode_anno$annotation = sub(" \\(.*" , "" , gencode_anno$annotation)
gencode_anno_df = as.data.frame(gencode_anno@elementMetadata )

gencode_annotation = data.frame(cbind( gencode$gene_id , gencode$gene_name, gencode$gene_type ), stringsAsFactors = F)
colnames(gencode_annotation) = c( "geneId" , "SYMBOL", "gene_type" )
gencode_annotation = gencode_annotation[ !duplicated(gencode_annotation$geneId),]
library(org.Hs.eg.db)
gencode_annotation$ENTREZID =  mapIds(org.Hs.eg.db, gencode_annotation$SYMBOL, 'ENTREZID', 'SYMBOL')
gencode_annotation$GENENAME =  mapIds(org.Hs.eg.db, gencode_annotation$SYMBOL, 'GENENAME', 'SYMBOL')

gencode_annotation = dplyr::left_join( gencode_anno_df , gencode_annotation , by = c("geneId", "geneId" ) )

table( gencode_annotation$name == f5_regions$name  )

rm(gencode_TxDb, gencode,gencode_anno, gencode_anno_df)

#####

#### UCSC.hg38.knownGene annotation #####
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
ucsc_anno <- ChIPseeker::annotatePeak( f5_regions, tssRegion=c(-500, 500),
                                       TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db",
                                       level = "transcript", sameStrand = TRUE)
ucsc_anno = ucsc_anno@anno
table(ucsc_anno$name == f5_regions$name)
ucsc_anno$annotation = sub(" \\(.*" , "" , ucsc_anno$annotation)


# ADDING THE ANNOTATION
table( gencode_annotation$name == f5_regions$name  )
table(f5_regions$name == ucsc_anno$name)

f5_regions@elementMetadata = DataFrame ( gencode_annotation[ , - (13:17) ])
f5_regions$ucsc_annotation = ucsc_anno$annotation
f5_regions$ucsc_ENTREZID = ucsc_anno$geneId
f5_regions$ucsc_distanceToTSS = ucsc_anno$distanceToTSS
f5_regions$ucsc_ENSEMBL  = ucsc_anno$ENSEMBL
f5_regions$ucsc_SYMBOL   = ucsc_anno$SYMBOL
f5_regions$ucsc_GENENAME = ucsc_anno$GENENAME

table(f5_regions$ENTREZID  == f5_regions$ucsc_ENTREZID , f5_regions$SYMBOL == f5_regions$ucsc_SYMBOL)
table(f5_regions$annotation, f5_regions$ucsc_annotation)

f5_regions = f5_regions[ order(f5_regions$name) , ]
save(f5_regions, file = "data-raw/F5_genomic_regions/raw_f5_regions_ChIPseeker_annotated_gencode31_UCSC_Apr25_2019.RData")
