
f5_dpi   = rtracklayer::import.bed("~/projects/resources/FANTOM5/hg38/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz" )
f5_dpi@elementMetadata$type = "DPI"

f5_enh   = rtracklayer::import.bed("~/projects/resources/FANTOM5/hg38/F5.hg38.enhancers.bed.gz" )
f5_enh@elementMetadata$blocks = NULL
f5_enh@elementMetadata$type = "enhancer"

f5_regions = c(f5_dpi , f5_enh)
rm(f5_dpi, f5_enh)

f5_regions$itemRgb = NULL; f5_regions$thick = NULL

gencode = rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")



f5_enh_nearest_promoter = gencode_nearest_promoter( regions = f5_enh , gencode = gencode )


gencode_gene        = gencode[ gencode$type == "gene"]

GenomicRanges::distanceToNearest( f5_regions , gencode_gene, ignore.strand=FALSE)

regions_tmp = f5_regions
strand(regions_tmp) = "+"
nearest_plus  = GenomicRanges::distanceToNearest( regions_tmp, gencode_gene, ignore.strand=FALSE)
regions_annotated_plus  = mcols(regions_annotated_plus)

strand(regions_tmp) = "-"
nearest_minus = GenomicRanges::distanceToNearest( regions_tmp, gencode_gene, ignore.strand=FALSE)


f5_regions_annotated_one_direction   = gencode_one_direction_annotator( regions = f5_regions, gencode = gencode)

f5_regions_annotated_two_direction   = gencode_two_direction_annotator( regions = f5_regions, gencode = gencode)

f5_regions_annotated_full_annotation = gencode_strandless_worker( regions = f5_regions_annotated_two_direction, gencode = gencode)

f5_regions_annotated_full_annotation$annotatation [ f5_regions_annotated_full_annotation$annotatation=="five_prime_UTR"] = "promoter"
f5_regions_annotated_full_annotation$opposite_strand_annotatation [ f5_regions_annotated_full_annotation$opposite_strand_annotatation=="five_prime_UTR"] = "promoter"

f5_regions_annotated_full_annotation$annotatation_mod = with (f5_regions_annotated_full_annotation, paste(annotatation, "(sense);", opposite_strand_annotatation, "(anti-sense)",sep = ""))

f5_regions_annotated_full_annotation$annotatation_mod  = sub(";\\(anti-sense\\)", ";" , f5_regions_annotated_full_annotation$annotatation_mod )
f5_regions_annotated_full_annotation$annotatation_mod  = sub("^\\(sense\\);", ";" , f5_regions_annotated_full_annotation$annotatation_mod )
sort(table(f5_regions_annotated_full_annotation$annotatation_mod[f5_regions_annotated_full_annotation$type=="DPI"]))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
uscs_knownGene = GenomicFeatures::asGFF(TxDb.Hsapiens.UCSC.hg38.knownGene)

uscs_knownGene_gene        = uscs_knownGene[ uscs_knownGene$type == "gene"]
uscs_knownGene_exon        = uscs_knownGene[ uscs_knownGene$type == "exon"]
uscs_knownGene_transcripts = uscs_knownGene[ uscs_knownGene$type == "transcript"]

hits = findOverlaps( uscs_knownGene_cds, uscs_knownGene_exon )
f5_regions$roadmap[hits@from] = roadmap$name[hits@to]
f5_regions$roadmap_type    = sub("\\|.*" , "" , f5_regions$roadmap)


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
