load("data-raw/F5_genomic_regions/raw_f5_regions_annotated_gencode31_UCSC_Apr25_2019.RData")

#adding miRNA annotations
load("data-raw/F5_genomic_regions/f5_shortRNA_microRNA_annotation.RData")

table(f5_regions$name == f5_shortRNA_microRNA_annotation$name)
f5_regions$premiRNA = as.character(f5_shortRNA_microRNA_annotation$premiRNA)

f5_regions$transcriptId = NULL
f5_regions$geneId = NULL
f5_regions$ucsc_ENSEMBL = NULL
f5_regions$ucsc_GENENAME = NULL
f5_regions$itemRgb = NULL
f5_regions$thick.width = NULL
f5_regions$thick.start = NULL
f5_regions$thick.end = NULL

promoters = f5_regions$annotation %in% c("5' UTR" , "Promoter")
gene_body = f5_regions$annotation %in% c("Intron", "3' UTR" , "Exon", "Downstream" )
distal    = f5_regions$annotation %in% c("Distal Intergenic"  )


f5_regions$overlapping_gene = ""
f5_regions$overlapping_gene[gene_body] = f5_regions$SYMBOL[gene_body]
f5_regions$SYMBOL[gene_body] = ""

f5_regions$nearest_gene = ""
f5_regions$nearest_gene[distal] = f5_regions$SYMBOL[distal]
f5_regions$SYMBOL[distal] = ""

f5_regions$SYMBOL    [ !promoters ] = ""
f5_regions$gene_type [ !promoters ] = NA
f5_regions$ENTREZID  [ !promoters ] = ""
f5_regions$GENENAME  [ !promoters ] = ""

f5_regions$ucsc_ENTREZID [ ! f5_regions$ucsc_annotation %in% c("5' UTR" , "Promoter") ] = ""
f5_regions$ucsc_SYMBOL [ ! f5_regions$ucsc_annotation %in% c("5' UTR" , "Promoter") ] = ""


table(f5_regions$ucsc_ENTREZID != "")
table(f5_regions$ucsc_SYMBOL != "")

table(f5_regions$ENTREZID != "")
table(f5_regions$SYMBOL != "")


sel = f5_regions$ucsc_ENTREZID != ""
table(f5_regions$ucsc_ENTREZID[sel] == f5_regions$ENTREZID[sel])

sel = f5_regions$ucsc_SYMBOL != ""
table(f5_regions$ucsc_SYMBOL[sel] == f5_regions$SYMBOL[sel])

length(unique(f5_regions$ucsc_SYMBOL))
length(unique(f5_regions$SYMBOL))

length(unique(f5_regions$ucsc_ENTREZID))
length(unique(f5_regions$ENTREZID))

f5_regions$roadmap_type = as.factor(f5_regions$roadmap_type)
f5_regions$type = as.factor(f5_regions$type)
f5_regions$annotation = as.factor(f5_regions$annotation)
f5_regions$gene_type = as.factor(f5_regions$gene_type)
f5_regions$ep300 = as.factor(f5_regions$ep300)

save( f5_regions , file = "data-raw/F5_genomic_regions/f5_regions.RData")

#png("dist.png")
#smoothScatter(log10(f5_regions$distanceToTSS+1), log10(f5_regions$ucsc_distanceToTSS+1))
#dev.off()
