#loading previously prepared Fantom5 annotations
# this is hg38 but hg19 coordinates can be extracted from the DPI name
load( "data-raw/F5_genomic_regions/raw_f5_regions_annotated_gencode31_UCSC_Apr25_2019.RData")
f5_regions$hg19_dpi = gsub( "hg19::|;.*", "",f5_regions$name )
dpi_annot = data.frame( name = f5_regions$name, hg19_dpi = f5_regions$hg19_dpi,stringsAsFactors = F)
dpi_annot$order = 1:nrow(dpi_annot)

# reading in FANTOM5 phase 1 annotation to get the promoter/dpi ID to "short_description" mapping.
# promoter ID provides genomic coordinates (hg19) and is stable, short_description provides link to genes and it "evolves" over time
# FANTOM5 short RNA data only provide the "short_description" (promoter), thus I need the fist columns fromt the PHASE1 annotation to get the mapping
f5_dpi_phase1_hg19_annot = read.table("~/projects/resources/FANTOM5/phase1.3/hg19.cage_peak_ann.txt.gz", header = T,sep = "\t",as.is = T)[,1:2]
dpi_id_to_promoter = data.frame( hg19_dpi = f5_dpi_phase1_hg19_annot$X00Annotation, promoter = f5_dpi_phase1_hg19_annot$short_description, stringsAsFactors = F)
rm(f5_dpi_phase1_hg19_annot)

#  merging the DPI/promoter annotations
dpi_annot = merge(dpi_annot, dpi_id_to_promoter, by = c( "hg19_dpi", "hg19_dpi"), all.x = T)
table(!is.na(dpi_annot$promoter))

# reading in F5 short RNA data
f5_shortRNA = read.table("~/projects/resources/F5_shortRNA/human.promoters.tsv", header = T,sep = "\t",as.is = T)
f5_shortRNA = f5_shortRNA[ f5_shortRNA$promoter != "Unknown" ,]
premiRNAs = tapply( X = f5_shortRNA$premiRNA ,INDEX = f5_shortRNA$promoter, FUN = paste, collapse = ";" )
miRNAs =    tapply( X = f5_shortRNA$miRNA ,INDEX = f5_shortRNA$promoter, FUN = paste, collapse = ";" )
table(names(premiRNAs) == names(miRNAs ))
miRNA_annot = data.frame( promoter = names(miRNAs ), miRNA =  miRNAs, premiRNA = premiRNAs)
rm(f5_shortRNA, miRNAs, premiRNAs)
miRNA_annot = miRNA_annot[ miRNA_annot$promoter %in% dpi_id_to_promoter$promoter,]
table(miRNA_annot$promoter %in% dpi_id_to_promoter$promoter)

# merging the F5 short RNA data with previously prepared DPI annotation table
dpi_annot = merge(dpi_annot, miRNA_annot, by = c( "promoter", "promoter"), all.x = T )
dpi_annot$miRNA [is.na(dpi_annot$miRNA) ] = ""
dpi_annot$premiRNA [is.na(dpi_annot$premiRNA) ] = ""
dpi_annot = dpi_annot [ order(dpi_annot$order), ]
dpi_annot$order = NULL

f5_shortRNA_microRNA_annotation = dpi_annot
# changing the object name to something for unique and descriptive
# so it doesn't overwrite anything after loading to R enviroment later and can be recognized
save(f5_shortRNA_microRNA_annotation, file = "data-raw/F5_genomic_regions/f5_shortRNA_microRNA_annotation.RData")

write.csv(dpi_annot[ dpi_annot$miRNA != "" , ], file = "data-raw/F5_genomic_regions/f5_shortRNA_microRNA_annotation.csv", row.names = F)
dpi_annot = read.csv("data-raw/F5_genomic_regions/f5_shortRNA_microRNA_annotation.csv")




