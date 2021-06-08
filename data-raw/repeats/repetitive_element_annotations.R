
# Bogumil Kaczkowski August 19, 2019
# A script is to make a BED file with repeats with Dfam database

dfam_nrph = data.table::fread("~/projects/resources/Dfam/hg38_dfam.3.1.nrph.hits.gz")

dfam_nrph_bed = data.frame( chr = dfam_nrph$`#seq_name` ,
                            start = dfam_nrph$`ali-st` -1 , # adjusting "start" to the BED format which is 0-based
                            end = dfam_nrph$`ali-en`,
                            name = dfam_nrph$family_name,
                            score = -log10(dfam_nrph$`e-value`),
                            strand = dfam_nrph$strand)

# adjusting "start" and "end" of the alighnments to the BED format for "-" strand
minus_strand = dfam_nrph_bed$strand == "-"
dfam_nrph_bed$start[ minus_strand ] = dfam_nrph$`ali-en` [ minus_strand] -1  # adjusting "start" to the BED format which is 0-based
dfam_nrph_bed$end  [ minus_strand ] = dfam_nrph$`ali-st` [ minus_strand]

# scaling the score to 1-1000 range
dfam_nrph_bed$score [ is.infinite(dfam_nrph_bed$score) ] = 320
dfam_nrph_bed$score = dfam_nrph_bed$score + 10 # offsetting the few negative values to positive teritory
dfam_nrph_bed$score = dfam_nrph_bed$score * 3 # multiplying by 3 to transform 1-300 range in to 1-1000 range
dfam_nrph_bed$score = round(dfam_nrph_bed$score)
summary(dfam_nrph_bed$score)
hist(dfam_nrph_bed$score)

options(scipen = 999)
write.table( dfam_nrph_bed, file = "data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed",
             col.names = F, row.names = F, quote = F, sep = "\t" )
#write.table( dfam_nrph_bed, file = "~/projects/resources/Dfam/hg38_dfam.3.1.nrph.hits.bed",
#             col.names = F, row.names = F, quote = F, sep = "\t" )
options(scipen = 15)

system( "gzip data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed")


