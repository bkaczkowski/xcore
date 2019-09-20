# Bogumil Kaczkowski
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions
# with ChipAtlas, Chip-Seq peak data-base

# Following steps are implemented"
# 1) Chip-seq peaks
#   1a) Download ChipAtlas (all) peaks from http://chip-atlas.org/peak_browser
#   1b) SHELL commands to process the downloaded BED file
#   1c) read in the peak into R
# 2) get F5 promoter region
# 3) create intersection matrix
# 4) Fix the meta data for the chipseq experiments

# 1a) Get ChipAtlas (all) peaks (start direclty from downloaded bed, chip_atlas.sh is not needed)
# Oth.ALL.05.AllAg.AllCell.bed was obtained by downloading all Transcription Factor peaks [TF and others(9868)]
# from http://chip-atlas.org/peak_browser the file is ~ 72G uncompressed and takes some hours to download
# the peaks coordinates are of the hg19 genome https://github.com/inutano/chip-atlas/wiki

# 1b) SHELL commands to process the downloaded BED file (also in the bash file)
# extract SRX and remove the unnecessary data from 4th column
#zcat Oth.ALL.05.AllAg.AllCell.bed.gz  | awk '{gsub(";.*|ID=","",$4)}1' |  awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > Oth.ALL.05.AllAg.AllCell_SRX_only.bed
# remove the first line (dummy header) and last row (faulty row)
#tail -n +2  Oth.ALL.05.AllAg.AllCell_SRX_only.bed > Oth.ALL.05.AllAg.AllCell_SRX_only_hg19.bed
# rm Oth.ALL.05.AllAg.AllCell_SRX_only.bed
#liftOver to hg38
# ~/bin/liftOver  Oth.ALL.05.AllAg.AllCell_SRX_only_hg19.bed  ~/projects/resources/liftOver_chains/hg19ToHg38.over.chain Oth.ALL.05.AllAg.AllCell_SRX_only_hg19ToHg38_liftOver.bed unMapped

# 1c) read in the peak into R
chip_atlas = rtracklayer::import("~/projects/resources/chip_atlas/Oth.ALL.05.AllAg.AllCell_SRX_only_hg19ToHg38_liftOver.bed")
chip_atlas$score = as.integer(chip_atlas$score)
unique_srxIDs = unique(chip_atlas$name)
unique_srxIDs = unique_srxIDs [ order( unique_srxIDs )]

# 2) get F5 promoter region
data("promoters", package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500

data("enhancers", package = "xcore")
enhancers_ext500 = enhancers
GenomicRanges::start(enhancers_ext500)   = GenomicRanges::start (enhancers)   - 500
GenomicRanges::end  (enhancers_ext500)   = GenomicRanges::end   (enhancers)   + 500

rm( enhancers, promoters)
gc()

# 3) create intersection matrises
overlap_mat_promoters = matrix (0L, nrow = length(promoters_ext500), ncol = length(unique_srxIDs) )
colnames(overlap_mat_promoters) = unique_srxIDs
rownames(overlap_mat_promoters) = promoters_ext500$name

overlap_mat_enhancers = matrix (0L, nrow = length(enhancers_ext500), ncol = length(unique_srxIDs) )
colnames(overlap_mat_enhancers) = unique_srxIDs
rownames(overlap_mat_enhancers) = enhancers_ext500$name

system( "mkdir ~/projects/resources/chip_atlas/split_by_SRX/")
export_dir = "~/projects/resources/chip_atlas/split_by_SRX/"

total_number_of_peaks = rep( NA,length(unique_srxIDs) )
names( total_number_of_peaks ) = unique_srxIDs
gc()

for ( j in 1:length(unique_srxIDs)){

  query_tf = chip_atlas[chip_atlas$name == unique_srxIDs [j] ]
  total_number_of_peaks [j] = length(query_tf)
  rtracklayer::export.bed(con = paste( export_dir, unique_srxIDs [j], ".bed",sep = ""), query_tf)

  hits = GenomicRanges::findOverlaps(promoters_ext500 ,query_tf)
  hits = S4Vectors::DataFrame(from =hits@from , to = hits@to,  score = query_tf$score[ hits@to])
  hits = hits [ order ( hits$score , decreasing = T) , ]
  hits = hits [ !duplicated(hits$from) ,]
  overlap_mat_promoters[ hits$from , j] =   query_tf$score[ hits$to ]

  hits = GenomicRanges::findOverlaps(enhancers_ext500 ,query_tf)
  hits = S4Vectors::DataFrame(from =hits@from , to = hits@to,  score = query_tf$score[ hits@to])
  hits = hits [ order ( hits$score , decreasing = T) , ]
  hits = hits [ !duplicated(hits$from) ,]
  overlap_mat_enhancers[ hits$from , j] =   query_tf$score[ hits$to ]

  print( j)
}

rm(hits, query_tf, j, chip_atlas, export_dir)
gc()

# 4) Meta data for the chipseq experiments
# grep -w hg19 experimentList.tab.txt | grep -w 'TFs and others'  > experimentList.tab_TF_hg19.txt
sxr_meta = read.delim("~/projects/resources/chip_atlas/experimentList.tab_TF_hg19.txt", sep ="\t", header = F, quote = "\"")
dim(sxr_meta)
sxr_meta = sxr_meta [ sxr_meta$V1 %in%  unique_srxIDs, ]
sxr_meta = sxr_meta [ , c( 1, 4:7)]
colnames(sxr_meta) = c("id", "TF", "origin", "cell", "description")
sxr_meta = data.frame(apply( sxr_meta, 2, as.character),stringsAsFactors = F)
sxr_meta = sxr_meta [ order( sxr_meta$id) ,]
table(sxr_meta$origin)

sxr_meta$name = paste(sxr_meta$TF, sxr_meta$origin, sxr_meta$cell, sxr_meta$id, sep = "_")
sxr_meta$name = gsub( " " , ".", sxr_meta$name)

table(sxr_meta$id == colnames(overlap_mat_enhancers))
table(sxr_meta$id == colnames(overlap_mat_promoters))

colnames(overlap_mat_enhancers) = sxr_meta$name
colnames(overlap_mat_promoters) = sxr_meta$name

save( overlap_mat_enhancers, file = "data-raw/chip_atlas/overlap_mat_enhancers.rda")
save( overlap_mat_promoters, file = "data-raw/chip_atlas/overlap_mat_promoters.rda")

sxr_meta$number_of_promoters_with_peak   =  colSums(overlap_mat_promoters >0 )
sxr_meta$number_of_enhancers_with_peak   =  colSums(overlap_mat_enhancers >0 )

# 5) Gene level summarization

empty = promoters_ext500$SYMBOL == ""
promoters_ext500_for_gene_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_gene_collapsing = overlap_mat_promoters[ !empty ,]

empty = promoters_ext500$ENTREZID == ""
promoters_ext500_for_entrez_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_entrez_collapsing = overlap_mat_promoters[ !empty ,]


rm(overlap_mat_promoters, overlap_mat_enhancers ); gc()

max_score_symbol = function(x) {
  tapply( x, INDEX = promoters_ext500_for_gene_collapsing$SYMBOL, FUN = max)
}
max_score_entrez = function(x) {
  tapply( x, INDEX = promoters_ext500_for_entrez_collapsing$ENTREZID, FUN = max)
}

overlap_mat_symbol = apply( overlap_mat_promoters_for_gene_collapsing,
                            2, max_score_symbol)
overlap_mat_entrez = apply( overlap_mat_promoters_for_entrez_collapsing,
                            2, max_score_entrez)

save( overlap_mat_symbol, file = "data-raw/chip_atlas/overlap_mat_symbol.rda")
save( overlap_mat_entrez, file = "data-raw/chip_atlas/overlap_mat_entrez.rda")
sxr_meta$number_of_geneSymbols_with_peak =  colSums(overlap_mat_symbol >0 )
sxr_meta$number_of_geneEntrez_with_peak  =  colSums(overlap_mat_entrez >0 )
sxr_meta$total_number_of_peaks = total_number_of_peaks

# Export data
chip_atlas_entrez = overlap_mat_entrez
chip_atlas_symbol = overlap_mat_symbol
chip_atlas_meta   = sxr_meta

chip_atlas_promoters   =overlap_mat_promoters
chip_atlas_enhancers   =overlap_mat_enhancers

usethis::use_data( chip_atlas_entrez   , internal = FALSE, overwrite = TRUE)
usethis::use_data( chip_atlas_symbol   , internal = FALSE, overwrite = TRUE)
usethis::use_data( chip_atlas_meta     , internal = FALSE, overwrite = TRUE)
usethis::use_data( chip_atlas_promoters, internal = FALSE, overwrite = TRUE)
usethis::use_data( chip_atlas_enhancers, internal = FALSE, overwrite = TRUE)

# 6) Plotting
number_of_chipExp_per_gene      =  rowSums(overlap_mat_symbol >0 )
number_of_chipExp_per_entrez    =  rowSums(overlap_mat_entrez >0 )

dev.off()
pdf("data-raw/chip_atlas/NumberOfPeaks_vs_NumberOfFeatures_with_peaks.pdf", width = 8, height = 6)
with( sxr_meta , smoothScatter(total_number_of_peaks /1000 , number_of_promoters_with_peak/1000  ,xlim = c(0,110 ), main = "DPIs" ))
with( sxr_meta , smoothScatter(total_number_of_peaks /1000 , number_of_enhancers_with_peak/1000  ,xlim = c(0,110 ), main = "enhancers"  ))
with( sxr_meta , smoothScatter(total_number_of_peaks /1000 , number_of_geneSymbols_with_peak/1000,xlim = c(0,110 ), main = "Genes(SYMBOLS)" ))
with( sxr_meta , smoothScatter(total_number_of_peaks /1000 , number_of_geneEntrez_with_peak/1000 ,xlim = c(0,110 ), main = "Genes(ENTREZID)" ))
dev.off()

pdf("data-raw/chip_atlas/NumberOfPeaks(Experiments)_per_gene.pdf", width = 8, height = 6)
plot(sort(number_of_chipExp_per_gene, decreasing = T), type = "l",
     ylab = "number of peaks", main = "Number of Chip-Seq experiments \nwith peak on given gene promoter ", xlab = "Genes index")
dev.off()

intersections_gene_level_sel = overlap_mat_symbol [ rowSums( overlap_mat_symbol>0 )>= 100,]
intersections_gene_level_sel = intersections_gene_level_sel [ rowSums( intersections_gene_level_sel>0 )< 3500,]
intersections_gene_level_sel = intersections_gene_level_sel [ , colSums(intersections_gene_level_sel > 0) > 500 ]
intersections_gene_level_sel = ifelse(intersections_gene_level_sel > 0, yes = 1L, no = 0L)

dev.off()
pdf("data-raw/chip_atlas/upset_top30_freq_chip.pdf", width = 15, height = 10)
UpSetR::upset(data.frame(intersections_gene_level_sel) , nsets = 30, order.by = "freq")
dev.off()

dev.off()
pdf("data-raw/chip_atlas/upset_top30_freq_genes.pdf", width = 15, height = 10)
UpSetR::upset(data.frame(t(intersections_gene_level_sel) ), nsets = 30, order.by = "freq")
dev.off()

intersections_gene_level_sel = overlap_mat_symbol [ rowSums(overlap_mat_symbol>500)>= 500,]
intersections_gene_level_sel = intersections_gene_level_sel [ rowSums( intersections_gene_level_sel>0 )< 3500,]
intersections_gene_level_sel = intersections_gene_level_sel [ , colSums(intersections_gene_level_sel >= 500) > 200 ]
intersections_gene_level_sel = intersections_gene_level_sel [ rowSums(intersections_gene_level_sel>500)>= 900,]
intersections_gene_level_sel = intersections_gene_level_sel [ , colSums(intersections_gene_level_sel >= 500) > 30 ]
pheatmap::pheatmap(t(intersections_gene_level_sel)[seq(1,860, 5),], cluster_rows = T,cellwidth = 10,cellheight = 10,
                file = "data-raw/chip_atlas/heatmap_most_dense_overlap.pdf")

