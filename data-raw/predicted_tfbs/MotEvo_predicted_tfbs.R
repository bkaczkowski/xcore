# Bogumil Kaczkowski
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions
# with predicted Transcribtion Factor Binding Sites (TFBS)

# 1c) read in the MotEvo tfbs into R
tfbs = rtracklayer::import("~/projects/resources/MotEvo_TFBS_predictions/hg38.sites.bed.bz2")
unique_factors = unique(tfbs$name)
unique_factors = unique_factors [ order( unique_factors )]

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
overlap_mat_promoters = matrix (0L, nrow = length(promoters_ext500), ncol = length(unique_factors) )
colnames(overlap_mat_promoters) = unique_factors
rownames(overlap_mat_promoters) = promoters_ext500$name

overlap_mat_enhancers = matrix (0L, nrow = length(enhancers_ext500), ncol = length(unique_factors) )
colnames(overlap_mat_enhancers) = unique_factors
rownames(overlap_mat_enhancers) = enhancers_ext500$name

total_number_of_peaks = rep( NA,length(unique_factors) )
names( total_number_of_peaks ) = unique_factors
gc()

for ( j in 1:length(unique_factors)){

  query_tf = tfbs[tfbs$name == unique_factors [j] ]
  total_number_of_peaks [j] = length(query_tf)

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

rm(hits, query_tf, j, tfbs)
gc()

save( overlap_mat_enhancers, file = "data-raw/chip_atlas/overlap_mat_enhancers.rda")
save( overlap_mat_promoters, file = "data-raw/chip_atlas/overlap_mat_promoters.rda")

MotEvo_meta = data.frame( total_number_of_peaks = total_number_of_peaks,
                          number_of_promoters_with_peak   =  colSums(overlap_mat_promoters >0),
                          number_of_enhancers_with_peak   =  colSums(overlap_mat_enhancers >0) )
  )

# 5) Gene level summarization

empty = promoters_ext500$SYMBOL == ""
promoters_ext500_for_gene_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_gene_collapsing = overlap_mat_promoters[ !empty ,]

empty = promoters_ext500$ENTREZID == ""
promoters_ext500_for_entrez_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_entrez_collapsing = overlap_mat_promoters[ !empty ,]

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


MotEvo_meta$number_of_geneSymbols_with_peak =  colSums(overlap_mat_symbol >0 )
MotEvo_meta$number_of_geneEntrez_with_peak  =  colSums(overlap_mat_entrez >0 )


# Export data
MotEvo_entrez = Matrix::drop0( overlap_mat_entrez )
MotEvo_symbol = Matrix::drop0( overlap_mat_symbol )
MotEvo_promoters = Matrix::drop0( overlap_mat_promoters )
MotEvo_enhancers = Matrix::drop0( overlap_mat_enhancers )

usethis::use_data( MotEvo_entrez   , internal = FALSE, overwrite = TRUE)
usethis::use_data( MotEvo_symbol   , internal = FALSE, overwrite = TRUE)
usethis::use_data( MotEvo_promoters, internal = FALSE, overwrite = TRUE)
usethis::use_data( MotEvo_enhancers, internal = FALSE, overwrite = TRUE)
#usethis::use_data( MotEvo_meta     , internal = FALSE, overwrite = TRUE)

# 6) Plotting
number_of_TFBS_per_gene      =  rowSums(overlap_mat_symbol >0 )
number_of_TFBS_per_entrez    =  rowSums(overlap_mat_entrez >0 )

dev.off()
pdf("data-raw/predicted_tfbs/MotEvo_NumberOfPeaks_vs_NumberOfFeatures_with_peaks.pdf", width = 8, height = 6)
with( MotEvo_meta , smoothScatter(total_number_of_peaks /1000 , number_of_promoters_with_peak/1000  ,xlim = c(0,110 ), main = "DPIs" ))
with( MotEvo_meta , smoothScatter(total_number_of_peaks /1000 , number_of_enhancers_with_peak/1000  ,xlim = c(0,110 ), main = "enhancers"  ))
with( MotEvo_meta , smoothScatter(total_number_of_peaks /1000 , number_of_geneSymbols_with_peak/1000,xlim = c(0,110 ), main = "Genes(SYMBOLS)" ))
with( MotEvo_meta , smoothScatter(total_number_of_peaks /1000 , number_of_geneEntrez_with_peak/1000 ,xlim = c(0,110 ), main = "Genes(ENTREZID)" ))
dev.off()

pdf("data-raw/predicted_tfbs/MotEvo_NumberOfTFBS_per_gene.pdf", width = 8, height = 6)
plot(sort(number_of_TFBS_per_gene, decreasing = T), type = "l",
     ylab = "number of TFBS", main = "Number of TF that have TFBS on given gene's promoter ", xlab = "Genes index")
dev.off()

intersections_gene_level_sel = overlap_mat_symbol [ rowSums( overlap_mat_symbol>0 )>= 10 ,]
intersections_gene_level_sel = intersections_gene_level_sel [ , colSums(intersections_gene_level_sel > 0) > 500 ]
intersections_gene_level_sel = ifelse(intersections_gene_level_sel > 0, yes = 1L, no = 0L)

dev.off()
pdf("data-raw/predicted_tfbs/MotEvo_upset_top30_freq_chip.pdf", width = 15, height = 10)
UpSetR::upset(data.frame(intersections_gene_level_sel) , nsets = 30, order.by = "freq")
dev.off()

intersections_gene_level_sel = overlap_mat_symbol [ rowSums( overlap_mat_symbol>0 )>= 20 ,]
intersections_gene_level_sel = intersections_gene_level_sel [ , colSums(intersections_gene_level_sel > 0) > 200 ]
intersections_gene_level_sel = intersections_gene_level_sel [ rowSums( intersections_gene_level_sel>0 )>= 35 ,]

pheatmap::pheatmap(t(intersections_gene_level_sel), cluster_rows = T,cellwidth = 10,cellheight = 10,
                file = "data-raw/predicted_tfbs/MotEvo_heatmap_most_dense_overlap.pdf")

