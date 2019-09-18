# Bogumil Kaczkowski
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions with ReMap Chip Seq data-base

# 1a) Get ReMap non redundant peaks
# wget http://tagc.univ-mrs.fr/remap/download/remap2018/hg38/MACS/remap2018_nr_macs2_hg38_v1_2.bed.gz
# gunzip remap2018_nr_macs2_hg38_v1_2.bed.gz
#remap = rtracklayer::import("~/projects/resources/remap/remap2018_remap2018_nr_macs2_hg38_v1_2.bed.gz")

# 1b) Get ReMap All peak set
# wget http://tagc.univ-mrs.fr/remap/download/remap2018/hg38/MACS/remap2018_all_macs2_hg38_v1_2.bed.gz
# zcat remap2018_all_macs2_hg38_v1_2.bed.gz | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > remap2018_all_macs2_hg38_v1_2.bed6.bed
remap = rtracklayer::import("~/projects/resources/remap/remap2018_all_macs2_hg38_v1_2.bed6.bed.gz")
trans_factors = unique(remap$name)
trans_factors = trans_factors [ order( trans_factors )]

# 2) get F5 promoter region
data("promoters", package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500

data("enhancers", package = "xcore")
enhancers_ext500 = enhancers
GenomicRanges::start(enhancers_ext500)   = GenomicRanges::start (enhancers)   - 500
GenomicRanges::end  (enhancers_ext500)   = GenomicRanges::end   (enhancers)   + 500

rm( enhancers, promoters) ; gc()

# 3) create intersection matrix

overlap_mat_promoters = matrix (0L, nrow = length(promoters_ext500), ncol = length(trans_factors) )
colnames(overlap_mat_promoters) = trans_factors
rownames(overlap_mat_promoters) = promoters_ext500$name

overlap_mat_enhancers = matrix (0L, nrow = length(enhancers_ext500), ncol = length(trans_factors) )
colnames(overlap_mat_enhancers) = trans_factors
rownames(overlap_mat_enhancers) = enhancers_ext500$name


for ( j in 1:length(trans_factors)){

  query_tf = remap[remap$name == trans_factors [j]]

  hits = GenomicRanges::findOverlaps(promoters_ext500 ,query_tf)
  overlap_mat_promoters [ hits@from , j] =  1L

  hits = GenomicRanges::findOverlaps(enhancers_ext500 ,query_tf)
  overlap_mat_enhancers [ hits@from , j] =  1L

  print( j)
}

rm(hits, query_tf, j, remap, export_dir) ; gc()

save( overlap_mat_enhancers, file = "data-raw/remap/overlap_mat_enhancers.rda")
save( overlap_mat_promoters, file = "data-raw/remap/overlap_mat_promoters.rda")

# 4) Gene level summarization
empty = promoters_ext500$SYMBOL == ""
promoters_ext500_for_gene_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_gene_collapsing = overlap_mat_promoters[ !empty ,]

empty = promoters_ext500$ENTREZID == ""
promoters_ext500_for_entrez_collapsing = promoters_ext500 [ !empty]
overlap_mat_promoters_for_entrez_collapsing = overlap_mat_promoters[ !empty ,]

rm(overlap_mat_promoters, overlap_mat_enhancers ); gc()

overlap_mat_symbol = rowsum(overlap_mat_promoters_for_gene_collapsing , group = promoters_ext500_for_gene_collapsing$SYMBOL)
overlap_mat_symbol[ overlap_mat_symbol2 > 0 ] =1

overlap_mat_entrez = rowsum(overlap_mat_promoters_for_entrez_collapsing , group = promoters_ext500_for_entrez_collapsing$ENTREZID )
overlap_mat_entrez[ overlap_mat_entrez > 0 ] =1

# 5) Exporting
save( overlap_mat_symbol, file = "data-raw/remap/overlap_mat_symbol.rda")
save( overlap_mat_entrez, file = "data-raw/remap/overlap_mat_entrez.rda")

remap_entrez = overlap_mat_entrez
remap_symbol = overlap_mat_symbol
usethis::use_data( remap_entrez , internal = FALSE, overwrite = TRUE)
usethis::use_data( remap_symbol , internal = FALSE, overwrite = TRUE)
