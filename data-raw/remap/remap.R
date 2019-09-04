# Bogumil Kaczkowski
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions with ReMap Chip Seq data-base

# 1a) Get ReMap non redundant peaks
# wget http://tagc.univ-mrs.fr/remap/download/remap2018/hg38/MACS/remap2018_nr_macs2_hg38_v1_2.bed.gz
# gunzip remap2018_nr_macs2_hg38_v1_2.bed.gz
#remap = rtracklayer::import("~/projects/resources/remap/remap2018_remap2018_nr_macs2_hg38_v1_2.bed.gz")

# 1b) Get ReMap All peak set
# wget http://tagc.univ-mrs.fr/remap/download/remap2018/hg38/MACS/remap2018_all_macs2_hg38_v1_2.bed.gz
# zcat remap2018_all_macs2_hg38_v1_2.bed.gz | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > remap2018_all_macs2_hg38_v1_2.bed6.bed
remap = rtracklayer::import("~/projects/resources/remap/remap2018_all_macs2_hg38_v1_2.bed6.bed")
trans_factors = unique(remap$name)

# 2) get F5 promoter region
data("promoters", package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500

# 3) create intersection matrix

intersection_matrix = matrix (0L, nrow = length(promoters), ncol = length( trans_factors ) )
colnames(intersection_matrix) = trans_factors
rownames(intersection_matrix) = promoters_ext500$name

for ( j in 1:length(trans_factors)){

  query_tf = remap[remap$name == trans_factors [j]]

  hits = GenomicRanges::findOverlaps(promoters_ext500 ,query_tf)

  intersection_matrix[ hits@from , j] =  1L
  print( j)
}

f5_n_remap_all = intersection_matrix
usethis::use_data(f5_n_remap_all,overwrite = TRUE)

# for the non redundant set
#f5_n_remap = intersection_matrix
#usethis::use_data(f5_n_remap,overwrite = TRUE)
