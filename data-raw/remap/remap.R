# Bogumil Kaczkowski August 27, 2019
# the purpose of the script is to overlap the FANTOM5 DPI/promoter regions with ReMap Chip Seq data-base

# 1) Get ReMap non redundant peaks
# wget http://tagc.univ-mrs.fr/remap/download/remap2018/hg38/MACS/remap2018_nr_macs2_hg38_v1_2.bed.gz
# gunzip remap2018_nr_macs2_hg38_v1_2.bed.gz

# to get the non redundant set
#remap = rtracklayer::import("~/projects/resources/remap/remap2018_remap2018_nr_macs2_hg38_v1_2.bed.gz")

# to get all peak set
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

UpSetR::upset(data.frame(intersection_matrix), order.by = "freq")

plot(sort(rowSums(intersection_matrix)))
barplot(sort(colSums(intersection_matrix)), las = 2)



UpSetR::upset(data.frame(f5_n_remap), order.by = "freq")

rownames (f5_n_remap) = promoters$name

# ALTERNATIVE WAY WITH THE SAME RESULT
# mkdir -p remap2018_nr_macs2_hg38_v1_2_split_by_tf
# for tf in `cut -f 4 remap2018_nr_macs2_hg38_v1_2.bed | sort | uniq`; do
# grep -w $tf remap2018_nr_macs2_hg38_v1_2.bed | awk 'OFS="\t" {print $1,$2,$3}' > remap2018_nr_macs2_hg38_v1_2_split_by_tf/$tf.bed
# done

peak_call_files = list.files("~/projects/resources/remap/remap2018_nr_macs2_hg38_v1_2_split_by_tf/", full.names = T, pattern = ".bed.gz$")
names(peak_call_files) = gsub( ".*/|.bed.gz" , "", peak_call_files )

create_intersection_wrapper = function(peak_call_file,  pre_defined_regions){
  gr = rtracklayer::import(peak_call_file)
  xcore::create_intersection_vector( query = gr, pre_defined_regions = pre_defined_regions )
}

gr = rtracklayer::import(peak_call_files[3])


f5_n_remap = lapply(peak_call_files, create_intersection_wrapper, promoters_ext500 )
f5_n_remap = do.call(cbind, f5_n_remap)
colnames(f5_n_remap) = names(peak_call_files)
rownames(f5_n_remap) = promoters$name
UpSetR::upset(data.frame(f5_n_remap), order.by = "freq")

#usethis::use_data(f5_n_remap)

#f5_n_remap_drop0 = Matrix::drop0 (f5_n_remap)
#f5_n_remap_logical = ifelse(f5_n_remap == 1, yes = TRUE , no = FALSE)

remap_all = data.table::fread("~/projects/resources/remap/remap2018_all_macs2_hg38_v1_2.bed")
remap_all$V9 = NULL
remap_all$V8 = NULL
remap_all$V7 = NULL
remap_all$cell = sub( ".*\\." , "", remap_all$V4)
remap_all$exp  = sub( "\\..*" , "", remap_all$V4)

remap_all$tf = do.call( rbind, strsplit(remap_all$V4 , split = ".", fixed = TRUE  ) )[,2]
table(remap_all$cell)

sort(table(remap_all$cell))

remap_all_unique = remap_all [ !duplicated(remap_all$V4) ,]


