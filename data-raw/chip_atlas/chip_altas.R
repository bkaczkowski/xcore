
data("promoters_detailed", package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500

peak_call_files = list.files("~/projects/resources/chip_atlas/split_by_experiment/", full.names = T, pattern = ".bed.gz$")
names(peak_call_files) = gsub( ".*/|.bed.gz" , "", peak_call_files )

create_intersection_wrapper = function(peak_call_file,  pre_defined_regions){
  gr = rtracklayer::import(peak_call_file)
  length(gr)
  int = xcore::create_intersection_vector( query = gr, pre_defined_regions = pre_defined_regions )
  return( list( len = length(gr) ,  int = int) )
}

dpi_chip_list = list()
for ( i in 1:length(peak_call_files)){
  peak_call_file = peak_call_files[i]
  gr = rtracklayer::import(peak_call_file)
  int = xcore::create_intersection_vector( query = gr, pre_defined_regions = promoters_ext500 )
  dpi_chip_list[[i]] = list( len = length(gr) ,  int = int)
  print(i)
}

extract_lenght = function(x) {x$len}
number_of_peaks = unlist(lapply( dpi_chip_list, extract_lenght))
hist(log10(number_of_peaks))

extract_int = function(x) {x$int}
intersection_list = lapply( dpi_chip_list, extract_int)
intersections = do.call(cbind, intersection_list)
colnames(intersections) = names(peak_call_files)
rownames(intersections) = promoters$name
rm(dpi_chip_list,intersection_list ) ; gc()

#f5dpi_n_chip_atlas = intersections
#save(f5dpi_n_chip_atlas, file = "data-raw/chip_atlas/f5dpi_n_chip_atlas.RData" )
#rm(f5dpi_n_chip_atlas ) ; gc()


peaks_overlapping_DPIs = colSums(intersections)
dev.off()
pdf("data-raw/chip_atlas/NumberOfPeaks_vs_PeaksOverlappingDPIs.pdf", width = 8, height = 6)
smoothScatter(number_of_peaks /1000 , peaks_overlapping_DPIs/1000 ,  xlim = c(0,110) )
dev.off()

intersections_gene_level = rowsum(intersections, group = promoters$SYMBOL_F5_annot)
not_annotated = intersections_gene_level[1,]
intersections_gene_level = intersections_gene_level[-1,]
intersections_gene_level[ intersections_gene_level > 0 ] =1
f5genes_n_chip_atlas = intersections_gene_level
save(f5genes_n_chip_atlas, file = "data-raw/chip_atlas/f5genes_n_chip_atlas.RData" )
rm(f5genes_n_chip_atlas ) ; gc()

intersections_gene_level = rowsum(intersections, group = promoters$SYMBOL_gencode)
not_annotated = intersections_gene_level[1,]
intersections_gene_level = intersections_gene_level[-1,]
intersections_gene_level[ intersections_gene_level > 0 ] =1
gencode_n_chip_atlas = intersections_gene_level
save(gencode_n_chip_atlas, file = "data-raw/chip_atlas/gencode_n_chip_atlas.RData" )
rm(gencode_n_chip_atlas ) ; gc()


peaks_overlapping_genes = colSums(intersections_gene_level)
dev.off()
pdf("data-raw/chip_atlas/NumberOfPeaks_vs_PeaksOverlappingGencodeGenes.pdf", width = 8, height = 6)
smoothScatter(number_of_peaks /1000 , peaks_overlapping_genes/1000 ,  xlim = c(0,110) )
dev.off()


peaks_at_genes = colSums(intersections_gene_level > 0)
plot(number_of_peaks /1000 , peaks_at_genes/1000 ,  xlim = c(0,110) , pch = ".")

intersections_gene_level_sel = intersections_gene_level [ rowSums(intersections_gene_level)> 50,]
intersections_gene_level_sel = intersections_gene_level [ , peaks_at_genes > 100 ]

dev.off()
pdf("data-raw/chip_atlas/upset_top30_freq_chip.pdf", width = 15, height = 10)
UpSetR::upset(data.frame(intersections_gene_level_sel) , nsets = 30, order.by = "freq")
dev.off()

dev.off()
pdf("data-raw/chip_atlas/upset_top30_freq_genes.pdf", width = 15, height = 10)
UpSetR::upset(data.frame(t(intersections_gene_level_sel) ), nsets = 30, order.by = "freq")
dev.off()


