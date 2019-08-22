#system("mkdir miseq_peaks_calls")
#system("scp dgt-ac4:/home/bogumil/brd_chipseq/* miseq_peaks_calls/")

library( xcore )
peak_call_files = list.files("miseq_peaks_calls", full.names = T, pattern = "ext300_peaks.narrowPeak$")
names(peak_call_files) = gsub( ".*/|.q20.*" , "", peak_call_files )

peak_calls = list(  )
for ( i in 1:length(peak_call_files )) { 
    peak_calls [[ names(peak_call_files)[i]  ]] = 
        rtracklayer::import( con =  peak_call_files[i],
                             format = "narrowPeak") }

peak_calls_merged = merge_list_of_granges( peak_calls )
peak_calls_reduced = GenomicRanges::reduce( peak_calls_merged)
peak_calls_reduced$name = name_from_granges( peak_calls_reduced )
peak_intersections = create_intersection_matrix( list_of_granges = peak_calls, pre_defined_regions = peak_calls_reduced)
# alternative
#peak_intersections = create_intersection_matrix( list_of_granges = peak_calls , pre_defined_regions = NULL)

keep = rowSums(peak_intersections) > 1

peak_calls_reduced = peak_calls_reduced[keep]
peak_intersections = peak_intersections[keep,]

par( mar = c(12,5,2,1))
barplot( colSums(peak_intersections), las = 2, main = "number of called peaks")

gencode = rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
peak_calls_reduced = xcore::gencode_strandless_annotator(regions = peak_calls_reduced ,  gencode = gencode, output_option = 1)
peak_calls_reduced$annotatation =  xcore::simplify_gencode_annotation(peak_calls_reduced$strandless_annotatation)

roadmap = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_hg38.bed", colnames = c("chrom", "start", "end", "name", "score", "strand"))
peak_calls_reduced$roadmap = ""
hits = GenomicRanges::findOverlaps(peak_calls_reduced , roadmap )
peak_calls_reduced$roadmap[hits@from] = roadmap$name[hits@to]
peak_calls_reduced$roadmap    =  as.factor( sub("\\|.*" , "" , peak_calls_reduced$roadmap) )

data("promoters", package = "xcore")
peak_calls_reduced$DPI = 0
hits = GenomicRanges::findOverlaps(peak_calls_reduced ,promoters)
peak_calls_reduced$DPI[hits@from] = 1

data("enhancers", package = "xcore")
peak_calls_reduced$f5_enhancer = 0
hits = GenomicRanges::findOverlaps(peak_calls_reduced ,enhancers)
peak_calls_reduced$f5_enhancer[hits@from] = 1

dfam  = rtracklayer::import.bed("~/r/xcore/data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed.gz")
peak_calls_reduced$repeat_dfam = 0
hits = GenomicRanges::findOverlaps(peak_calls_reduced , dfam )
peak_calls_reduced$repeat_dfam[hits@from] = 1

annotation_mat = as.matrix(model.matrix( ~ 0  +  peak_calls_reduced$annotatation ))
colnames(annotation_mat)  = sub( "peak_calls_reduced\\$annotatation", "gencode_", colnames(annotation_mat) )
roadmap_mat = as.matrix(model.matrix( ~ peak_calls_reduced$roadmap - 1))
colnames(roadmap_mat)  = sub( "peak_calls_reduced\\$roadmap", "roadmap_", colnames(roadmap_mat) )
roadmap_mat = roadmap_mat[,-1]

intersections = data.frame(cbind( peak_intersections , annotation_mat,roadmap_mat , GenomicRanges::mcols(peak_calls_reduced) [ , 7:9]))

dpi = intersections$DPI == 1
non_repeat = intersections$repeat_dfam == 0

intersections_chosen = intersections [ dpi & non_repeat ,]
chosen_regions = peak_calls_reduced[ dpi & non_repeat]
library(GenomicRanges)
mcols(chosen_regions) = cbind( mcols(chosen_regions)  , 
                               intersections_chosen)
chosen_regions$times_called = rowSums(intersections_chosen [,1:11])

chosen_regions = chosen_regions[ order(chosen_regions$times_called , decreasing = T )]
library(xlsx)
write.xlsx(mcols(chosen_regions), file = "BRD2_BRD4_reproduced_peaks.xlsx")

library(UpSetR)
dev.off()
pdf( file = "BRD2_BRD4_reproduced_peaks_UPSET_plot.pdf" , width = 10, height = 6)
upset(intersections_chosen, nsets = 13,order.by = "freq" , nintersects = 20)
dev.off()
