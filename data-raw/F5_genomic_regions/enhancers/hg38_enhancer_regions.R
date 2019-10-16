# Bogumil Kaczkowski August 5, 2019
# The purpose of the script is to combine FANTOM5 hg19 and hg38 enhnacer
# and to create a non redundant set
# Additionally the enhancers are overlapped with ROADMAP enhancer, promoter and dyadic regios
# as well as with EP300 peaks.

library(GenomicRanges)
#liftOver chain
hg19_hg38_chain = rtracklayer::import.chain("~/projects/resources/liftOver_chains/hg19ToHg38.over.chain")


# HG38 ENHANCERS (from reprocessed/remapped data)
hg38_enh   = rtracklayer::import.bed("~/projects/resources/FANTOM5/hg38/F5.hg38.enhancers.bed.gz" )
hg38_enh$itemRgb = NULL ; hg38_enh$thick = NULL ; hg38_enh$blocks = NULL; hg38_enh$score = NULL;


# ORIGINAL hg19 FANTOM5 enhancers (based on mapping to hg19)
hg19_enh = rtracklayer::import.bed("~/projects/resources/FANTOM5/human_permissive_enhancers_phase_1_and_2.bed.gz" )
hg19_enh$itemRgb = NULL ; hg19_enh$thick = NULL ; hg19_enh$blocks = NULL; hg19_enh$score = NULL;
# lifting over to hg38
hg19_enh = rtracklayer::liftOver(hg19_enh, chain = hg19_hg38_chain )

#making a copy for later use (to find overlaps/annotate with final enhancer set)
hg19_enh_copy = unlist(hg19_enh)
hg19_enh_copy$name2 = paste( "hg19_" , hg19_enh_copy$name , sep = "")

#removing hg19 enhancers that overlap the hg38 enhancers
hits = findOverlaps(hg19_enh ,hg38_enh , maxgap = 5)
hg19_enh = hg19_enh [- hits@from , ]

# we unlist, reduce the overlapping and close by enhancers, give them name, sort them by length and if an enhancer was lifted over to
# two different regions we remove the shorter one (to have one to one mapping)
hg19_enh = unlist(hg19_enh)
hg19_enh_reduced = reduce(hg19_enh, min.gapwidth=25 )
hg19_enh_reduced = hg19_enh_reduced [ order( width(hg19_enh_reduced), decreasing = TRUE )]
hits = findOverlaps(hg19_enh_reduced , hg19_enh )
hg19_enh_reduced$name = ""
hg19_enh_reduced$name [ hits@from  ] = hg19_enh$name [ hits@to ]
hg19_enh_reduced = hg19_enh_reduced [ ! duplicated(hg19_enh_reduced$name)]
hg19_enh_reduced$name = paste( "hg19_" , hg19_enh_reduced$name , sep = "")


# merging hg38 enhancers with lifted-over and adjusted hg19-only enhancer
all_enh = c(hg38_enh ,hg19_enh_reduced )
all_enh <- sortSeqlevels(all_enh)
all_enh <- sort(all_enh)

seqinfo_hg38  =  rtracklayer::SeqinfoForUCSCGenome("hg38")
seqinfo_hg38 = seqinfo_hg38[
  names(seqinfo_hg38) [names(seqinfo_hg38) %in% seqnames(all_enh)]
  ]
seqinfo(all_enh) = seqinfo_hg38

#.... hg19 enhancers...
all_enh$hg19_enhancer = ""
hits = findOverlaps(all_enh ,hg19_enh_copy )
all_enh$hg19_enhancer[hits@from] = hg19_enh_copy$name[hits@to]

#.... hg38 enhancers...
all_enh$hg38_enhancer = ""
hits = findOverlaps(all_enh ,hg38_enh )
all_enh$hg38_enhancer [hits@from] = hg38_enh$name[hits@to]

#.... EP300 peaks...
all_enh$ep300 = ""
ep300  = rtracklayer::import.bed("~/projects/resources/remap/remap2018_EP300_nr_macs2_hg38_v1_2.bed.gz" )
hits = findOverlaps(all_enh ,ep300)
all_enh$ep300[hits@from] = "EP300"
all_enh$ep300 = as.factor(all_enh$ep300)

# ROADMAP
roadmap_hg19 = rtracklayer::import.bed("~/projects/resources/roadmap/EpigenomeRoadmapDHS_bed6.bed", colnames = c("chrom", "start", "end", "name", "score"))
roadmap = rtracklayer::liftOver(roadmap_hg19, chain = hg19_hg38_chain )
roadmap = unlist(roadmap)

all_enh$roadmap = ""
hits = findOverlaps(all_enh , roadmap )
all_enh$roadmap[hits@from] = roadmap$name[hits@to]
all_enh$roadmap_type    =  as.factor( sub("\\|.*" , "" , all_enh$roadmap) )
rm(roadmap_hg19, roadmap ); gc()

dfam  = rtracklayer::import.bed("data-raw/repeats/hg38_dfam.3.1.nrph.hits.bed.gz" )
all_enh$repeat_dfam = ""
hits = findOverlaps(all_enh , dfam )
all_enh$repeat_dfam [hits@from] = dfam$name[hits@to]

all_enh$genome_call = ""
all_enh$genome_call[ all_enh$hg19_enhancer != "" & all_enh$hg38_enhancer != ""] = "hg19&hg38"
all_enh$genome_call[ all_enh$hg19_enhancer != "" & all_enh$hg38_enhancer == ""] = "hg19"
all_enh$genome_call[ all_enh$hg19_enhancer == "" & all_enh$hg38_enhancer != ""] = "hg38"

#EXPORTING the bed file
rtracklayer::export.bed(object = all_enh ,
                        con = "data-raw/F5_genomic_regions/enhancers/FANTOM5_hg38_and_hg19_enhancers_2019_10_05.bed",
                        format = "bed")


save(all_enh, file ="data-raw/F5_genomic_regions/enhancers/FANTOM5_hg38_and_hg19_enhancers_2019_10_05_GRanges.RData")

enhancers = all_enh
usethis::use_data( enhancers , internal = FALSE, overwrite = TRUE)
