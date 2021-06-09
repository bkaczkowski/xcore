
jaspar = PWMEnrich::readMotifs("~/projects/resources/PWMs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

data(promoters, package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500
promoters_ext500 = GenomicRanges::trim(promoters_ext500)
promoters_ext500 = GenomicRanges::reduce(promoters_ext500)

library(BSgenome.Hsapiens.UCSC.hg38)
promoters_ext500_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , promoters_ext500)

afreq  = alphabetFrequency(promoters_ext500_seq)
afreq  = afreq/rowSums(afreq)
promoters_ext500_seq = promoters_ext500_seq [ afreq[,"N"] <0.1  ]

bg_jaspar_dpi_reduced_10 = PWMEnrich::makeBackground(motifs = jaspar , bg.seq = promoters_ext500_seq[ seq( 1, 96088, 10)] ,algorithm="human" )

save( bg_jaspar_dpi_reduced_10 , file = "data-raw/motif_enrichment/bg_jaspar_dpi_reduced_10.rda")

#bg_jaspar_dpi_reduced_test = PWMEnrich::makeBackground(motifs = jaspar , bg.seq = promoters_ext500_seq[ seq( 1,96088, 100)] ,algorithm="human" )

#save( bg_jaspar_dpi_reduced_test , file = "data-raw/motif_enrichment/bg_jaspar_dpi_reduced_test.rda")

# nohup Rscript data-raw/motif_enrichment/F5_dpi_reduced_background_jaspar.R > data-raw/motif_enrichment/F5_dpi_reduced_background_jaspar.R.nohup.out &


