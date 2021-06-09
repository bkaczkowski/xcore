#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

data(promoters, package = "xcore")
promoters_ext500 = promoters
GenomicRanges::start(promoters_ext500)   = GenomicRanges::start (promoters) - 500
GenomicRanges::end  (promoters_ext500)   = GenomicRanges::end (promoters)   + 500
promoters_ext500 = GenomicRanges::trim(promoters_ext500)
promoters_ext500_reduced = GenomicRanges::reduce(promoters_ext500)

promoters_ext500_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , promoters_ext500)
table( width(promoters_ext500_seq) == width( promoters_ext500)) # checking if the lenght of the sequences is the same (same order)
names(promoters_ext500_seq) = promoters_ext500$name

promoters_ext500_reduced = GenomicRanges::reduce(promoters_ext500)
promoters_ext500_reduced_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , promoters_ext500_reduced )

afreq  = alphabetFrequency(promoters_ext500_reduced_seq)
afreq  = afreq/rowSums(afreq)
promoters_ext500_reduced_seq = promoters_ext500_reduced_seq [ afreq[,"N"] <0.1  ]

afreq  = alphabetFrequency(promoters_ext500_seq)
afreq  = afreq/rowSums(afreq)
promoters_ext500_seq = promoters_ext500_seq [ afreq[,"N"] <0.1  ]
promoters_ext500_reduced_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , promoters_ext500_reduced )



jaspar = PWMEnrich::readMotifs("~/projects/resources/PWMs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

bg_jaspar_dpi_reduced_test = PWMEnrich::makeBackground(motifs = jaspar , bg.seq = promoters_ext500_reduced_seq [seq(1,105240, 100 )] ,algorithm="human" )



data(enhancers, package = "xcore")
enhancers_ext500 = enhancers
GenomicRanges::start(enhancers_ext500)   = GenomicRanges::start (enhancers)   - 500
GenomicRanges::end  (enhancers_ext500)   = GenomicRanges::end   (enhancers)   + 500

enhancers_ext500_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , enhancers_ext500)
table( width(enhancers_ext500_seq) == width( enhancers_ext500)) # checking if the lenght of the sequences is the same (same order)
names(enhancers_ext500_seq) = enhancers_ext500$name
Biostrings::writeXStringSet(enhancers_ext500_seq, "data-raw/renomic_regions_fasta/enhancers_ext500_seq.fasta")

save(enhancers_ext500_seq , file = "data-raw/renomic_regions_fasta/enhancers_ext500_seq.rda" )
save(promoters_ext500_seq , file = "data-raw/renomic_regions_fasta/promoters_ext500_seq.rda" )

usethis::use_data( enhancers_ext500_seq   , internal = FALSE, overwrite = TRUE)
usethis::use_data( promoters_ext500_seq   , internal = FALSE, overwrite = TRUE)

