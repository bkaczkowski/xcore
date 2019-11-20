library(BSgenome.Hsapiens.UCSC.hg38)

#FANTOM5 PROMOTERS
data(promoters, package = "xcore")
promoters_1.1kb    = promoters
GenomicRanges::start(promoters_1.1kb)    = IRanges::mid(promoters) -550
GenomicRanges::end(promoters_1.1kb)      = IRanges::mid(promoters) +549
promoters_1.1kb = GenomicRanges::trim(promoters_1.1kb)
promoters_1.1kb_seq    = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , promoters_1.1kb)
names(promoters_1.1kb_seq) = promoters_1.1kb$name

Biostrings::writeXStringSet(promoters_1.1kb_seq, "data-raw/genomic_regions_fasta/promoters_1.1kb_seq.fasta")
save(promoters_1.1kb_seq    , file = "data-raw/genomic_regions_fasta/promoters_1.1kb_seq.rda" )
#usethis::use_data( promoters_1.1kb_seq   , internal = FALSE, overwrite = TRUE)

#FANTOM5 ENHANCERS
data(enhancers, package = "xcore")
enhancers_1.1kb    = enhancers
GenomicRanges::start(enhancers_1.1kb)    = IRanges::mid(enhancers) -550
GenomicRanges::end(enhancers_1.1kb)      = IRanges::mid(enhancers) +549
enhancers_1.1kb = GenomicRanges::trim(enhancers_1.1kb)

enhancers_1.1kb_seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38 , enhancers_1.1kb)
table( width(enhancers_1.1kb_seq) == width( enhancers_1.1kb)) # checking if the lenght of the sequences is the same (same order)
names(enhancers_1.1kb_seq) = enhancers_1.1kb$name
Biostrings::writeXStringSet(enhancers_1.1kb_seq, "data-raw/genomic_regions_fasta/enhancers_1.1kb_seq.fasta")

save(enhancers_1.1kb_seq , file = "data-raw/genomic_regions_fasta/enhancers_1.1kb_seq.rda" )
