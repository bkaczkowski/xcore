library(GenomicRanges)
library(PWMEnrich)
registerCoresPWMEnrich(30)

load("promoters_ext500_seq.rda")
load("regulon.rda")

jaspar = PWMEnrich::readMotifs("JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

promoters_regulon = PWMEnrich::motifScores( sequences = promoters_ext500_seq ,motifs = regulon)

save( promoters_regulon, file = "promoters_regulon.rda")

