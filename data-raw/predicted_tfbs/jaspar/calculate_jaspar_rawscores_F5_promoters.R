library(GenomicRanges)
library(PWMEnrich)
library( parallel)
registerCoresPWMEnrich(20)

load("promoters_1.1kb_seq.rda")
jaspar = PWMEnrich::readMotifs("JASPAR2020_CORE_vertebrates_redundant_pfms_jaspar.txt")
names(jaspar) = sub("\t", "_" ,names(jaspar) )
names(jaspar) = sub("::", "-" ,names(jaspar) )
names(jaspar) = sub("\\)$", "" ,names(jaspar) )
names(jaspar) = sub("\\(", "." ,names(jaspar) )

files = list.files("/work/bogumil/jaspar_raw_scores_promoters", full.names = T)
done  = gsub( ".*/|.csv.gz$|.csv$", "", files)

table( names(jaspar) %in% done)

jaspar = jaspar [!names(jaspar) %in% done ]


#load("data-raw/genomic_regions_fasta/promoters_1.1kb_seq.rda")
#jaspar = PWMEnrich::readMotifs("~/projects/resources/PWMs/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

promoters_1.1kb_seq = promoters_1.1kb_seq [ - grep( "chrM", names(promoters_1.1kb_seq))]

#checking for the lenght to be the same
min(width(promoters_1.1kb_seq)) == max(width(promoters_1.1kb_seq))

bins = sort( rep( seq(1, 220, 1), 10) )
collapse_bins = function( x ) { tapply( x, INDEX =bins, FUN = max , na.rm = TRUE)}

for ( i in 1:length(jaspar)){
  jaspar_scores = PWMEnrich::motifScores( sequences = promoters_1.1kb_seq ,motifs = jaspar[i],  raw.scores = TRUE )
  jaspar_scores_binned = mclapply( jaspar_scores, FUN = collapse_bins, mc.cores = 20)
  rm(jaspar_scores) ; gc()
  jaspar_mat = t(do.call(cbind , jaspar_scores_binned))
  jaspar_mat = jaspar_mat[ , c(6:105,116:215) ]
  write.csv( jaspar_mat , file = paste( "/work/bogumil/jaspar_raw_scores_promoters/", names(jaspar)[i], ".csv", sep = ""))
}

