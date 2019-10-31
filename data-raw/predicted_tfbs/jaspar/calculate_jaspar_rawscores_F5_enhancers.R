library(GenomicRanges)
library(PWMEnrich)
library( parallel)
registerCoresPWMEnrich(10)
#system("mkdir jaspar_raw_scores")

load("enhancers_1.1kb_seq.rda")
jaspar = PWMEnrich::readMotifs("JASPAR2020_CORE_vertebrates_redundant_pfms_jaspar.txt")
names(jaspar) = sub("\t", "_" ,names(jaspar) )
names(jaspar) = sub("::", "-" ,names(jaspar) )
names(jaspar) = sub("\\)$", "" ,names(jaspar) )
names(jaspar) = sub("\\(", "." ,names(jaspar) )

files = list.files("/work/bogumil/jaspar_raw_scores_enhancers", full.names = T)
done  = gsub( ".*/|.csv.gz$|.csv$", "", files)

table( names(jaspar) %in% done)

jaspar = jaspar [!names(jaspar) %in% done ]

#checking for the lenght to be the same
min(width(enhancers_1.1kb_seq)) == max(min(width(enhancers_1.1kb_seq)))

bins = sort( rep( seq(1, 220, 1), 10) )
collapse_bins = function( x ) { tapply( x, INDEX =bins, FUN = max , na.rm = TRUE)}

for ( i in 1:length(jaspar)){
  enhancers_jaspar_scores = PWMEnrich::motifScores( sequences = enhancers_1.1kb_seq ,motifs = jaspar[i],  raw.scores = TRUE )
  enhancers_jaspar_scores_binned = mclapply( enhancers_jaspar_scores, FUN = collapse_bins, mc.cores = 10)
  rm(enhancers_jaspar_scores) ; gc()
  enhancers_jaspar_mat = t(do.call(cbind , enhancers_jaspar_scores_binned))
  enhancers_jaspar_mat = enhancers_jaspar_mat[ , c(6:105,116:215) ]
  write.csv( enhancers_jaspar_mat , file = paste( "/work/bogumil/jaspar_raw_scores_enhancers/", names(jaspar)[i], ".csv", sep = ""))
}
