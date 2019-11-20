#FANTOM5 PROMOTERS
data(promoters, package = "xcore")
promoters_1kb    = promoters
GenomicRanges::start(promoters_1kb)    = IRanges::mid(promoters) -500
GenomicRanges::end  (promoters_1kb)    = IRanges::mid(promoters) +499
promoters_1kb = GenomicRanges::trim(promoters_1kb)

phastCons_file = "~/projects/resources/conservation/hg38.phastCons100way.bw"
dpi_cons_mat = matrix(0, nrow =length(promoters_1kb), ncol = 1000 )
rownames(dpi_cons_mat) = promoters_1kb$name

for ( i in 1: length(promoters_1kb)){
  try(expr = (dpi_cons_mat[ i , ] = (rtracklayer::import.bw(phastCons_file ,selection = as(promoters_1kb[i], "BigWigSelection")))$score ), silent = TRUE)
  print(i)
}
save( dpi_cons_mat, file = "data-raw/conservation/dpi_cons_mat.rda")

#FANTOM5 ENHANCERS
data(enhancers, package = "xcore")
enhancers_1kb    = enhancers
GenomicRanges::start(enhancers_1kb)    = IRanges::mid(enhancers) -500
GenomicRanges::end  (enhancers_1kb)    = IRanges::mid(enhancers) +499
enhancers_1kb = GenomicRanges::trim(enhancers_1kb)

enh_cons_mat = matrix(0, nrow =length(enhancers_1kb), ncol = 1000 )
rownames(enh_cons_mat) = enhancers_1kb$name
phastCons_file = "~/projects/resources/conservation/hg38.phastCons100way.bw"

for ( i in 1: length(enhancers_1kb)){
  try(expr = (enh_cons_mat[ i , ] = (rtracklayer::import.bw(phastCons_file ,selection = as(enhancers_1kb[i], "BigWigSelection")))$score ), silent = TRUE)
  print(i)
}
save( enh_cons_mat, file = "data-raw/conservation/enh_cons_mat.rda")

load( "data-raw/conservation/enh_cons_mat.rda" )
load( "data-raw/conservation/dpi_cons_mat.rda" )

bins = sort( rep( seq(1, 100, 1), 10) )

conservation_promoters_1kb = t(rowsum( x = t(dpi_cons_mat), group = bins)) / 10
conservation_enhancers_1kb = t(rowsum( x = t(enh_cons_mat), group = bins)) / 10

plot ( colMeans(conservation_promoters_1kb) , ylim = c(0,.35), type = "l", xaxt = 'n')
axis(side=1,at=seq( 0, 100, 10),labels=seq( 0, 1000, 100) - 500 )

plot ( colMeans(conservation_enhancers_1kb) , ylim = c(0,.35), type = "l", xaxt = 'n')
axis(side=1,at=seq( 0, 100, 10),labels=seq( 0, 1000, 100) - 500 )

plot( x = 1:1000 - 500 , y = colMeans(dpi_cons_mat),
      xlab = "position", col = "black",
      ylab = "Mean Conservation", type = "l")
