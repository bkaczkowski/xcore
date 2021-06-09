
gencode = rtracklayer::import.gff(con = "~/projects/resources/gencode/gencode31/gencode.v31.annotation.gff3.gz")
gencode_promoters_full   = GenomicRanges::promoters( gencode[ gencode$type == "transcript"] , upstream = 500, downstream = 500 , use.names=TRUE)
gencode_exons_full       = gencode[ gencode$type == "exon"]
gencode_gene_full        = gencode[ gencode$type == "gene"]


gencode_promoters       = GenomicRanges::reduce ( gencode_promoters_full )
gencode_exons           = GenomicRanges::reduce ( gencode_exons_full )
gencode_gene            = GenomicRanges::reduce ( gencode_gene_full )

gencode_exons           = GenomicRanges::setdiff( gencode_exons, gencode_promoters)
gencode_intron          = GenomicRanges::setdiff( gencode_gene , gencode_promoters )
gencode_intron          = GenomicRanges::setdiff( gencode_intron , gencode_exons )

length( GenomicRanges::findOverlaps(gencode_promoters_full, gencode_promoters) )
length( GenomicRanges::findOverlaps(gencode_exons_full, gencode_exons) )
length( GenomicRanges::findOverlaps(gencode_gene_full, gencode_gene) )

length( GenomicRanges::findOverlaps(gencode_promoters, gencode_exons) )
length( GenomicRanges::findOverlaps(gencode_exons, gencode_intron) )
length( GenomicRanges::findOverlaps(gencode_promoters, gencode_intron) )

gencode_promoters$name = "promoter"
gencode_exons$name = "exon"
gencode_intron$name = "intron"

gencode_simple = c(gencode_promoters, gencode_exons, gencode_intron)

usethis::use_data( gencode_simple , internal = FALSE, overwrite = TRUE)



