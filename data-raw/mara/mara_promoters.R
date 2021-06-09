# Bogumil Kaczkowski
# the purpose of the script is to obtain a TFBS (Transcription Factor Binding Site ) table
# to perform MARA and MARA-like analysis


##### PROMOTER LEVEL #####
load("data/promoters.rda")
promoters = promoters [GenomicRanges::seqnames(promoters) != "chrM"]
rtracklayer::export(promoters, "data-raw/mara/promoters.bed")

# shell script below
"
cd ~/r/xcore/data-raw/mara/
ln ~/projects/resources/mara/hg38.sites_no_conservation.bed.bz2 hg38.sites_no_conservation.bed.bz2
ln ~/projects/resources/mara/hg38.sites.bed.bz2 hg38.sites.bed.bz2

python3 ~/python/mara/make_profile.py --output=profile.txt --input=promoters.bed --tfbs=hg38.sites.bed.bz2
python3 ~/python/mara/make_profile.py --output=profile_no_conservation.txt --input=promoters.bed --tfbs=hg38.sites_no_conservation.bed.bz2

python3 ~/python/mara/associate_tfbs.py --output=tfbs.txt --input=promoters.bed --tfbs=hg38.sites.bed.bz2 --profile=profile.txt
python3 ~/python/mara/associate_tfbs.py --output=tfbs_no_conservation.txt --input=promoters.bed --tfbs=hg38.sites_no_conservation.bed.bz2 --profile=profile_no_conservation.txt
")


##### GENE LEVEL #####
load("data/promoters.rda")
promoters = promoters [ order( promoters$score , decreasing = TRUE)]
promoters = promoters[ promoters$SYMBOL != ""]
promoters = promoters[GenomicRanges::seqnames(promoters) != "chrM"]
promoters = promoters[ !duplicated(promoters$SYMBOL)]
promoters$name = promoters$SYMBOL

rtracklayer::export(promoters, "data-raw/mara/gene_promoters.bed")

# shell script below
"
cd ~/r/xcore/data-raw/mara/
python3 ~/python/mara/make_profile.py --output=gene_profile.txt --input=gene_promoters.bed --tfbs=hg38.sites.bed.bz2
python3 ~/python/mara/make_profile.py --output=gene_profile_no_conservation.txt --input=gene_promoters.bed --tfbs=hg38.sites_no_conservation.bed.bz2

python3 ~/python/mara/associate_tfbs.py --output=gene_tfbs.txt --input=gene_promoters.bed --tfbs=hg38.sites.bed.bz2 --profile=gene_profile.txt
python3 ~/python/mara/associate_tfbs.py --output=gene_tfbs_no_conservation.txt --input=gene_promoters.bed --tfbs=hg38.sites_no_conservation.bed.bz2 --profile=gene_profile_no_conservation.txt
")

gene_tfbs = data.matrix(read.table("data-raw/mara/gene_tfbs.txt"))
mara_gene_tfbs = Matrix::drop0( gene_tfbs )
save(object = mara_gene_tfbs, file = "data-raw/mara/mara_gene_tfbs.rda", compress = "bzip2")
usethis::use_data( mara_gene_tfbs   , internal = FALSE, overwrite = TRUE)

gene_tfbs_no_conservation = data.matrix(read.table("data-raw/mara/gene_tfbs_no_conservation.txt"))
mara_gene_tfbs_no_conservation = Matrix::drop0( gene_tfbs_no_conservation )
save(object = mara_gene_tfbs_no_conservation, file = "data-raw/mara/mara_gene_tfbs_no_conservation.rda", compress = "bzip2")
usethis::use_data( mara_gene_tfbs_no_conservation   , internal = FALSE, overwrite = TRUE)

tfbs = data.matrix(read.table("data-raw/mara/tfbs.txt"))
mara_promoter_tfbs = Matrix::drop0( tfbs )
save(object = mara_promoter_tfbs, file = "data-raw/mara/mara_promoter_tfbs.rda", compress = "bzip2")
usethis::use_data( mara_promoter_tfbs   , internal = FALSE, overwrite = TRUE)

tfbs_no_conservation = data.matrix(read.table("data-raw/mara/tfbs_no_conservation.txt"))
mara_promoter_tfbs_no_conservation = Matrix::drop0( tfbs )
save(object = mara_promoter_tfbs_no_conservation, file = "data-raw/mara/mara_promoter_tfbs_no_conservation.rda", compress = "bzip2")
usethis::use_data( mara_promoter_tfbs_no_conservation   , internal = FALSE, overwrite = TRUE)


gene_profile = data.matrix(read.table("data-raw/mara/gene_profile.txt", header = TRUE, row.names = 1))
gene_profile_no_conservation = data.matrix(read.table("data-raw/mara/gene_profile_no_conservation.txt", header = TRUE, row.names = 1))
promoter_profile = data.matrix(read.table("data-raw/mara/profile.txt", header = TRUE, row.names = 1))

gene_profile = gene_profile [ rownames(gene_profile) %in% rownames(gene_profile_no_conservation),]
promoter_profile = promoter_profile [ rownames(promoter_profile) %in% rownames(gene_profile_no_conservation),]

plot( rowSums(gene_profile), rowSums(gene_profile_no_conservation) )
plot( rowSums(gene_profile), rowSums(promoter_profile) )


plot( rowSums(gene_profile), rowSums(gene_profile_no_conservation) )

gene_profile_no_conservation = xcore::order_by_variance(gene_profile_no_conservation)

i = 2
plot(gene_profile_no_conservation[i,], type = "l", main = rownames(gene_profile_no_conservation)[i],
     ylab = "P")

plot(gene_profile[i,], type = "l", main = rownames(gene_profile)[i],
     ylab = "P")

gene_profile = xcore::order_by_variance(gene_profile)
i = 1
plot(gene_profile[i,], type = "l", main = rownames(gene_profile)[i],
     ylab = "P")
hist( rowSums(gene_profile), breaks = 20, xlim = c(0,1))

promoter_profile = data.matrix(read.table("data-raw/mara/profile.txt", header = TRUE, row.names = 1))
promoter_profile = xcore::order_by_variance(promoter_profile)
plot(promoter_profile[i,], type = "l", main = rownames(gene_profile)[i],
     ylab = "P")

hist( rowSums(promoter_profile), breaks = 15, xlim = c(0,1))


enhancer_profile = data.matrix(read.table("data-raw/mara/profile_enhancers.txt", header = TRUE, row.names = 1))
enhancer_profile = xcore::order_by_variance(enhancer_profile)
plot(enhancer_profile[i,], type = "l", main = rownames(enhancer_profile)[i],
     ylab = "P")

