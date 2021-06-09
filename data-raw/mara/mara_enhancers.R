# Bogumil Kaczkowski
# the purpose of the script is to obtain a TFBS (Transcription Factor Binding Site ) table
# to perform MARA and MARA-like analysis


##### ENHANCER LEVEL #####
load("data/enhancers.rda")
GenomicRanges::strand(enhancers) = "+"
rtracklayer::export.bed(enhancers, "data-raw/mara/enhancers.bed")

# shell script below
"
cd ~/r/xcore/data-raw/mara/
ln ~/projects/resources/mara/hg38.sites_no_conservation.bed.bz2 hg38.sites_no_conservation.bed.bz2
ln ~/projects/resources/mara/hg38.sites.bed.bz2 hg38.sites.bed.bz2

python3 ~/python/mara/make_profile.py --output=profile_enhancers.txt --input=enhancers.bed --tfbs=hg38.sites.bed.bz2 --upstream=500 --downstream=500 --symmetric
python3 ~/python/mara/make_profile.py --output=profile_enhancers_no_conservation.txt --input=enhancers.bed --tfbs=hg38.sites_no_conservation.bed.bz2 --upstream=500 --downstream=500 --symmetric

python3 ~/python/mara/associate_tfbs.py --output=tfbs_enhancers.txt --input=enhancers.bed --tfbs=hg38.sites.bed.bz2 --profile=profile_enhancers.txt
python3 ~/python/mara/associate_tfbs.py --output=tfbs_enhancers_no_conservation.txt --input=enhancers.bed --tfbs=hg38.sites_no_conservation.bed.bz2 --profile=profile_enhancers_no_conservation.txt
"


enhancer_tfbs = data.matrix(read.table("data-raw/mara/tfbs_enhancers.txt"))
mara_enhancer_tfbs = Matrix::drop0( enhancer_tfbs )
save(object = mara_enhancer_tfbs, file = "data-raw/mara/mara_enhancer_tfbs.rda", compress = "bzip2")
usethis::use_data( enhancer_tfbs   , internal = FALSE, overwrite = TRUE)

enhancer_tfbs_no_conservation = data.matrix(read.table("data-raw/mara/tfbs_enhancers_no_conservation.txt"))
mara_enhancer_tfbs_no_conservation = Matrix::drop0( enhancer_tfbs_no_conservation )
save(object = mara_enhancer_tfbs_no_conservation, file = "data-raw/mara/mara_enhancer_tfbs.rda", compress = "bzip2")
usethis::use_data( mara_enhancer_tfbs_no_conservation   , internal = FALSE, overwrite = TRUE)


