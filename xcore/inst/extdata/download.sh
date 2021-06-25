#!/usr/bin/bash
# Download external dataset including ReMap peaks, 

## Check if the dataset, given as the URL, has been already downloaded and download it if it has not.
function downloadDataset {
  dataset=`basename $1`
  if [ ! -f ${dataset} ]; then
    wget $1
  fi
}

# Download FANTOM5 annotations
# enhancers bed hg38
downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/enhancer/F5.hg38.enhancers.bed.gz
# promoters bed hg38
downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz
# downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks/hg38_liftover+new_CAGE_peaks_phase1and2.bed.gz
# promoters annotation hg38
downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz
# downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks_annotation/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz

# Download ReMap
# all peaks hg38
downloadDataset http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_all_macs2_hg38_v1_0.bed.gz
# merged (non redundant) peaks hg38
downloadDataset http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_nr_macs2_hg38_v1_0.bed.gz
# EP300 hg38
downloadDataset https://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/TF/EP300/remap2020_EP300_nr_macs2_hg38_v1_0.bed.gz

# Download ChIP-Atlas
# hg38 TFs and others (15217) in All cell types (61679) with 50 Threshold for Significance
downloadDataset http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/assembled/Oth.ALL.05.AllAg.AllCell.bed

# Download GENCODE annotation
downloadDataset http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz
