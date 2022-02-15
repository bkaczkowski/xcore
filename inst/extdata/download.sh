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
downloadDataset https://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_all_macs2_hg38_v1_0.bed.gz
# merged (non redundant) peaks hg38
downloadDataset https://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/remap2020_nr_macs2_hg38_v1_0.bed.gz
# EP300 hg38
downloadDataset https://remap.univ-amu.fr/storage/remap2020/hg38/MACS2/TF/EP300/remap2020_EP300_nr_macs2_hg38_v1_0.bed.gz

# Download ChIP-Atlas
# hg38 TFs and others (15217) in All cell types (61679) with 50 Threshold for Significance
downloadDataset http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/assembled/Oth.ALL.05.AllAg.AllCell.bed
# process ChIP-Atlas data to extract SRX field and remove unncecessary data from 4th column
if [ ! -f chip_atlas_hg38.Oth.ALL.05.AllAg.AllCell_SRX_only.bed.gz ]; then
  cat Oth.ALL.05.AllAg.AllCell.bed | \
    awk '{gsub(";.*|ID=","",$4)}1' |  \
    awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' | \
    tail -n +2 | \
    gzip > chip_atlas_hg38.Oth.ALL.05.AllAg.AllCell_SRX_only.bed.gz
fi
# download experiments metadata
downloadDataset http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab
if [ ! -f experimentList_TF_hg38.txt ]; then
  grep -w hg38 experimentList.tab | grep -w 'TFs and others' | awk '{if($2 == "hg38"){print $0}}' > experimentList_TF_hg38.txt
fi

# Download GENCODE annotation
downloadDataset http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz

# Download ENCODE Blacklist for hg38
downloadDataset https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz

# Protein Atlas data
downloadDataset https://www.proteinatlas.org/download/proteinatlas.tsv.zip

# CAGE peaks (promoters) expression tables
# The format is weird so I reprocessed it a bit
if [ ! -f hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz ]; then
  downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
  zcat hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz | grep -v -e "^##" -e "^01STAT" -e "^02STAT" | gzip > t && mv t hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz
fi

# CAGE peaks (promoters) counts tables
downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz

# enhancers expression tables
# The format is weird so I reprocessed it a bit
downloadDataset https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v8/extra/enhancer/F5.hg38.enhancers.expression.tpm.matrix.gz

# Download FANTOM CAT data
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv.gz

downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.CAGE_cluster.bed.gz
zcat FANTOM_CAT.lv3_robust.CAGE_cluster.bed.gz | awk '{printf "%s\t%s\t%s\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, $6}' | gzip > t && mv t FANTOM_CAT.lv3_robust.CAGE_cluster.bed.gz
if [ ! -f FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.bed.gz ]; then
  ./liftOver FANTOM_CAT.lv3_robust.CAGE_cluster.bed.gz hg19ToHg38.over.chain.gz FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.bed FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.unmapped.bed
  gzip FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.bed
  # remove record with outlayer width
  zcat FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.bed.gz | grep -v 'chr1:145176389..145176406,+' | gzip > t && mv t FANTOM_CAT.lv3_robust.CAGE_cluster.liftOver.hg38.bed.gz
fi

downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.only_divergent_p_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.only_e_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.only_intergenic_p_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.only_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/FANTOM_CAT.lv3_robust.only_mRNA.gtf.gz

downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/expression/expression_atlas/FANTOM_CAT.expression_atlas.CAGE_cluster.lv1_raw.rle_cpm.tsv.gz

# lv4_stringet
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.info_table.ID_mapping.tsv.gz

downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.CAGE_cluster.bed.gz
zcat FANTOM_CAT.lv4_stringent.CAGE_cluster.bed.gz | awk '{printf "%s\t%s\t%s\t%s\t%.0f\t%s\n", $1, $2, $3, $4, $5, $6}' | gzip > t && mv t FANTOM_CAT.lv4_stringent.CAGE_cluster.bed.gz
if [ ! -f FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.bed.gz ]; then
  ./liftOver FANTOM_CAT.lv4_stringent.CAGE_cluster.bed.gz hg19ToHg38.over.chain.gz FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.bed FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.unmapped.bed
  gzip FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.bed
  # remove record with outlayer width
  zcat FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.bed.gz | grep -v 'chr1:145176389..145176406,+' | gzip > t && mv t FANTOM_CAT.lv4_stringent.CAGE_cluster.liftOver.hg38.bed.gz
fi

downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.only_divergent_p_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.only_e_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.only_intergenic_p_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.only_lncRNA.gtf.gz
downloadDataset https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv4_stringent/FANTOM_CAT.lv4_stringent.only_mRNA.gtf.gz

# jaspar hg38 tfbs -- after download I have it intersected with F5 to lower the size down
downloadDataset http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/JASPAR2020_hg38.bb

# swissregulon tfbs
downloadDataset https://swissregulon.unibas.ch/data/hg38_f5/hg38_sites_v1.gff.gz

# GSE17708_Keshamouni_TGFB1_logs.xls
downloadDataset https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE17708&format=file&file=GSE17708%5FKeshamouni%5FTGFB1%5Flogs%2Exls%2Egz
if [ ! -f GSE17708_Keshamouni_TGFB1_logs.xls ]; then
	        gunzip GSE17708_Keshamouni_TGFB1_logs.xls.gz
fi
