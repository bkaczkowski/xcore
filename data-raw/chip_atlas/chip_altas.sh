
# 1)  Oth.ALL.05.AllAg.AllCell.bed was obtained by downloading all Transcription Factor peaks [TF and others(9868)]
# from http://chip-atlas.org/peak_browser the file is ~ 72G uncompressed and takes some hours to download
# the peaks coordinates are of the hg19 genome https://github.com/inutano/chip-atlas/wiki

# 2) "grep"ing the 4th collumn on the huge file to extract the TF and cell type/sample information,
# keeping the 5th SCORE column for future quality filtering
# grep in awk in split into seperate chunk to quicker and clearere processing

cat Oth.ALL.05.AllAg.AllCell.bed | awk '{gsub("\\);.*","",$4)}1' |  awk '{gsub(";Name=","___Tf_",$4)}1' |  awk '{gsub("\\%20\\(\\@\\%20","__",$4)}1' |  awk '{gsub("ID=","",$4)}1' | awk '{gsub("%20","-",$4)}1'| awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > Oth.ALL.05.AllAg.AllCell_extracted.bed

# remove first line (header)
tail -n +2 Oth.ALL.05.AllAg.AllCell_extracted.bed > Oth.ALL.05.AllAg.AllCell_extracted_hg19.bed
rm Oth.ALL.05.AllAg.AllCell_extracted.bed

# liftOver to hg38
./liftOver  Oth.ALL.05.AllAg.AllCell_extracted_hg19.bed hg19ToHg38.over.chain Oth.ALL.05.AllAg.AllCell_extracted_hg19ToHg38_liftOver.bed unMapped

# create SRX_to_GSM_map.txt
cut -f 4 Oth.ALL.05.AllAg.AllCell.bed |awk '{gsub(":%20.*","",$1)}1' | awk '{gsub(";Name.*Title=","_",$1)}1' |  awk '{gsub("ID=","",$1)}1' | grep _GSM | sort | uniq > SRX_to_GSM_map.txt &

# create Transcription_Factors_list.txt
cut -f 4 Oth.ALL.05.AllAg.AllCell_extracted_hg19.bed |awk '{gsub(".*___","",$1)}1' | awk '{gsub("__.*","_",$1)}1' | sort | uniq >  Transcription_Factors_list.txt &

mkdir -p split_by_ft
for tf in `cat Transcription_Factors_list.txt`; do
  grep  $tf Oth.ALL.05.AllAg.AllCell_extracted_hg19ToHg38_liftOver.bed  > split_by_ft/$tf.bed
done

mkdir -p split_by_experiment
for tf_bed in `ls split_by_ft/*.bed`; do
  for exp in `cut -f 4 $tf_bed | sort | uniq`; do
    grep $exp $tf_bed  > split_by_experiment/$exp.bed
  done
done

gzip split_by_experiment/*bed

# The above was done on the server so we dowload the data to local computer, the copy is on the server
scp -r dgt-ac4:/work/bogumil/chip_atlas/split_by_experiment ~/projects/resources/chip_atlas/
