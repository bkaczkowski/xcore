# OBSOLETE script

cat remap2018_all_macs2_hg38_v1_2.bed | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > remap2018_all_macs2_hg38_v1_2.bed6

cut -f 4 remap2018_all_macs2_hg38_v1_2.bed6  | sort | uniq >  Experiment_list.txt &
cut -f 4 remap2018_all_macs2_hg38_v1_2.bed6  | awk '{gsub(".*\\.","",$1)}1'|  sort | uniq >  Cell_list.txt &


mkdir -p split_by_cell
for cell in `cat Cell_list.txt`; do
  grep  $cell remap2018_all_macs2_hg38_v1_2.bed6  > split_by_cell/$cell.bed
done


mkdir -p split_by_experiment
for cell_bed in `ls split_by_cell/*.bed`; do
  for exp in `cut -f 4 cell_bed | sort | uniq`; do
    grep -w $exp $cell_bed > split_by_experiment/$exp.bed
  done
done
