# extract SRX and remove the unnecessary data from 4th column
zcat Oth.ALL.05.AllAg.AllCell.bed.gz  | awk '{gsub(";.*|ID=","",$4)}1' |  awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6}' > Oth.ALL.05.AllAg.AllCell_SRX_only.bed
# remove the first line (dummy header)
tail -n +2  Oth.ALL.05.AllAg.AllCell_SRX_only.bed > Oth.ALL.05.AllAg.AllCell_SRX_only_hg19.bed
#liftOver to hg38
 ~/bin/liftOver  Oth.ALL.05.AllAg.AllCell_SRX_only_hg19.bed  ~/projects/resources/liftOver_chains/hg19ToHg38.over.chain Oth.ALL.05.AllAg.AllCell_SRX_only_hg19ToHg38_liftOver.bed unMapped
