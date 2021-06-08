# Bogumil Kaczkowski August 22, 2019
# bash commands to
# 1) map fastq files to human genome
# 2) uplaad the data to ZENBU for visualization
# 3) do peak calling with MACS2


# MAPPING THE SEQUENCES (FATSTQ) TO THE GENOME
/home/scratch/moirai/bin/symbolic_link_files.pl /home/scratch/moirai/genome/hg38_male.f*
for file in `ls  *.fastq.gz`; do
/home/scratch/moirai/bin/bwa_0.7.4/bwa aln -t 4  hg38_male.fa $file > $file.sai
/home/scratch/moirai/bin/bwa_0.7.4/bwa samse -n 10000 hg38_male.fa $file.sai $file > $file.sam
#samtools view -bS  $file.sam | samtools sort - -o $file.bam
samtools view -b -q 20 $file.sam  | samtools sort - -o $file.bam # wit filtering for reads with low mapping quality
samtools index $file.bam
rm $file.sam
rm $file.sai
done

# ZENBU UPLOAD
ls *.bam > bam_files.txt
zenbu_upload -url "http://fantom.gsc.riken.jp/zenbu"  -assembly hg38  -platform Chipseq  -filelist bam_files.txt


# PEAK CALLING (default)
for file in `ls   *q20.bam`; do
    macs2 callpeak -t $file -f BAM -g hs -p 0.001  -n $file
done

# PEAK CALLING (no-model, fixed extension size)
for file in `ls   *.bam`; do
    macs2 callpeak -t $file -f BAM -g hs -p 0.001  -n $file.ext300 --nomodel --extsize 300
done

# ALTERNATIVE, if q20 filter was not apply to bam files, this for loop can do it
# fitering for single mapping reads (Q20 filter)
# should be performed prior to peak calling
#for file in `ls   *.bam`; do
#    # removing emply spaces from the headers
#    file="$(echo $file | sed 's/.bam//g')"
#    echo $file
#    samtools view -b -q 20 $file.bam  > $file.q20.bam
# done
