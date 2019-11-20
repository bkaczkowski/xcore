# declare a name for this job to be sample_job
#PBS -N enh_jaspar
# request the queue (enter the possible names, if omitted, serial is the default)
#PBS -q serial
# request 1 node
#PBS -l nodes=1:ppn=20:lm4.q
# request 240 hours and 30 minutes of cpu time
#PBS -l cput=240:30:00

cd /home/bogumil/motif_enrichment/
~/R-3.6.0/bin/Rscript calculate_jaspar_rawscores_F5_promoters.R

# qsub -q bigmem.q  calculate_jaspar_rawscores_F5_promoters_qsub.sh

# scp data-raw/predicted_tfbs/jaspar/calculate_jaspar_rawscores_F5_promoters* dgt-ac4:~/motif_enrichment/
