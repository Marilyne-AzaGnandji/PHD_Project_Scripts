#!/bin/bash

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads
#create a file named samples.txt
for R1 in *R1.fastq.gz ; do
echo -e "${R1/_*/}\t${R1}\t${R1/R1/R2}"
done > samples.txt
#headers add
sed -i $'1 i\\\nsample\tr1\tr2' samples.txt
cat samples.txt
#create config files for illumina-utils 
iu-gen-configs samples.txt -o 01_QC/
#run minoche quality filtering for each sample
for ini in 01_QC/*.ini; do
    iu-filter-quality-minoche $ini
done

exit 0
