#!/bin/bash -

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads

# Test if iu-filter-quality-minoche is available on your computer
which iu-filter-quality-minoche && \
    echo "iu-filter_quality-minoche is installed!" || \
   { echo "Error: iu-filter-quality-minoche is missing because illumina-utils is not installed"; exit 1 ; }
# To check the version
iu-filter-quality-minoche -v # the answer in my case is Illumina-utils v2.6

# We need first to create a file named samples.txt
(printf "sample\tr1\tr2\n"
for R1 in *R1.fastq; do
  echo -e "${R1/_*/}\t${R1}\t${R1/R1/R2}"
done)> samples.txt

# cat samples.txt

# create config files for illumina-utils 
iu-gen-configs samples.txt -o 01_QC/

# run minoche quality filtering for each sample
for ini in 01_QC/*.ini; do
  iu-filter-quality-minoche $ini
done

exit 0
