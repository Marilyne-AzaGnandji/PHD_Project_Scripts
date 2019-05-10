#!/bin/bash -

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads
#Test to detect potential adapters in my paired-en-reads
for f in *_R1.fastq.gz;
do
atropos detect -pe1 "${f}" -pe2 "${f/R1/R2}"
done>AdapterDetection.txt
#Perform quality trimming
#SampleNames=
 # OUTPUT_FILE1=*_R1.fastq.gz
#  for f in *_R1.fastq.gz;
 # do
  # atropos trim -q 20,20 -o ${OUTPUT_FILE1}.fastq -p ${OUTPUT_FILE2}.fastq -pe1 "${f}" -pe2 "${f/R1/R2}"
 # done

exit 0

