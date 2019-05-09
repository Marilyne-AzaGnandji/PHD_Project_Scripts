#!/bin/bash -

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads

  #Perform quality trimming
  for f in *_R1.fastq.gz;
  do
   atropos trim -q 20,20 -o output1.fastq -p output2.fastq -pe1 "${f}" -pe2 "${f/R1/R2}"
  done

exit 0
