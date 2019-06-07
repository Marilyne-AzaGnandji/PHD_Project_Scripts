#!/bin/bash -

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads

# Test if atropos is available on your computer
which atropos && echo "atropos is installed!" || if echo "!!! atropos is not installed"; then exit 1; fi

# Detect potential adapters in my paired-end-reads
for f in *_R1.fastq.gz;
do
  atropos detect -pe1 "${f}" -pe2 "${f/R1/R2}"
done

# Perform quality trimming
for f in *_R1.fastq.gz;
do
  OUTPUT_FILE1="${f:0:-9}"
  OUTPUT_FILE2="${OUTPUT_FILE1/R1/R2}"
  atropos trim -q 20,20 -o ${OUTPUT_FILE1}.fastq -p ${OUTPUT_FILE2}.fastq -pe1 "${f}" -pe2 "${f/R1/R2}"
done

exit 0

