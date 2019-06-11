#!/bin/bash -

cd ${HOME}/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/

# Test if atropos is available on your computer
which atropos && \
    echo "atropos is installed!" || \
    { echo "Error: atropos is not installed" ; exit 1 ; }

# Perform quality trimming (use full length option names)
for f in *_R1.fastq.gz ; do
    #OUTPUT_FILE1="${f:0:-9}"
    OUTPUT_FILE1="${f/_*/}"
    OUTPUT_FILE2="$OUTPUT_FILE1"
    
## objective: file_R1.fastq.gz -> file_trimmed_R1.fastq.gz   
  atropos \
      trim \
      --quality-cutoff 20,20 \
      --output ${OUTPUT_FILE1}_trimmed_R1.fastq.gz \
      --paired-output ${OUTPUT_FILE2}_trimmed_R2.fastq.gz \
      --input1 "${f}" \
      --input2 "${f/R1/R2}"
done

exit 0

