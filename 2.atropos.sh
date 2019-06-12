#!/bin/bash -

cd ${HOME}/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/

# Test if atropos is available on your computer
which atropos && \
    echo "atropos is installed!" || \
        { echo "Error: atropos is not installed" ; exit 1 ; }

# Perform quality trimming (use full length option names)
for f in *_R1.fastq.gz ; do
    OUTPUT_FILE1=${f/R1/trimmed_R1}
    OUTPUT_FILE2=${f/R1/trimmed_R2}
   #echo "$OUTPUT_FILE1"
    
## objective: file_R1.fastq.gz -> file_trimmed_R1.fastq.gz   
  atropos \
      trim \
      --quality-cutoff 20,20 \
      --output ${OUTPUT_FILE1} \
      --paired-output ${OUTPUT_FILE2} \
      --input1 "${f}" \
      --input2 "${f/R1/R2}"
done

exit 0

