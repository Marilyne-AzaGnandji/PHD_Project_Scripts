#!/bin/bash -

cd ${HOME}/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/

# Test if atropos is available on your computer
which atropos && \
    echo "atropos is installed!" || \
        { echo "Error: atropos is not installed" ; exit 1 ; }

# Perform quality trimming (use full length option names)
<<<<<<< HEAD
for f in *_R1.fastq.gz ; do  
=======
for f in *_R1.fastq.gz ; do    
>>>>>>> b45f2b93119a169823c4f96d40844c095ff767da
## objective: file_R1.fastq.gz -> file_trimmed_R1.fastq.gz   
  atropos \
      trim \
      --quality-cutoff 20,20 \
      --output ${f/R1/trimmed_R1} \
      --paired-output ${f/R1/trimmed_R2} \
      --input1 "${f}" \
      --input2 "${f/R1/R2}"
done

exit 0

