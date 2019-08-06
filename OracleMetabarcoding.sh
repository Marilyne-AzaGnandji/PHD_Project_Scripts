#!/bin/bash -

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/test
## Check quality encoding (33 or 64?)

# Test if vsearch is available on your computer
VSEARCH=$(which vsearch) && \
    echo "vsearch is installed!" || \
	{ echo "Error: vsearch is not installed" ; exit 1 ; }


# Check quality encoding (33 or 64?)

for f in *fastq.gz ; do
"${VSEARCH}" \
      --fastq_chars ${f} 2> ${f/fastq.gz/log}
done


# Now we know how quality is encoded (offset of 33), we can merge paired-reads (R1 and R2).


#VSEARCH=$(which vsearch)
#THREADS=4
#ENCODING=33
#OUTPUT=$HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/test/*.fastq.gz
# Merge read pairs
#for f in *_R1_001.fastq.gz; do
#"${VSEARCH}" \
 #   --threads ${THREADS} \
  #  --fastq_mergepairs ${f} \
   # --reverse ${f/R1/R2} \
   # --fastq_ascii ${ENCODING} \
   # --fastqout ${OUTPUT} \
   # --fastq_allowmergestagger \
   # --quiet 2>> ${OUTPUT/.fastq.gz/.log}
#done

exit 0
