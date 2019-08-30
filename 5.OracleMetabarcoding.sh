#!/bin/bash -

### We want to use Fred's metabarcoding pipeline to obtain OTU table from FASTQ files (with swarm)

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/Metabarcoding
## Check quality encoding (33 or 64?) : Our data are ARN 16S metadata

# Test if vsearch is available on your computer
VSEARCH=$(which vsearch) && \
    echo "vsearch is installed!" || \
	{ echo "Error: vsearch is not installed" ; exit 1 ; }


# Check quality encoding (33 or 64?)

for f in *fastq.gz ; do
"${VSEARCH}" \
      --fastq_chars ${f} 2> ${f/fastq.gz/log}
done


## Now we know how quality is encoded (offset of 33), we can merge paired-reads (R1 and R2)

# Merge read pairs
VSEARCH=$(which vsearch)
THREADS=4
ENCODING=33
for f in *R1_001.fastq.gz ; do
    FORWARD=$f
    REVERSE=${f/R1/R2}
    OUTPUT=${f/_*/_assembled.fastq}
    echo "${VSEARCH}" \
	 --threads ${THREADS} \
	 --fastq_mergepairs ${FORWARD} \
	 --reverse ${REVERSE} \
	 --fastq_ascii ${ENCODING} \
	 --fastqout ${OUTPUT} \
	 --fastq_allowmergestagger \
	 --quiet 2> ${OUTPUT/.fastq/.log}
done


## Demultiplexing , primer clipping, sample dereplication and quality extraction

exit 0


