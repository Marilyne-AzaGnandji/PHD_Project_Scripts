#!/bin/bash -

# quality control with fastqc

# how to install?
#
# indicate here or in another script the command lines used to
# install fastqc and multiqc
#
OUTPUT_FOLDER="FastqcResults"
mkdir ${OUTPUT_FOLDER}  # what if the folder already exists?

for f in *fastq.gz ; do
    fastqc --outdir ./${OUTPUT_FOLDER}/ "${f}"
done

# aggregate quality control results with multiqc
cd ./${OUTPUT_FOLDER}/
multiqc .

exit 0
