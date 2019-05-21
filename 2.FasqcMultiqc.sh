#!/bin/bash -
#You will see all informations about fastqc and multiqc installation in FastqcMultiqc_installation.sh
#Test if fastqc and multiqc are available on your computer
which fastq && echo "OK" || "not installed"
which multiqc && echo "OK" || "not installed"
OUTPUT_FOLDER="FastqcResults"
#we can check this folder and its content like this:
sudo find / -type d -name "${OUTPUT_FOLDER}"
#In case this directory is found, do:
rm -r ${OUTPUT_FOLDER} #to delete
#In case this folder is not available;we can create easily one:
mkdir ${OUTPUT_FOLDER}
# quality control with fastqc
for f in *fastq.gz ; do
    fastqc --outdir ./${OUTPUT_FOLDER}/ "${f}"
done

# aggregate quality control results with multiqc
cd ./${OUTPUT_FOLDER}/
multiqc .

exit 0
