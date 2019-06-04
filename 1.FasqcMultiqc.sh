#!/bin/bash -

cd /home/u082-f048/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads

# You will see all informations about fastqc and multiqc installation in Softwares_installation.sh

# Test if fastqc and multiqc are available on your computer
which fastq && echo "fastqc is installed!" || exit 1
which multiqc && echo "multiqc is installed" || exit 1
OUTPUT_FOLDER="FastqcResults"

# we can check this folder and its content like this:
sudo find / -type d -name "${OUTPUT_FOLDER}"

# In case this directory is found, do:
ls ${OUTPUT_FOLDER} && rm -r ${OUTPUT_FOLDER} #check the content and delete the old directory;we need an empty available directory

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
