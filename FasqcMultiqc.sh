
#!/bin/bash
# how to install?
 ##You will see all informations about fastqc and multiqc installation in FastqcMultiqc_installation.sh
#Test if fastqc and multiqc are available on your computer
fastq #fastqc --version or fastqc -h or fastqc --help for more details
multiqc #multiqc --version or multiqc -h or multiqc --help for more details
<<<<<<< HEAD
=======

>>>>>>> 36f8b2bfa598d430566c645839a3b77c9baf971f
OUTPUT_FOLDER="FastqcResults"
# what if the folder already exists?
##we can check this folder and its content like this:
cat ${OUTPUT_FOLDER}
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
