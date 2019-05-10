#!/bin/bash -
#Test if fastqc and multiqc are available on your computer
fastq #or fastqc -h or fastqc --help for more details
multiqc #or multiqc -h or multiqc --help for more details

# how to install?
#
# indicate here or in another script the command lines used to
# install fastqc and multiqc
#To install fastqc
apt-get install fastqc
#To install multiqc, i use this alternative:
git clone https://github.com/ewels/Multiqc.git
#then
cd MultiQC
#then to install its setup.py
python setup.py install

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
