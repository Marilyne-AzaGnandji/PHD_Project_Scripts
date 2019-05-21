#!/bin/bash -

# All informations about fastqc downloading and installation are available here: 
   ##https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# install fastqc and multiqc
#To install fastqc, i used this alternative
apt-get install fastqc
#All informations about multiqc downloading and installation are available here:
 ##https://github.com/ewels/MultiQC
#To install multiqc, i used this alternative:
git clone https://github.com/ewels/Multiqc.git
#then
cd MultiQC
#then to install its setup.py
python setup.py install



exit 0
