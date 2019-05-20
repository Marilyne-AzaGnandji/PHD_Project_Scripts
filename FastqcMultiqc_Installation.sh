#!/bin/bash

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



exit 0
