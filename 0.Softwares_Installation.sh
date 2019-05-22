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
#Atropos installation
#Python3 at least and cython0.25.2+ are required to install atropos
sudo apt-get update
sudo apt-get install python3.6 #python3.6 installation
sudo apt-get cython #cython installation 
#In case there is an old version of python, you can make available python3.6:
$update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1
$update-alternatives --install /usr/bin/python python /usr/bin/python3.6 2 #when you put "2",you give priority to python version 3.6 that is currently needed
#To install atropos
sudo pip install atropos #sudo apt-get install pip if not installed
#To install illumina-utils
pip install illumina-utils



exit 0
