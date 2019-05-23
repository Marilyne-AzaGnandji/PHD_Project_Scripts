#!/bin/bash -
# All informations about fastqc downloading and installation are available here: 
   ##https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# install fastqc and multiqc
#To install fastqc, i used this alternative
sudo apt-get install fastqc
#All informations about multiqc downloading and installation are available here:
 ##https://github.com/ewels/MultiQC
#To install multiqc, i used this alternative:
git clone https://github.com/ewels/Multiqc.git
#python(version 2 at least),pip(to get for instance setuptools packages) are required
sudo apt-get install python
sudo apt install python-pip
#then
cd Multiqc/
#then to install its setup.py
sudo python setup.py install
#Atropos installation
#Python3.3+ at least and cython0.25.2+ are required to install atropos
sudo apt-get update
sudo apt-get install python3 #python3.7 installation(my current update version)
sudo apt-get install cython #cython installation 
#In case there is an old version of python, you can make available python3.7:
update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1
update-alternatives --install /usr/bin/python python /usr/bin/python3.7 2 #when you put "2",you give priority to python version 3.7 that is currently needed
#To install atropos
sudo apt install python3-pip #sudo apt install python3-pip if not installed
sudo pip install atropos 
#To install illumina-utils
sudo pip install illumina-utils



exit 0
