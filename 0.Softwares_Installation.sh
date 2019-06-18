#!/bin/bash -

## install fastqc
# (see https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
sudo apt-get update
sudo apt-get install fastqc

## install multiqc
# (see https://github.com/ewels/MultiQC)
# requires python 2+, pip (for setuptools packages)
sudo apt-get update
sudo apt-get install python2 python-pip
git clone https://github.com/ewels/Multiqc.git
cd ./Multiqc/
sudo python2 setup.py install

## install atropos
# (see https://github.com/jdidion/atropos)
# requires python3.3+ and cython0.25.2+
sudo apt-get update
sudo apt-get install python3 cython python3-pip # currently python 3.7 and cython 0.28
sudo pip3 install atropos
# give priority to python 3 over python 2:
sudo update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.7 2 # 2 has priority over 1

## install illumina-utils
sudo pip3 install illumina-utils

## install simka/simkaMin
# (see https://github.com/GATB/simka/tree/master/simkaMin)
# requires cmake2.6+ and gcc4.4.7+
sudo apt-get update
sudo apt-get -y install cmake
sudo apt-get -y install gcc
git clone https://github.com/GATB/simka.git
cd simka
sh INSTALL
cd example
./simple_test.sh # To test the sofware on your computer
cd ..
sudo python simkaMin/simkaMin.py # To see simka in-line help


exit 0
