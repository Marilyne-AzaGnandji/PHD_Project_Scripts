#!/bin/bash -

# We want to create 3 directories containing our data (paired and reads)
cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads 

# SAMA_12k_reads_subsampling_seed_1: first directory
mkdir SAMA_12k_reads_subsampling_seed_1
cd SAMA_12k_reads_subsampling_seed_1 && cp ../*_R1.fastq.gz ../*_R2.fastq.gz .

# SAMA_12k_reads_subsampling_seed_2: second directory
cd ..
mkdir SAMA_12k_reads_subsampling_seed_2
cd SAMA_12k_reads_subsampling_seed_2 && cp ../*_R1.fastq.gz ../*_R2.fastq.gz .

# SAMA_12k_reads_subsampling_seed_3: third directory
cd ..
mkdir SAMA_12k_reads_subsampling_seed_3
cd  SAMA_12k_reads_subsampling_seed_3 && cp ../*_R1.fastq.gz ../*_R2.fastq.gz .

# The directory Data contains the 3 directories created
cd ..
mkdir Data && cd Data
mv ../SAMA_12k_reads_subsampling_seed_1 .
mv ../SAMA_12k_reads_subsampling_seed_2 .
mv ../SAMA_12k_reads_subsampling_seed_3 .

cd
cd $HOME/Bureau/marilyne/PhD_Thesis/simka
cp -R ../SAMA_12_first_10k_reads/Data . && rm -r ../SAMA_12_first_10k_reads/Data

cd $HOME/Bureau/marilyne/PhD_Thesis/simka/Data

## simka folder
SIMKA_FOLDER="$HOME/Bureau/marilyne/PhD_Thesis/simka"

## create a list of sets for simka
FOLDER="SAMA_12k_reads_subsampling_seed"
for i in 1 2 3 ; do
    for f in ./${FOLDER}_${i}/*_R1.fastq.gz ; do
        SAMPLE="${f/.*BU/BU}"
        SAMPLE="${SAMPLE/_*/}"
        echo ${SAMPLE}" : "$(readlink -f $f)" ; "$(readlink -f ${f/_R1/_R2})
    done > sets_${i}.config
done

## create output folders
rm -rf simkamin_{1..3}
mkdir simkamin_{1..3}

# run the pairwise sample comparisons
for i in 1 2 3 ; do
    python2 \
        ${SIMKA_FOLDER}/simkaMin/simkaMin.py \
        -bin ${SIMKA_FOLDER}/build/bin/simkaMin \
        -in sets_${i}.config \
        -out ./simkamin_${i}/
done 

exit 0
