#!/bin/bash -

# We want to create in a repertory named Data 3 folders: SAMA_12k_reads_subsampling_seed_{1..3}

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads
# Create repertory Data and move in it
mkdir Data && cd ./Data
# Copy trimmed .fastq.gz files ("trimmed" extension is due to performing sequences quality with atropos & illumina-utils)
cp ../*trimmed_* . | gunzip *
mv Data ../../simka

cd $HOME/Bureau/marilyne/PhD_Thesis/simka/Data
 #Performing subsampling at 1% with vsearch
 
for f in *fastq; do
       vsearch \
           --fastx_subsample $f \
           --fastqout "$f"_subsampling \
           --sample_pct 0.01 /
done
 
mkdir SAMA_12_one_percent_subsampling_seed_1 SAMA_12_one_percent_subsampling_seed_2 SAMA_12_one_percent_subsampling_seed_3
cp -r *_subsampling ./SAMA_12_one_percent_subsampling_seed_1 ./SAMA_12_one_percent_subsampling_seed_2 ./SAMA_12_one_percent_subsampling_seed_3
rm *_subsampling

cd $HOME/Bureau/marilyne/PhD_Thesis/simka

## simka folder
SIMKA_FOLDER="$HOME/Bureau/marilyne/PhD_Thesis/simka"

## create a list of sets for simka
FOLDER="SAMA_12_one_percent_subsampling_seed"
for i in 1 2 3 ; do
    for f in ./${FOLDER}_${i}/*_R1.fastq_subsampling ; do
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
