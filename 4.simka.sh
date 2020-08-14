#!/bin/bash -

## We want to create in a repertory named Data 3 folders: SAMA_12k_reads_subsampling_seed_{1..3}

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads
# Create the folder Data and move in it
mkdir Data

# Copy trimmed .fastq.gz files ("trimmed" extension is due to performing sequences quality filtering with atropos & illumina-utils)

cp *_trimmed_R1.fastq.gz *_trimmed_R2.fastq.gz ./Data
cd ./Data && gunzip *
cd ..

cp -R Data ./simka
rm -r Data


cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka/Data/
#Performing subsampling at 1% with vsearch
mkdir SAMA_12_one_percent_subsampling_seed_1 SAMA_12_one_percent_subsampling_seed_2 SAMA_12_one_percent_subsampling_seed_3
cp *.fastq ./SAMA_12_one_percent_subsampling_seed_1
cp *.fastq ./SAMA_12_one_percent_subsampling_seed_2
cp *.fastq ./SAMA_12_one_percent_subsampling_seed_3

FOLDER="SAMA_12_one_percent_subsampling_seed"
# Here i need to review the subsampling test taking into account the zipped files (just need to add --gzip_decompress option)
for i in 1 2 3 ; do
    for f in ./${FOLDER}_${i}/*.fastq ; do
         vsearch --gzip_decompress\
           --fastx_subsample "$f" \
           --fastqout "$f"_subsampling \
           --sample_pct 0.01
    done
done

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka
cp -R ./Data/SAMA_12_one_percent_subsampling_seed_* . && rm ./Data/SAMA_12_one_percent_subsampling_seed_*/*.fastq
rm -r Data

cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka

# Rename files extensions: *subsampling in .fastq
FOLDER="SAMA_12_one_percent_subsampling_seed"
for i in 1 2 3 ; do
    for f in ./${FOLDER}_${i}/*.fastq_subsampling ; do
	newfiles=${f/.fastq_subsampling/.fastq}
	mv $f $newfiles
    done
done




#cd $HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka
## simka folder
SIMKA_FOLDER="$HOME/Bureau/marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka"

## create a list of sets for simka
FOLDER="SAMA_12_one_percent_subsampling_seed"
for i in 1 2 3 ; do
    for f in ./${FOLDER}_${i}/*_R1.fastq ; do
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

#cd $HOME/Documents/Marilyne/PhD_Thesis/SAMA_12_first_10k_reads/simka
#set -x
# simka results visualization
for i in 1 2 3; do
 python2 scripts/visualization/run-visualization.py -in simkamin_${i} -out simkamin${i} -pca -heatmap -tree 
done

exit 0



