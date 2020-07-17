#!/bin/bash -

### We want to use Fred's metabarcoding pipeline to obtain OTU table from FASTQ files (with swarm)

#cd $HOME/Documents/Marilyne/PhD_Thesis/SAMA_12_first_10k_reads/Metabarcoding
## Check quality encoding (33 or 64?) : Our data are ARN 16S metadata
cd $HOME/Documents/marilyne/Metabarcoding/
# Test if vsearch is available on your computer
VSEARCH=$(which vsearch) && \
    echo "vsearch is installed!" || \
	{ echo "Error: vsearch is not installed" ; exit 1 ; }


# Check quality encoding (33 or 64?)

for f in *fastq.gz ; do
    "${VSEARCH}" \
	--fastq_chars ${f} 2> ${f/fastq.gz/log}
done


## Now we know how quality is encoded (offset of 33), we can merge paired-reads (R1 and R2)

# Merge read pairs
VSEARCH=$(which vsearch)
THREADS=4
ENCODING=33
for f in *R1_001.fastq.gz ; do
    FORWARD=$f
    REVERSE=${f/R1/R2}
    OUTPUT=${f/_*/_assembled.fastq}
     "${VSEARCH}" \
	 --threads ${THREADS} \
	 --fastq_mergepairs ${FORWARD} \
	 --reverse ${REVERSE} \
	 --fastq_ascii ${ENCODING} \
	 --fastqout ${OUTPUT} \
	 --fastq_allowmergestagger \
	 --quiet 2> ${OUTPUT/.fastq/.log}
done 

#cd $HOME/Documents/Marilyne/PhD_Thesis/SAMA_12_first_10k_reads/Metabarcoding
cd $HOME/Documents/marilyne/Metabarcoding/
##  primer clipping, sample dereplication and quality extraction
#set -x

# Define binaries, temporary files and output files
MIN_LENGTH=32
CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH}"
VSEARCH=$(which vsearch)
INPUT_REVCOMP=$(mktemp)
TMP_FASTQ=$(mktemp)
TMP_FASTQ1=$(mktemp)
TMP_FASTQ2=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT=$(mktemp)
PRIMER_F="CCTACGGGNGGCWGCAG"
ANTI_PRIMER_R="GGATTAGATACCCBDGTAGTC" 
#PRIMER_R="GACTACHVGGGTATCTAATCC"

for INPUT in *_assembled.fastq ; do    
    MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))  # primer match is >= 2/3 of primer length
    MIN_R=$(( ${#ANTI_PRIMER_R} * 2 / 3 ))
    QUALITY_FILE="${INPUT/.fastq/.qual}"
    FINAL_FASTA="${INPUT/_assembled.fastq/.fas}"
    LOG="${INPUT/_*/.log}"
    LOG1="${INPUT/_*/_final.log}"

    # # Reverse complement fastq file
    "${VSEARCH}" --quiet \
    		 --fastx_revcomp "${INPUT}" \
    		 --fastqout "${INPUT_REVCOMP}"

    # Trim forward & reverse primers (search normal and antisens)
    cat "$INPUT" "${INPUT_REVCOMP}" | \
        ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" -o "${TMP_FASTQ1}" - > "${LOG}"
    
    
    cat "${TMP_FASTQ1}" | \
    	${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" -o "${TMP_FASTQ}" - >> "${LOG}"
    #cat "${LOG}"

    
    # Discard sequences containing Ns, add expected error rates
    "${VSEARCH}" \
    	--quiet \
    	--fastq_filter "${TMP_FASTQ}" \
    	--fastq_maxns 0 \
    	--relabel_sha1 \
    	--eeout \
        --log /dev/stdout \
    	--fastqout "${TMP_FASTQ2}" > "${LOG1}"
    #cat "${LOG1}" 
    # Discard sequences containing Ns, convert to fasta
    "${VSEARCH}" \
    	--quiet \
    	--fastq_filter "${TMP_FASTQ}" \
    	--fastq_maxns 0 \
    	--fastaout "${TMP_FASTA}" >> "${LOG1}"
    
    # Dereplicate at the study level
    "${VSEARCH}" \
        --quiet \
        --derep_fulllength "${TMP_FASTA}" \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${FINAL_FASTA}" >> "${LOG1}"
    # cat "${LOG1}"


    #Discard quality lines, extract hash, expected error rates and read length
    sed 'n;n;N;d' "${TMP_FASTQ2}" | \
        awk 'BEGIN {FS = "[;=]"}
             {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
    		 tr -d "@" > "${OUTPUT}"
    cat "${OUTPUT}"

    # Produce the final quality file
    sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
	uniq --check-chars=40 > "${QUALITY_FILE}"

done

# Clean
rm -f "${INPUT_REVCOMP}" "${TMP_FASTQ}" "${TMP_FASTA}" "${TMP_FASTQ2}" "${OUTPUT}"

#cd $HOME/Documents/Marilyne/PhD_Thesis/SAMA_12_first_10k_reads/Metabarcoding
## Global dereplication, clustering and chimera detection
cd $HOME/Documents/marilyne/Metabarcoding/
VSEARCH=$(which vsearch)
SWARM=$(which swarm)
TMP_FASTA=$(mktemp --tmpdir=".")
FINAL_FASTA="combined_samples.fas"

# Pool sequences
#cat *.fas > "${TMP_FASTA}" 

# Dereplicate (vsearch)
cat *.fas | \
    "${VSEARCH}" \
        --derep_fulllength - \
        --sizein \
	--sizeout \
	--fasta_width 0 \
	--output "${FINAL_FASTA}" > /dev/null

# "${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
    #              --sizein \
    #              --sizeout \
    #              --fasta_width 0 \
    #              --output "${FINAL_FASTA}" > /dev/null

#rm -f "${TMP_FASTA}"

#cd $HOME/Documents/Marilyne/PhD_Thesis/SAMA_12_first_10k_reads/Metabarcoding
cd $HOME/Documents/marilyne/Metabarcoding/
# Clustering
SWARM=$(which swarm) #
VSEARCH=$(which vsearch) #
THREADS=16
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
FINAL_FASTA="combined_samples.fas"
"${SWARM}" \
    -d 1 --fastidious -y 16 -b 3 -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s ${FINAL_FASTA/.fas/_1f.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}


# Sort representatives
"${VSEARCH}" --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/_1f_representatives.fas}
rm ${TMP_REPRESENTATIVES}
# Chimera checking
REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
"${VSEARCH}" --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"

# Define variables, temporary files and output files
cd $HOME/Documents/marilyne/Metabarcoding/
RELEASE=132
URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${RELEASE}/Exports"
INPUT="SILVA_${RELEASE}_SSURef_Nr99_tax_silva.fasta.gz"
PRIMER_F_NAME="341F_785R"
OUTPUT="${INPUT/.fasta.gz/}_${PRIMER_F_NAME}.fas"
LOG="${OUTPUT/.fas/.log}"
PRIMER_F="CCTACGGGNGGCWGCAG"
PRIMER_R="GGATTAGATACCCBDGTAGTC"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH}"

# Download
[[ -f "${INPUT}" ]] || wget \
     --no-clobber \
     --continue \
     "${URL}/${INPUT}" "${URL}/${INPUT}.md5"

# Check
md5sum -c "${INPUT}.md5"

# Trim forward & reverse primers
zcat "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
     ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
     ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_F}" - 2>> "${LOG}" | \
     sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${OUTPUT}"

# Stats
grep -c "^>" "${OUTPUT}"
set -x
cd $HOME/Documents/marilyne/Metabarcoding/
FINAL_FASTA="combined_samples_1f_representatives.fas"
DATABASE="SILVA_132_SSURef_Nr99_tax_silva_341F_785R.fas"
FILTER=2
THREADS=2
VSEARCH=$(which vsearch)
set -x "${VSEARCH}" \
     --fastx_filter "${FINAL_FASTA}" \
     --minsize "${FILTER}" \
     --quiet \
     --fastaout - | \
 set -x  "${VSEARCH}" \
     --usearch_global - \
     --threads ${THREADS} \
     --dbmask none \
     --qmask none \
     --rowlen 0 \
     --notrunclabels \
     --userfields query+id1+target \
     --maxaccepts 0 \
     --maxrejects 32 \
     --top_hits_only \
     --output_no_hits \
     --db ${DATABASE} \
     --id 0.5 \
     --iddef 1 \
     --userout - > hits.representatives
set -x
 # in case of multi-best hit, find the last-common ancestor
cd $HOME/Documents/marilyne/Metabarcoding/
python3 ./stampa_merge.py $(pwd)

# sort by decreasing abundance
cd $HOME/Documents/marilyne/Metabarcoding/
RESULTS="tailhitsrepresentatives.txt"
sort -k2,2nr -k1,1d results.representatives > "${RESULTS}"

# clean
rm {hits,results}.representatives

# Build the OTU table
set -x
cd $HOME/Documents/marilyne/Metabarcoding/
FASTA="combined_samples.fas"
SCRIPT="OTU_contingency_table.py"
STATS="combined_samples_1f.stats"
SWARMS="combined_samples_1f.swarms"
REPRESENTATIVES="combined_samples_1f_representatives.fas"
UCHIME="combined_samples_1f_representatives.uchime"
ASSIGNMENTS="results.representatives"
QUALITY="combined_samples.qual"
OTU_TABLE="${FASTA/.fas/.OTU.table}"

 python \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${SWARMS}" \
    "${UCHIME}" \
    "${QUALITY}" \
    "${ASSIGNMENTS}" \
  OR-*.fas [T]-*.fas > "${OTU_TABLE}"

 # Filter the OTU table
cd $HOME/Documents/marilyne/Metabarcoding/
TABLE="combined_samples.OTU.table"
FILTERED="${TABLE/.table/.filtered.table}"

head -n 1 "${TABLE}" > "${FILTERED}" 
cat "${TABLE}"| awk '$7 == "N" && $9 <= 0.0002 && ($2 >= 3 || $8 >= 2)' >> "${FILTERED}"

exit 0


