#!/usr/bin/bash

# For reliable, robust and maintainable bash scripts, start with following commands
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
# set -euo pipefail
# IFS=$"\n\t"

#########################################################################################
# HEADER
#########################################################################################
#% 
#% DESCRIPTION
#% This bash file that stiches together steps to analyze Cut&Run experimental data.
#% -------------------
#% USAGE and EXAMPLE
#% bash ${SCRIPT_NAME} [-h] /path/to/input/folder/ File name handle? [number] Run QC? [y | n] Run preprocessing? [y | n] Run BAM processing? [y | n] Run Peakcalling? [y | n]
#% ./${SCRIPT_NAME} [-h] /path/to/input/folder/ File name handle? [number] Run QC? [y | n] Run preprocessing? [y | n] Run BAM processing? [y | n] Run Peakcalling? [y | n]
#% -------------------
#% OPTIONS
#% 1) /path/to/input/folder/: Directory containing all raw sequence FASTQ files
#% 2) File name handle? [number]: For filenames starting with "TF_TAG_REPLICATE-NUM" like
#% "TEC1_IGG_1", then enter 3. For filenames starting with
#% "TF_TAG_CONDITION_REPLICATE-NUM" like "DPB4_MC_IGG_2", then enter 4. This is useful to
#% associate each GFP with associated IGG file. Files will be identified from beginning
#% with TF name up until REPLICATE-NUM value.
#% 3) Run QC? [y | n]: Should FASTQ QC metric be run? Please write in y for Yes or n for No
#% 4) Run preprocessing? [y | n]: Do you want to run adapter trimming and alignment steps
#% on raw FASTQ files? Please write in y for Yes or n for No
#% 5) Run BAM processing? [y | n]: Do you want to run spike-in calibration and PCR
#% deduplication of alignment BAM files? Please write in y for Yes or n for No
#% 6) Run Peakcalling? [y | n]: Do you want to run peakcalling on all calibrated BAM
#% files? Please write in y for Yes or n for No
#% 7) Smart compare? [y | n] Do you want to compare GFP file processing to IGG file with
#% least E. coli spike-in sequences, if its IGG sample is not present? Please write in y
#% for Yes or n for No. If not skipping, other IGG replicate file will be used, whichever
#% IGG contains lower spike-in read count.
#% 
#+ ---------------------------------------------------------------------------------------
#+ SCRIPT INFORMATION
#+ 
#+ VERSION: 0.5.9
#+ AUTHOR:  Akshay Paropkari
#+ LICENSE: GPL v3
#+ 
#+ ---------------------------------------------------------------------------------------
#+ REQUIRED TOOLS - which can be installed via Conda package management software
#+ 
#+ 1. FastQC - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#+ 2. MultiQC - https://multiqc.info/docs/
#+ 3. Trimmomatic - https://github.com/usadellab/Trimmomatic
#+ 4. Picard - https://github.com/broadinstitute/picard
#+ 5. Macs2 tutorial - https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2021/ChIPSeq/practicals/ChIP_Practical1_peakcalling_2021.html#check-the-input-data
#+ 6. deepTools - https://deeptools.readthedocs.io/en/develop/index.html
#+ 7. UCSC faCount - https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
#+ 8. USCS bedGraphToBigWig - https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
#+ 
#########################################################################################
# END_HEADER
#########################################################################################

# Help output
#== needed variables ==#
SCRIPT_HEADSIZE=$(grep -n "^#\ END_HEADER" "${0}" | cut -d: -f1)
SCRIPT_NAME=$(basename "${0}")

#== usage functions ==#
usagefull() { head -"${SCRIPT_HEADSIZE}" "${0}"| grep -e "^#%" | sed -e "s/^#[%+]\ //g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g" ; }
scriptinfo() { head -"${SCRIPT_HEADSIZE}" "${0}" | grep -e "^#+" | sed -e "s/^#+\ //g" -e "s/\${SCRIPT_NAME}/${SCRIPT_NAME}/g"; }

if [ "${1}" = "-h" ] || [ "${1}" = "--help" ] ; then
    usagefull && scriptinfo
    exit 0
fi


# Testing for input
if [ -z "${1}" ]; then
  echo -e "\nInput directory not supplied. Please supply a directory with raw FASTQ files."
  exit 0
elif [ -z "${2}" ]; then
  echo -e "\nPlease supply a second numerical argument to identify file name handles [3/4/5]"
  exit 0
elif [ -z "${3}" ]; then
  echo -e "\nPlease supply a third argument to run preprocessing on raw FASTQ files to create BigWig files [y/yes for Yes or n/no for No]."
  exit 0
elif [ -z "${4}" ]; then
  echo -e "\nPlease supply a fourth argument to run adapter trimming and alignment on raw FASTQ files [y/yes for Yes or n/no for No]."
  exit 0
elif [ -z "${5}" ]; then
  echo -e "\nPlease supply a fifth argument to run run spike-in calibration and PCR deduplication of alignment BAM files [y/yes for Yes or n/no for No]."
  exit 0
elif [ -z "${6}" ]; then
  echo -e "\nPlease supply a sixth argument to run peakcalling on all calibrated BAM files via Macs2 [y/yes for Yes or n/no for No]."
  exit 0
elif [ -z "${7}" ]; then
  echo -e "\nPlease supply a seventh argument to process GFP files without their control IGG file [y/yes for Yes or n/no for No]."
  exit 0
else
  echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Running Cut&Run pipeline with the following arguments -\nInput folder:\t$(realpath ${1})\nFile ID:\t${2}\nRun QC?\t${3}\nPreprocess FASTQ?\t${4}\nProcess BAM?\t${5}\nRun Peakcalling?\t${6}\nSkip GFP processing if IGG not present?\t${7}"
fi

# Changing working directory to input directory
PROJECT_PATH=$(realpath "${1}")
echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Input directory: ${PROJECT_PATH}"  # print the directory name which is being processed
# cd into each directory in the directory
cd "${PROJECT_PATH}" || { echo "cd into input directory failed! Please check your working directory." ; exit 1 ; }


if [ "${3}" = "y" ] || [ "${3}" == "Y" ] || [ "${3}" == "yes" ] || [ "${3}" == "Yes" ]; then
  echo -e "\n############"
  echo "# FASTQ QC #"
  echo "############"

  for INPUTFILE in $(find . -type f -name "*.fastq*" -exec realpath {} +)
  do
    SAMPLEID=$(basename "$INPUTFILE" | grep -o "^.*_R[12]")
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Processing ${SAMPLEID}"

    # Data Pre-processing - Quality Control using FastQC
    [ ! -d "$SAMPLEID"_raw_qc ] && mkdir -p "$SAMPLEID"_raw_qc
    fastqc -t 12 --nogroup -o "$SAMPLEID"_raw_qc "$INPUTFILE"
  done

  # Run MultiQC report on FastQC output
  multiqc --no-data-dir -f ./ > "${PROJECT_PATH}"/multiqc.log 2>&1

  # Clean up FastQC outputs
  mkdir -p fastqc
  mv -t fastqc/ ./*_raw_qc/
fi

echo -e "\n#########################################################"
echo "# Creating necessary folders and define input constants #"
echo "#########################################################"


# Uncomment to use Assembly 22, and comment out Assembly 21 lines below
# ASSEMBLY="22"
# CA_GENOME=$(find ~ -type f -name "ca22\.genome" -exec realpath {} +)
# [[ ! -f "${CA_GENOME}" ]] && echo -e "\nC. albicans genome file not found. Please execute this command - grep \"^>\" path/to/C_albicans_SC5314_A22_current_chromosomes.fasta | sed 's/_C_albicans_SC5314\ (/\ /g' | sed 's/>//g' | cut -d' ' -f1,2 | column -t > ~/ca22.genome\n"
# CA_REF="/home/aparopkari/ca22_bt2/ca22"
# Effective genome size = total number of base pairs minus the total number of ‘N’
# CA_GENOME_SIZE=28588771  # calculated via len - N from /home/aparopkari/ca22_genome/C_albicans_SC5314_A22_current_chromosomes_tidy_effective_genome_summary.txt

ASSEMBLY="21"
CA_GENOME=$(find ~ -type f -name "ca21\.genome" -exec realpath {} +)
CA_REF="/home/aparopkari/ca21_genome/ca21_bt2/ca21"
[[ ! -f "${CA_GENOME}" ]] && echo -e "\nC. albicans genome file not found. Please execute this command - grep "^>" /path/to/C_albicans_SC5314_A21_chromosomes.fasta | sed 's/_C_albicans_SC5314\ (/\ /g' | sed 's/>//g' | cut -d' ' -f1,2 | column -t > ~/ca21.genome\n"
# Effective genome size = total number of base pairs minus the total number of ‘N’
CA_GENOME_SIZE=14280188 # calculated via awk 'BEGIN{OFS="\t"} NR == 2 {print $2 - $7}' C_albicans_SC5314_A21_chromosomes_cleaned_nochrMT_counts.txt

# Dont edit this next block of code
EC_REF="/home/aparopkari/ecoli_genome/ecoli_bt2/ecoli_k12"
FIX_BDG=$(find ~ -type f -name "fix_bedgraph_coordinates.py" -exec realpath {} +)
SUBTRACT_IGG_SIGNAL=$(find ~ -type f -name "subtract_control_bed.py" -exec realpath {} +)
BDG_TO_BW=$(find ~ -type f -name "bedGraphToBigWig" -exec realpath {} +)
[[ ! -f "${BDG_TO_BW}" ]] && echo -e "\nUCSC bedGraphToBigWig script not found. Please execute this command - wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && chmod 755 bedGraphToBigWig - to download and make the script executable before re-running Cut&Run analysis pipeline.\n"
THREADS=12


# Create necessary folders

mkdir -p "${PROJECT_PATH}"/trimmed_reads/
[ -d "${PROJECT_PATH}"/trimmed_reads ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath trimmed_reads) directory successfully created."

mkdir -p "${PROJECT_PATH}"/alignment/bam/bowtie2_summary alignment/bam
[ -d "${PROJECT_PATH}"/alignment/bam/bowtie_summary ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath alignment/*) directory successfully created."

mkdir -p "${PROJECT_PATH}"/alignment/removeDuplicate/picard_summary
[ -d "${PROJECT_PATH}"/alignment/removeDuplicate/picard_summary ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath alignment/removeDuplicate) directory successfully created."

mkdir -p "${PROJECT_PATH}"/alignment/bam/fragment_len
[ -d "${PROJECT_PATH}"/alignment/bam/fragment_len ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath alignment/bam/fragment_len) directory successfully created."

mkdir -p "${PROJECT_PATH}"/peakcalling/macs2
[ -d "${PROJECT_PATH}"/peakcalling/macs2 ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath peakcalling/macs2) directory successfully created."

mkdir -p "${PROJECT_PATH}"/alignment/bigwig
[ -d "${PROJECT_PATH}"/alignment/bigwig ] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: $(realpath alignment/bigwig) directory successfully created."


if [ "${4}" = "y" ] || [ "${4}" == "Y" ] || [ "${4}" == "yes" ] || [ "${4}" == "Yes" ]; then
  for R1 in $(find ./ -type f -name "*R1_001.fastq.gz" -exec realpath {} +)
  do
    CDS=$(basename "$R1" .fastq.gz | cut -d "_" -f$(seq -s, "${2}"))
    echo -e "\n########################################################################################"
    echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: PROCESSING $CDS sample"
    echo "########################################################################################"

    R2="${R1//_R1/_R2}"
    R1_PAIRED="${PROJECT_PATH}/trimmed_reads/$(basename ${R1} .fastq.gz)_paired_trimmed.fastq.gz"
    R2_PAIRED="${PROJECT_PATH}/trimmed_reads/$(basename ${R2} .fastq.gz)_paired_trimmed.fastq.gz"

    echo -e "\n####################"
    echo "# ADAPTER TRIMMING #"
    echo "####################"

    CMD1="cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 0 -m 30 --report minimal -o ${R1_PAIRED} -p ${R2_PAIRED} ${R1} ${R2}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Adapter Trimming\n${CMD1}"
    $CMD1


    echo -e "\n#######################################"
    echo "# Aligning to Candida albicans genome #"
    echo "#######################################"

    BAM_OUT="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2.bam"
    BAM_LOG="${PROJECT_PATH}/alignment/bam/bowtie2_summary/${CDS}_bowtie2_alignment.log"

    # Align to Candida albicans A22 genome
    bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 --dovetail -p "${THREADS}" -x "${CA_REF}" -1 "${R1_PAIRED}" -2 "${R2_PAIRED}" 2> "${BAM_LOG}" | samtools view -@ "${THREADS}" -bh - > "${BAM_OUT}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Aligned to C. albicans ${ASSEMBLY} genome and saved $(samtools view -@ ${THREADS} ${BAM_OUT} | wc -l) aligned reads file to ${BAM_OUT}"


    echo -e "\n#####################"
    echo "# Sorting BAM file #"
    echo "####################"

    # Sort by coordinate
    SORTED_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted.bam"
    samtools sort -l 9 -O BAM -@ "${THREADS}" "${BAM_OUT}" -o "${SORTED_BAM}" || exit 0
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Sorted BAM file by coordinates and saved to ${SORTED_BAM}"
    BOWTIE2_SORTED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${SORTED_BAM} | wc -l) / 2" | bc)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Number of sorted alignments: ${BOWTIE2_SORTED_ALN}"

    echo -e "\n#####################################################"
    echo "# Remove PCR duplicates from C. albicans alignments #"
    echo "#####################################################"

    # Mark duplicates
    DUP_BAM="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_bowtie2_sorted_dup_marked.bam"
    DUP_BAM_LOG="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_bowtie2_sorted_dup_marked.log"
    DUP_METRIC="${PROJECT_PATH}/alignment/removeDuplicate/picard_summary/$(basename ${DUP_BAM} .bam).txt"
    picard MarkDuplicates I="${SORTED_BAM}" O="${DUP_BAM}" METRICS_FILE="${DUP_METRIC}" VALIDATION_STRINGENCY=SILENT 2> "${DUP_BAM_LOG}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Marked PCR duplicates and saved marked BAM file to ${DUP_BAM}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Duplications metrics saved to ${DUP_METRIC}"

    # Remove duplicates
    DEDUP_BAM="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_bowtie2_sorted_dedup.bam"
    DEDUP_BAM_LOG="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_dedup.log"
    DEDUP_METRIC="${PROJECT_PATH}/alignment/removeDuplicate/picard_summary/$(basename ${DEDUP_BAM} .bam).txt"
    picard MarkDuplicates I="${DUP_BAM}" O="${DEDUP_BAM}" REMOVE_DUPLICATES=true METRICS_FILE="${DEDUP_METRIC}" VALIDATION_STRINGENCY=SILENT 2> "${DEDUP_BAM_LOG}"
    ALN_DEDUP_PCT=$(echo "($(samtools view -@ "${THREADS}" "${DEDUP_BAM}" | wc -l) / $(samtools view -@ "${THREADS}" "${BAM_OUT}" | wc -l) * 100)" | bc -l | xargs printf "%'.3f\n")
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Removed PCR duplicates from BAM file and saved BAM file WITHOUT PCR duplicates to ${DEDUP_BAM}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Deduplications metrics saved to ${DEDUP_METRIC}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Number of alignments after removing PCR duplicates: $(samtools view -@ "${THREADS}" ${DEDUP_BAM} | wc -l)"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: PCR duplicate removal retained ${ALN_DEDUP_PCT}% of all C. albicans alignments"


    echo -e "\n#######################################"
    echo "# Aligning to Escherichia coli genome #"
    echo "#######################################"

    # Align to Escherichia coli K12 genome
    SpikeIn_BAM_OUT="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_spikein.bam"
    SpikeIn_BAM_LOG="${PROJECT_PATH}/alignment/bam/bowtie2_summary/$(basename ${R1} .fastq.gz | cut -d_ -f1,2,3)_bowtie2_spikeIn_alignment.log"

    bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-dovetail -p "${THREADS}" -x "${EC_REF}" -1 "${R1_PAIRED}" -2 "${R2_PAIRED}" 2> "${SpikeIn_BAM_LOG}" | samtools view -@ "${THREADS}" -bh - > "${SpikeIn_BAM_OUT}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Aligned to E. coli K12 genome and saved $(samtools view -@ ${THREADS} ${SpikeIn_BAM_OUT} | wc -l) aligned reads file to ${SpikeIn_BAM_OUT}"


    echo -e "\n#####################"
    echo "# Sorting BAM file #"
    echo "####################"

    # Sort by coordinate
    SpikeIn_SORTED_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_spikein_sorted.bam"
    samtools sort -l 9 -O BAM -@ "${THREADS}" "${SpikeIn_BAM_OUT}" -o "${SpikeIn_SORTED_BAM}" || exit 0
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Sorted BAM file by coordinates and saved to ${SpikeIn_SORTED_BAM}"
    BOWTIE2_SPIKEIN_SORTED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${SpikeIn_SORTED_BAM} | wc -l) / 2" | bc)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Number of sorted E. coli alignments: ${BOWTIE2_SPIKEIN_SORTED_ALN}"

    echo -e "\n#################################################"
    echo "# Remove PCR duplicates from E. coli alignments #"
    echo "#################################################"

    # Mark duplicates
    SPIKEIN_DUP_BAM="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_bowtie2_spikein_sorted_dup_marked.bam"
    SPIKEIN_DUP_BAM_LOG="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_spikein_dup_marked.log"
    SPIKEIN_DUP_METRIC="${PROJECT_PATH}/alignment/removeDuplicate/picard_summary/$(basename ${DUP_BAM} .bam)_spikein.txt"
    picard MarkDuplicates I="${SpikeIn_SORTED_BAM}" O="${SPIKEIN_DUP_BAM}" METRICS_FILE="${SPIKEIN_DUP_METRIC}" VALIDATION_STRINGENCY=SILENT 2> "${SPIKEIN_DUP_BAM_LOG}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Marked PCR duplicates and saved marked BAM file to ${SPIKEIN_DUP_BAM}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Duplications metrics saved to ${SPIKEIN_DUP_METRIC}"

    # Remove duplicates
    SPIKEIN_DEDUP_BAM="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_bowtie2_spikein_sorted_dedup.bam"
    SPIKEIN_DEDUP_BAM_LOG="${PROJECT_PATH}/alignment/removeDuplicate/${CDS}_spikein_dedup.log"
    SPIKEIN_DEDUP_METRIC="${PROJECT_PATH}/alignment/removeDuplicate/picard_summary/$(basename ${SPIKEIN_DEDUP_BAM} .bam)_spikein.txt"
    picard MarkDuplicates I="${SPIKEIN_DUP_BAM}" O="${SPIKEIN_DEDUP_BAM}" REMOVE_DUPLICATES=true METRICS_FILE="${SPIKEIN_DEDUP_METRIC}" VALIDATION_STRINGENCY=SILENT 2> "${SPIKEIN_DEDUP_BAM_LOG}"
    SPIKEIN_ALN_DEDUP_PCT=$(echo "($(samtools view -@ "${THREADS}" "${SPIKEIN_DEDUP_BAM}" | wc -l) / $(samtools view -@ ${THREADS} "${SpikeIn_BAM_OUT}" | wc -l) * 100)" | bc -l | xargs printf "%'.3f\n")
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Removed PCR duplicates from BAM file and saved BAM file WITHOUT PCR duplicates to ${SPIKEIN_DEDUP_BAM}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Deduplications metrics saved to ${SPIKEIN_DEDUP_METRIC}"
    BOWTIE2_SPIKEIN_SORTED_DEDUP_ALN=$(echo "$(samtools view -@ "${THREADS}" ${SPIKEIN_DEDUP_BAM} | wc -l) / 2" | bc)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Number of alignments after removing PCR duplicates: ${BOWTIE2_SPIKEIN_SORTED_DEDUP_ALN}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: PCR duplicate removal retained ${SPIKEIN_ALN_DEDUP_PCT}% of all E. coli alignments"
  done
fi


if  [ "${5}" = "y" ] || [ "${5}" == "Y" ] || [ "${5}" == "yes" ] || [ "${5}" == "Yes" ]; then
  for DEDUP_BAM in $(find ${PROJECT_PATH}/alignment/removeDuplicate -type f -name "*_bowtie2_sorted_dedup.bam" -exec realpath {} +)
  do
    [[ $(basename "${DEDUP_BAM}" .bam | cut -d"_" -f1) = "SN250" ]] && continue
    CDS=$(basename "${DEDUP_BAM}" .bam | cut -d "_" -f$(seq -s, "${2}"))
    SPIKEIN_DEDUP=$(find ${PROJECT_PATH}/alignment/removeDuplicate -type f -name "${CDS}_bowtie2_spikein_sorted_dedup.bam" -exec realpath {} +)

    echo -e "\n########################################################################################"
    echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: PROCESSING ${CDS} sample"
    echo "########################################################################################"

    echo -e "\n########################"
    echo "# Spike-in Calibration #"
    echo "########################"

    CALIBRATED_DEDUP_ALN="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dedup_spikein_normalized.bam"
    CALIBRATED_DEDUP_ALN_LOG="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dedup_spikein_normalized.log"

    if [[ "${DEDUP_BAM}" =~ .*GFP.* ]]; then
      GFP_reads=$(echo "$(samtools view -@ ${THREADS} ${SPIKEIN_DEDUP} | wc -l) / 2" | bc)

      # if IGG file not found, skip processing this GFP file and continue processing next file GFP file
      if [ "${7}" = "y" ] || [ "${7}" == "Y" ] || [ "${7}" == "yes" ] || [ "${7}" == "Yes" ]; then
        IGG_A="${SPIKEIN_DEDUP//GFP_3/IGG_1}"
        IGG_A_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_A} | wc -l) / 2" | bc)
        IGG_B="${SPIKEIN_DEDUP//GFP_3/IGG_2}"
        IGG_B_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_B} | wc -l) / 2" | bc)
        [[ "${IGG_A_reads}" -gt "${IGG_B_reads}" ]] && IGG_reads="${IGG_B_reads}" || IGG_reads="${IGG_A_reads}"
      else
        IGG="${SPIKEIN_DEDUP//GFP/IGG}"
        [[ ! -f "${IGG}" ]] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG} file DOES NOT exist, skipping E.coli spike-in calibration for ${SPIKEIN_DEDUP}" && continue || IGG_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG} | wc -l) / 2" | bc)
      fi

      if [[ "${GFP_reads}" -gt "${IGG_reads}" ]]; then
        SubSampleFactor=$(echo "1 + (${IGG_reads} / ${GFP_reads})" | bc -l)
        echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP_reads} (GFP E. coli spike-in reads) > ${IGG_reads} (IGG E. coli spike-in reads). Subsampling factor for ${CDS} (without PCR duplicates) = ${SubSampleFactor}"
        samtools view -@ "${THREADS}" -hbs "${SubSampleFactor}" ${DEDUP_BAM} > ${CALIBRATED_DEDUP_ALN} 2> ${CALIBRATED_DEDUP_ALN_LOG}
        BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DEDUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITHOUT PCR duplicates to ${CALIBRATED_DEDUP_ALN}"
      else
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG_reads} (IGG E. coli spike-in reads) > ${GFP_reads} (GFP E. coli spike-in reads)\n$(basename ${DEDUP_BAM}) will be saved to ${CALIBRATED_DEDUP_ALN}\n"
        cp "${DEDUP_BAM}" "${CALIBRATED_DEDUP_ALN}"
        BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DEDUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITHOUT PCR duplicates to ${CALIBRATED_DEDUP_ALN}"
      fi

    else
      IGG_reads=$(echo "$(samtools view -@ ${THREADS} ${SPIKEIN_DEDUP} | wc -l) / 2" | bc)

      # if GFP file not found, skip processing this IGG file and continue processing next file IGG file
      if [ "${7}" = "y" ] || [ "${7}" == "Y" ] || [ "${7}" == "yes" ] || [ "${7}" == "Yes" ]; then
        GFP_A="${SPIKEIN_DEDUP//IGG_3/GFP_1}"
        GFP_A_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP_A} | wc -l) / 2" | bc)
        GFP_B="${SPIKEIN_DEDUP//IGG_3/GFP_2}"
        GFP_B_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP_B} | wc -l) / 2" | bc)
        [[ "${GFP_A_reads}" -gt "${GFP_B_reads}" ]] && GFP_reads="${GFP_B_reads}" || GFP_reads="${GFP_A_reads}"
      else
        GFP="${SPIKEIN_DEDUP//IGG/GFP}"
        [[ ! -f "${GFP}" ]] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP} file DOES NOT exist, skipping E. coli spike-in calibration for ${SPIKEIN_DEDUP}" && continue || GFP_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP} | wc -l) / 2" | bc)
      fi

      if [[ "${IGG_reads}" -gt "${GFP_reads}" ]]; then
        SubSampleFactor=$(echo "1 + (${GFP_reads} / ${IGG_reads})" | bc -l)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG_reads} (IGG E. coli spike-in reads) > ${GFP_reads} (GFP E. coli spike-in reads). Subsampling factor for ${CDS} (without PCR duplicates) = ${SubSampleFactor}"
        samtools view -@ "${THREADS}" -hbs "${SubSampleFactor}" "${DEDUP_BAM}" > ${CALIBRATED_DEDUP_ALN} 2> ${CALIBRATED_DEDUP_ALN_LOG}
        BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DEDUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITHOUT PCR duplicates to ${CALIBRATED_DEDUP_ALN}"
      else
        echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP_reads} (GFP E. coli spike-in reads) > ${IGG_reads} (IGG E. coli spike-in reads)\n$(basename ${DEDUP_BAM}) will be saved to ${CALIBRATED_DEDUP_ALN}\n"
        cp "${DEDUP_BAM}" "${CALIBRATED_DEDUP_ALN}"
        BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DEDUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITHOUT PCR duplicates to ${CALIBRATED_DEDUP_ALN}"
      fi
    fi


    DUP_BAM=$(find "${PROJECT_PATH}"/alignment/bam -type f -name "${CDS}_bowtie2_sorted.bam")
    SPIKEIN_DUP=$(find "${PROJECT_PATH}"/alignment/bam -type f -name "${CDS}_bowtie2_spikein_sorted.bam")
    CALIBRATED_DUP_ALN="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dup_spikein_normalized.bam"
    CALIBRATED_DUP_ALN_LOG="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dup_spikein_normalized.log"

    if [[ "${DUP_BAM}" =~ .*GFP.* ]]; then
      GFP_reads=$(echo "$(samtools view -@ ${THREADS} ${SPIKEIN_DUP} | wc -l) / 2" | bc)

      # if IGG file not found, skip processing this GFP file and continue processing next file GFP file
      if [ "${7}" = "y" ] || [ "${7}" == "Y" ] || [ "${7}" == "yes" ] || [ "${7}" == "Yes" ]; then
        IGG_A="${SPIKEIN_DUP//GFP_3/IGG_1}"
        IGG_A_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_A} | wc -l) / 2" | bc)
        IGG_B="${SPIKEIN_DUP//GFP_3/IGG_2}"
        IGG_B_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_B} | wc -l) / 2" | bc)
        [[ "${IGG_A_reads}" -gt "${IGG_B_reads}" ]] && IGG_reads="${IGG_B_reads}" || IGG_reads="${IGG_A_reads}"
      else
        IGG="${SPIKEIN_DUP//GFP/IGG}"
        [[ ! -f "${IGG}" ]] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG} file DOES NOT exist, skipping E.coli spike-in calibration for ${SPIKEIN_DUP}" && continue || IGG_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG} | wc -l) / 2" | bc)
      fi

      if [[ "${GFP_reads}" -gt "${IGG_reads}" ]]; then
        SubSampleFactor=$(echo "1 + (${IGG_reads} / ${GFP_reads})" | bc -l)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP_reads} (GFP E. coli spike-in reads) > ${IGG_reads} (IGG E. coli spike-in reads). Subsampling factor for ${CDS} (with PCR duplicates) = ${SubSampleFactor}"
        samtools view -@ "${THREADS}" -hbs "${SubSampleFactor}" ${DUP_BAM} > ${CALIBRATED_DUP_ALN} 2> ${CALIBRATED_DUP_ALN_LOG}
        BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITH PCR duplicates to ${CALIBRATED_DUP_ALN}"
      else
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG_reads} (IGG E. coli spike-in reads) > ${GFP_reads} (GFP E. coli spike-in reads)\n$(basename ${DUP_BAM}) will be saved to ${CALIBRATED_DUP_ALN}\n"
        cp "${DUP_BAM}" "${CALIBRATED_DUP_ALN}"
        BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITH PCR duplicates to ${CALIBRATED_DUP_ALN}"
      fi

    else
      IGG_reads=$(echo "$(samtools view -@ ${THREADS} ${SPIKEIN_DUP} | wc -l) / 2" | bc)

      # if GFP file not found, skip processing this IGG file and continue processing next file IGG file
      if [ "${7}" = "y" ] || [ "${7}" == "Y" ] || [ "${7}" == "yes" ] || [ "${7}" == "Yes" ]; then
        GFP_A="${SPIKEIN_DUP//IGG_3/GFP_1}"
        GFP_A_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP_A} | wc -l) / 2" | bc)
        GFP_B="${SPIKEIN_DUP//IGG_3/GFP_2}"
        GFP_B_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP_B} | wc -l) / 2" | bc)
        [[ "${GFP_A_reads}" -gt "${GFP_B_reads}" ]] && GFP_reads="${GFP_B_reads}" || GFP_reads="${GFP_A_reads}"
      else
        GFP="${SPIKEIN_DUP//IGG/GFP}"
        [[ ! -f "${GFP}" ]] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP} file does not exist, skipping E. coli spike-in calibration for ${SPIKEIN_DUP}" && continue || GFP_reads=$(echo "$(samtools view -@ ${THREADS} ${GFP} | wc -l) / 2" | bc)
      fi

      if [[ "${IGG_reads}" -gt "${GFP_reads}" ]]; then
        SubSampleFactor=$(echo "1 + (${GFP_reads} / ${IGG_reads})" | bc -l)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG_reads} (IGG E. coli spike-in reads) > ${GFP_reads} (GFP E. coli spike-in reads). Subsampling factor for ${CDS} (with PCR duplicates) = ${SubSampleFactor}"
        samtools view -@ "${THREADS}" -hbs "${SubSampleFactor}" "${DUP_BAM}" > ${CALIBRATED_DUP_ALN} 2> ${CALIBRATED_DUP_ALN_LOG}
        BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITH PCR duplicates to ${CALIBRATED_DUP_ALN}"
      else
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${GFP_reads} (GFP E. coli spike-in reads) > ${IGG_reads} (IGG E. coli spike-in reads)\n$(basename ${DUP_BAM}) will be saved to ${CALIBRATED_DUP_ALN}\n"
        cp "${DUP_BAM}" "${CALIBRATED_DUP_ALN}"
        BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${CALIBRATED_DUP_ALN} | wc -l) / 2" | bc)
        echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved ${BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_ALN} mapped spike-in calibrated alignments WITH PCR duplicates to ${CALIBRATED_DUP_ALN}"
      fi
    fi

    echo -e "\n###########################################################"
    echo "# Size filtering BAM file with and without PCR duplicates #"
    echo "###########################################################"

    # Filter sequences >= 120 bp without PCR duplicates
    SIZE_FILTERED_DEDUP_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dedup_spikein_normalized_size_filtered.bam"
    samtools view -h@ "${THREADS}" "${CALIBRATED_DEDUP_ALN}" | awk 'BEGIN{OFS="\t"; FS=OFS} {if (/^@/ || $9*$9 <= 14400) print $0}' > "${SIZE_FILTERED_DEDUP_BAM}"
    NUMERATOR=$(echo "$(samtools view -@ ${THREADS} ${SIZE_FILTERED_DEDUP_BAM} | wc -l) / 2" | bc)
    DENOMINATOR=$(echo "$(samtools view -@ ${THREADS} ${CALIBRATED_DEDUP_ALN} | wc -l) / 2" | bc)
    DEDUP_ALN_FIL_PCT=$(echo "scale=3; ${NUMERATOR} / ${DENOMINATOR} * 100" | bc -l)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Retained ${NUMERATOR} (${DEDUP_ALN_FIL_PCT}%) 120bp alignments WITHOUT PCR duplicates and saved them to ${SIZE_FILTERED_DEDUP_BAM}"

    # Filter sequences >= 120 bp with PCR duplicates
    SIZE_FILTERED_DUP_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dup_spikein_normalized_size_filtered.bam"
    samtools view -h@ "${THREADS}" "${CALIBRATED_DUP_ALN}" | awk 'BEGIN{OFS="\t"; FS=OFS} {if (/^@/ || $9*$9 <= 14400) print $0}' > "${SIZE_FILTERED_DUP_BAM}"
    NUMERATOR=$(echo "$(samtools view -@ ${THREADS} ${SIZE_FILTERED_DUP_BAM} | wc -l) / 2" | bc)
    DENOMINATOR=$(echo "$(samtools view -@ ${THREADS} ${CALIBRATED_DUP_ALN} | wc -l) / 2" | bc)
    DUP_ALN_FIL_PCT=$(echo "scale=3; ${NUMERATOR} / ${DENOMINATOR} * 100" | bc -l)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Retained ${NUMERATOR} (${DUP_ALN_FIL_PCT}%) 120bp alignments WITH PCR duplicates and saved them to ${SIZE_FILTERED_DUP_BAM}"


    echo -e "\n############################################"
    echo "# Assess mapped fragment size distribution #"
    echo "############################################"

    # Extract the 9th column from the alignment bam file which is the fragment length
    DEDUP_FRAG_LEN="${PROJECT_PATH}/alignment/bam/fragment_len/${CDS}_bowtie2_sorted_dedup_spikein_normalized_size_filtered_frag_len.txt"
    samtools view -@ "${THREADS}" -h -F 0x04 "${SIZE_FILTERED_DEDUP_BAM}" | awk -F"\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > "${DEDUP_FRAG_LEN}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved mapped fragment size distribution of alignments without PCR duplicates to ${DEDUP_FRAG_LEN}"

    DUP_FRAG_LEN="${PROJECT_PATH}/alignment/bam/fragment_len/${CDS}_bowtie2_sorted_dup_spikein_normalized_size_filtered_frag_len.txt"
    samtools view -@ "${THREADS}" -h -F 0x04 "${SIZE_FILTERED_DUP_BAM}" | awk -F"\t" 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > "${DUP_FRAG_LEN}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved mapped fragment size distribution of alignments with PCR duplicates to ${DUP_FRAG_LEN}"


    echo -e "\n###########################################"
    echo "# File format conversion for Peak Calling #"
    echo "###########################################"

    # Filter and keep the mapped read pairs
    SORTED_DEDUP_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dedup_spikein_normalized_size_filtered_sorted.bam"
	SORTED_DUP_BAM="${PROJECT_PATH}/alignment/bam/${CDS}_bowtie2_sorted_dup_spikein_normalized_size_filtered_sorted.bam"
    samtools sort -l 9 --write-index -O BAM -@ "${THREADS}" "${SIZE_FILTERED_DEDUP_BAM}" -o "${SORTED_DEDUP_BAM}" || exit 0
    samtools index "${SORTED_DEDUP_BAM}"
    samtools sort -l 9 --write-index -O BAM -@ "${THREADS}" "${SIZE_FILTERED_DUP_BAM}" -o "${SORTED_DUP_BAM}" || exit 0
    samtools index "${SORTED_DUP_BAM}"
    BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_SORTED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${SORTED_DEDUP_BAM} | wc -l) / 2" | bc)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved Macs2 input BAM file with ${BOWTIE2_SORTED_DEDUP_SPIKEIN_NORMALIZED_SORTED_ALN} alignments WITHOUT PCR duplicates to ${SORTED_DEDUP_BAM}"
    BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_SORTED_ALN=$(echo "$(samtools view -@ "${THREADS}" ${SORTED_DUP_BAM} | wc -l) / 2" | bc)
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved Macs2 input BAM file with ${BOWTIE2_SORTED_DUP_SPIKEIN_NORMALIZED_SORTED_ALN} alignments WITH PCR duplicates to ${SORTED_DUP_BAM}"
  done
fi


if [ "${6}" = "y" ] || [ "${6}" == "Y" ] || [ "${6}" == "yes" ] || [ "${6}" == "Yes" ]; then
    echo -e "\n#################"
    echo "# Calling peaks #"
    echo "#################"

  for INPUT_BAM in $(find "${PROJECT_PATH}"/alignment/bam -type f -name "*_spikein_normalized_size_filtered_sorted.bam" -exec realpath {} +)
  do
    CDS=$(basename "${INPUT_BAM}" .bam | cut -d "_" -f$(seq -s, "${2}"))
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Calling peaks for ${CDS}"
    PEAK_CALLED="${PROJECT_PATH}/peakcalling/macs2/$(basename "${INPUT_BAM}" .bam)_macs2_fdr"
    PEAK_CALLED_LOG="${PROJECT_PATH}/peakcalling/macs2/$(basename "${INPUT_BAM}" .bam)_macs2_fdr.log"
    macs2 callpeak -t "${INPUT_BAM}" -g "${CA_GENOME_SIZE}" -f BAMPE -n "${PEAK_CALLED}" --outdir "$(realpath ${PROJECT_PATH}/peakcalling/macs2)" -q 0.01 -B --nolambda --nomodel --keep-dup all 2> "${PEAK_CALLED_LOG}"
    grep -v -e 'chrM' "${PEAK_CALLED}"_peaks.narrowPeak > "${PEAK_CALLED}"_nochrM_peaks.narrowPeak
    echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved $(wc -l ${PEAK_CALLED}_nochrM_peaks.narrowPeak | cut -d' ' -f1) Macs2 peaks with 1% FDR to ${PEAK_CALLED}_nochrM_peaks.narrowPeak"

    BDG_IGV="${PROJECT_PATH}/peakcalling/macs2/$(basename ${PEAK_CALLED}_treat_pileup.bdg .bdg)_for_IGV.bdg"
    BDG_IGV_LOG="${PROJECT_PATH}/peakcalling/macs2/$(basename ${PEAK_CALLED}_treat_pileup.bdg .bdg)_for_IGV.log"
    echo -e "\n$(ls ${PEAK_CALLED}_treat_pileup.bdg)\n${BDG_IGV}\n${BDG_IGV_LOG}\n"
    python "${FIX_BDG}" "${PEAK_CALLED}_treat_pileup.bdg" "${ASSEMBLY}" "${BDG_IGV}" > "${BDG_IGV_LOG}"
    sort -k1,1 -k2,2n -o "${BDG_IGV}" "${BDG_IGV}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Subtracted ${CDS} IGG signal from ${CDS} GFP BedGraph file and fixed the resulting BedGraph file"

    BW_IGV_OUT="${PROJECT_PATH}/alignment/bigwig/$(basename ${BDG_IGV} .bdg).bw"
    "${BDG_TO_BW}" "${BDG_IGV}" "${CA_GENOME}" "${BW_IGV_OUT}" || echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Error during BigWig file creation, troubleshoot this command: ${BDG_TO_BW} ${BDG_IGV} ${CA_GENOME} ${BW_IGV_OUT}" && echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved IGG subtracted and sorted BedGraph file as BigWig to ${BW_IGV_OUT}"
  done

  echo -e "\n#######################################################"
  echo "# Subtracting IGG from GFP and obtaining BigWig files #"
  echo "#######################################################"

  for GFP in $(find "${PROJECT_PATH}"/peakcalling/macs2/ -type f -name "*GFP*_macs2_fdr_treat_pileup.bdg" -exec realpath {} +)
  do
    CDS=$(basename "${GFP}" .bdg | cut -d "_" -f$(seq -s, "${2}"))
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Subtracting IGG signal from GFP signal for ${CDS}"
    GFP_PEAKS="${GFP//_treat_pileup.bdg/_nochrM_peaks.narrowPeak}"
    IGG_PEAKS="${GFP_PEAKS//GFP/IGG}"
    IGG=${GFP//GFP/IGG}
    if [[ ! -f "${IGG}" ]] && [ "${7}" = "y" ] || [ "${7}" == "Y" ] || [ "${7}" == "yes" ] || [ "${7}" == "Yes" ]; then
      [[ "${GFP}" =~ .*_dup_.* ]] && BAM=$(find "${PROJECT_PATH}"/alignment/bam -type f -name "${CDS}_bowtie2_spikein_sorted.bam") || BAM=$(find ${PROJECT_PATH}/alignment/removeDuplicate -type f -name "${CDS}_bowtie2_spikein_sorted_dedup.bam" -exec realpath {} +)
      IGG_A="${BAM//GFP_3/IGG_1}"
      IGG_A_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_A} | wc -l) / 2" | bc)
      IGG_B="${BAM//GFP_3/IGG_2}"
      IGG_B_reads=$(echo "$(samtools view -@ ${THREADS} ${IGG_B} | wc -l) / 2" | bc)
      [[ "${IGG_A_reads}" -gt "${IGG_B_reads}" ]] && IGG="${GFP//GFP_3/IGG_2}" || IGG="${GFP//GFP_3/IGG_1}"
      IGG_PEAKS="${IGG//_treat_pileup.bdg/_nochrM_peaks.narrowPeak}"
    else
      [[ ! -f "${IGG}" ]] && echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: ${IGG} file DOES NOT exist, skipping signal subtraction ${GFP}" && continue
    fi

    PEAK_COMPARED="${PROJECT_PATH}/peakcalling/macs2/$(basename ${GFP} .bdg)_final.bdg"
    PEAK_COMPARED_LOG="${PROJECT_PATH}/peakcalling/macs2/$(basename ${GFP} .bdg)_final.log"
    macs2 bdgcmp -t "${GFP}" -c "${IGG}" -m subtract -o "${PEAK_COMPARED}" --outdir "$(realpath ${PROJECT_PATH}/peakcalling/macs2)" 2> "${PEAK_COMPARED_LOG}"
    sort -k1,1 -k2,2n -o "${PEAK_COMPARED}" "${PEAK_COMPARED}"

    BDG_FIXED="${PROJECT_PATH}/peakcalling/macs2/$(basename ${PEAK_COMPARED} .bdg)_fixed.bdg"
    BDG_FIXED_LOG="${PROJECT_PATH}/peakcalling/macs2/$(basename ${PEAK_COMPARED} .bdg)_fixed.log"
    python "${FIX_BDG}" "${PEAK_COMPARED}" "${ASSEMBLY}" "${BDG_FIXED}" > "${BDG_FIXED_LOG}"
    sort -k1,1 -k2,2n -o "${BDG_FIXED}" "${BDG_FIXED}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Subtracted ${CDS} IGG signal from ${CDS} GFP BedGraph file and fixed the resulting BedGraph file"

    # Subtract IGG narrowPeak from GFP narrowPeak
    FINAL_PEAKS="${PROJECT_PATH}/peakcalling/macs2/$(basename ${GFP_PEAKS} .narrowPeak)_gfp_minus_igg.bed"
    FINAL_PEAKS_LOG="${PROJECT_PATH}/peakcalling/macs2/$(basename ${GFP_PEAKS} .narrowPeak)_gfp_minus_igg.log"
    subtractBed -a "${GFP_PEAKS}" -b "${IGG_PEAKS}" | mergeBed -c 4,5,6,7,8,9,10 -o distinct,sum,distinct,mean,mean,mean,sum > "${FINAL_PEAKS}" 2> "${FINAL_PEAKS_LOG}"
    echo -e "\n\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved $(wc -l ${FINAL_PEAKS} | cut -d' ' -f1) IGG subtracted and sorted peaks to ${FINAL_PEAKS}"

    BW_OUT="${PROJECT_PATH}/alignment/bigwig/$(basename ${BDG_FIXED} .bdg).bw"
    "${BDG_TO_BW}" "${BDG_FIXED}" "${CA_GENOME}" "${BW_OUT}" || echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Error during BigWig file creation, troubleshoot this command: ${BDG_TO_BW} ${GFP_PEAKS}_treat_pileup.bdg ${CA_GENOME} ${BW_OUT}" && echo -e "\033[1;4m$(date '+%A %B %d, %Y %T %Z')\033[0m: Saved IGG subtracted and sorted BedGraph file as BigWig to ${BW_OUT}"
  done
fi
