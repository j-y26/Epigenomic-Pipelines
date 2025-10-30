# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  fastq raw file directory: ${fastqDir}"
echo "  Trimmed fastq directory: ${trimmedDir}"
echo "  Adapter fasta file: ${adapterFastaFile}"
echo "  threads: ${threads}"

# [Main]

# Trimming of Small RNA-seq reads using fastp
# Small RNAseq are commonly single end reads, so this script is designed for # processing single end reads
# We also trim polyG and polyX tails, which could be artifacts from the sequencing process/library preparation

# Since miRNA smaller than 15 bases have been discovered in the literature,
# The minimum allowed read length is set to 10 bases

# Create output directory if it doesn't exist

if [ ! -d ${trimmedDir} ]; then
    mkdir ${trimmedDir}
    mkdir ${trimmedDir}/fastp_report
fi

if [ ! -d ${trimmedDir}/fastp_report ]; then
    mkdir ${trimmedDir}/fastp_report
fi

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and perform trimming using fastp

for file in $(find ${fastqDir} -type f -name '*_R1_*.fastq.gz'); do
    sample=$(basename $file _${lane}_R1_001.fastq.gz)
    echo "Trimming ${sample}"
    forward_file="${sample}_${lane}_R1_001.fastq.gz"
    fastp -i ${fastqDir}/${forward_file} \
    -o ${trimmedDir}/${sample}_R1_trimmed.fastq.gz \
    --trim_poly_g \
    --trim_poly_x \
    --length_required 10 \
    --adapter_fasta ${adapterFastaFile} \
    -w ${threads} \
    -h ${trimmedDir}/fastp_report/${sample}_fastp_report.html \
    -j ${trimmedDir}/fastp_report/${sample}_fastp_report.json 2>&1 | \
    tee ${trimmedDir}/fastp_report/${sample}_fastp_summary.txt

    echo "Done trimming ${sample}"
done

# [END]