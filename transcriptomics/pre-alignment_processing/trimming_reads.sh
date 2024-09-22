# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  fastq raw file directory: ${fastqDir}"
echo "  Trimmed fastq file directory: ${trimmedDir}"
echo "  threads: ${threads}"
echo "  TRIMMOMATIC_JAR: ${TRIMMOMATIC_JAR}"
echo "  TRIMMOMATIC_ADAPTER_FILE: ${TRIMMOMATIC_ADAPTER_FILE}"

# [Main]

# Trimming of RNA-seq reads using Trimmomatic
# Since Trimmomatic is a java tool, it requires a compatible version of java to be installed

# Check that java is installed
java -version

# Create output directory if it doesn't exist
if [ ! -d ${trimmedDir} ]; then
    mkdir ${trimmedDir}
    mkdir ${trimmedDir}/unpaired
fi

if [ ! -d ${trimmedDir}/unpaired ]; then
    mkdir ${trimmedDir}/unpaired
fi

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and perform trimming using Trimmomatic

for file in $(find ${fastqDir} -type f -name '*_R1_*.fastq.gz'); do
    sample=$(basename $file _${lane}_R1_001.fastq.gz)
    echo "Trimming ${sample}"
    forward_file="${sample}_${lane}_R1_001.fastq.gz"
    reverse_file="${sample}_${lane}_R2_001.fastq.gz"
    java -jar ${TRIMMOMATIC_JAR} PE -threads ${threads} -phred33 \
        ${fastqDir}/${forward_file} ${fastqDir}/${reverse_file} \
        ${trimmedDir}/${sample}_trimmed_R1.fastq.gz ${trimmedDir}/unpaired/${sample}_unpaired_R1.fastq.gz \
        ${trimmedDir}/${sample}_trimmed_R2.fastq.gz ${trimmedDir}/unpaired/${sample}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC_ADAPTER_FILE}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:36
done

# Remove unpaired files
rm -r ${trimmedDir}/unpaired

# [END]