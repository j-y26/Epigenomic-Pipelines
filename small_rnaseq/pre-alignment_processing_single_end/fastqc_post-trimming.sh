# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Trimmed fastq file directory: ${trimmedDir}"
echo "  fastqcOutDir: ${fastqcOutDir}"
echo "  threads: ${threads}"

# [Main]

# Quality control of fastq.gz files with FastQC for post-trimming reads
# This script will run FastQC on all fastq.gz files in a directory and output 
# the results to a specified directory
# Note this code only supports paired-end reads

# Create output directory if it doesn't exist
if [ ! -d ${fastqcOutDir} ]; then
    mkdir ${fastqcOutDir}
    mkdir ${fastqcOutDir}/post_trimming
fi

if [ ! -d ${fastqcOutDir}/post_trimming ]; then
    mkdir ${fastqcOutDir}/post_trimming
fi

# Run fastqc on all files ending in R1_trimmed.fastq.gz

echo "Running FastQC on all fastq.gz files in ${trimmedDir}..."
echo "Analysis prior to trimming begins now."
fastqc -t ${threads} \
    -o ${fastqcOutDir}/post_trimming \
    ${trimmedDir}/*_R1_trimmed.fastq.gz
echo "FastQC analysis complete."

# [END]