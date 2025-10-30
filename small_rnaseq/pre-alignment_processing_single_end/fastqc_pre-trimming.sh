# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  fastq raw file directory: ${fastqDir}"
echo "  FastQC output directory: ${fastqcOutDir}"
echo "  threads: ${threads}"

# [Main]

# Quality control of fastq.gz files with FastQC
# This script will run FastQC on all fastq.gz files in a directory and output 
# the results to a specified directory
# Note this code only supports paired-end reads

# Create output directory if it doesn't exist
if [ ! -d ${fastqcOutDir} ]; then
    mkdir ${fastqcOutDir}
    mkdir ${fastqcOutDir}/pre_trimming
fi

if [ ! -d ${fastqcOutDir}/pre_trimming ]; then
    mkdir ${fastqcOutDir}/pre_trimming
fi

# Run fastqc on all files ending in _R1_*.fastq.gz

echo "Running FastQC on all fastq.gz files in ${fastqDir}..."
echo "Analysis prior to trimming begins now."
fastqc -t ${threads} \
    -o ${fastqcOutDir}/pre_trimming \
    ${fastqDir}/*_R1_*.fastq.gz
echo "FastQC analysis complete."

# [END]