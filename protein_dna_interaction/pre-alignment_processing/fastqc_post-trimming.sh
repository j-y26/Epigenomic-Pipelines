# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  Trimmed fastq directory: ${trimmedDir}"
echo "  FastQC output directory: ${fastqcOutDir}"
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

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and run FastQC on the paired files

for file in $(find ${trimmedDir} -type f -name '*_R1_*.fastq.gz'); do
    sample=$(basename $file _R1_trimmed.fastq.gz)
    echo "Running FastQC on ${sample} after trimming"
    forward_file="${sample}_R1_trimmed.fastq.gz"
    reverse_file="${sample}_R2_trimmed.fastq.gz"
    fastqc -t ${threads} -o ${fastqcOutDir}/post_trimming ${trimmedDir}/${forward_file} ${trimmedDir}/${reverse_file}
    echo "FastQC for ${sample} complete"
done

# [END]