# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  fastqDir: ${fastqDir}"
echo "  fastqcOutDir: ${fastqcOutDir}"
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

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and run FastQC on the paired files

for file in $(find ${fastqDir} -type f -name '*_R1_*.fastq.gz'); do
    sample=$(basename $file _${lane}_R1_001.fastq.gz)
    echo "Running FastQC on ${sample} before trimming"
    forward_file="${sample}_${lane}_R1_001.fastq.gz"
    reverse_file="${sample}_${lane}_R2_001.fastq.gz"
    fastqc -t ${threads} -o ${fastqcOutDir}/pre_trimming ${fastqDir}/${forward_file} ${fastqDir}/${reverse_file}
    echo "FastQC for ${sample} complete"
done

# [END]