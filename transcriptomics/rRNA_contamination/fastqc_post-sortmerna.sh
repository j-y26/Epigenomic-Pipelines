# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  rRNAcleanDir: ${rRNAcleanDir}"
echo "  fastqcOutDir: ${fastqcOutDir}"
echo "  threads: ${threads}"

# [Main]

# Quality control of fastq.gz files with FastQC for post-cleaning rRNA reads
# This script will run FastQC on all fastq.gz files in a directory and output 
# the results to a specified directory
# Note this code only supports paired-end reads

# Create output directory if it doesn't exist
if [ ! -d ${fastqcOutDir} ]; then
    mkdir ${fastqcOutDir}
    mkdir ${fastqcOutDir}/post_rRNA_cleaning
fi

if [ ! -d ${fastqcOutDir}/post_rRNA_cleaning ]; then
    mkdir ${fastqcOutDir}/post_rRNA_cleaning
fi

# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# and run FastQC on the paired files

for file in $(find ${trimmedDir} -type f -name '*R1_cleaned.fastq.gz'); do
    sample=$(basename $file _R1_cleaned.fastq.gz)
    echo "Running FastQC on ${sample} after cleaning rRNA reads"
    forward_file="${sample}_R1_cleaned.fastq.gz"
    reverse_file="${sample}_R2_cleaned.fastq.gz"
    fastqc -t ${threads} -o ${fastqcOutDir}/post_rRNA_cleaning ${rRNAcleanDir}/${forward_file} ${rRNAcleanDir}/${reverse_file}
    echo "FastQC for ${sample} complete"
done

# [END]