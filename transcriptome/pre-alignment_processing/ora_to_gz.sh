# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "  fastq raw file directory: ${fastqDir}"
echo "  ora reference file directory: ${oraRefDir}"
echo "  threads: ${threads}"

# [Main]

# Convert fastq.ora files to fastq.gz
# The .ora files are the compression format usd by Illumina
# This compression format results in a smaller file size than .gz
# However, this format is not as widely supported as .gz
# Therefore, we will convert the .ora files to .gz files

# The output will be stored in the same directory as the input files
for file in $(find ${fastqDir} -type f -name '*.fastq.ora'); do
    sample=$(basename $file .fastq.ora)
    echo "Converting ${sample}.fastq.ora to ${sample}.fastq.gz"
    orad --ora-reference ${oraRefDir} \
         --gz ${fastqDir}/${sample}.fastq.ora
    echo "Conversion of ${sample}.fastq.ora to ${sample}.fastq.gz complete"
done

# [END]