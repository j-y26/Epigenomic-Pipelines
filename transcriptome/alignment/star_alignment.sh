# [Config]
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_script.sh>"
    exit 1
fi

config_script=$1
source ${config_script}
echo "Running with config:"
echo "STARIndexDir: ${STARIndexDir}"
echo "STAROutputDir: ${STAROutputDir}"

# [Main]

# Check to ensure that the reference genome directory exists
if [ ! -d ${STARIndexDir} ]; then
    echo "Error: STAR index directory does not exist." >&2
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d ${STAROutputDir} ]; then
    mkdir ${STAROutputDir}
fi

# Align the reads to the reference genome
# Loop through the directory to match patterns ending in _R1_*.fastq.gz
# We use trimmed files for alignment, so we need to loop through the trimmed directory

for file in $(find ${trimmedDir} -type f -name '*_R1.fastq.gz'); do
    sample=$(basename $file _trimmed_R1.fastq.gz)
    echo "Aligning ${sample} to reference genome..."
    echo "This may take a while..."
    forward_file="${sample}_trimmed_R1.fastq.gz"
    reverse_file="${sample}_trimmed_R2.fastq.gz"
    STAR \
    --runThreadN ${threads} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${STAROutputDir}/${samplename} \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ${STARIndexDir} \
    --readFilesIn ${trimmedDir}/${forward_file} ${trimmedDir}/${reverse_file}
    echo "Alignment for ${sample} complete"
done

# [END]
